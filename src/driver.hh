#ifndef DRIVER_HH
#define DRIVER_HH

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh> // added
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/pdelab/gridoperator/onestep.hh>
#include <string>
#include <stdexcept>

#include "bctype.hh"
#include "wavefem.hh"

template<typename GV>
void driver (const GV& gv, double dt, double c, double T_end, std::string name)
{
  double time = 0.0;

  const int dim = GV::dimension;
  using RF = double;                   // type for computations
  const int k = 2;
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV, double, double, k>;
  //using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, double, double, k>;
  FEM fem(gv);

  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS0 = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0>;
  GFS0 gfs0(gv,fem);

  using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
  using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
  using GFS = Dune::PDELab::PowerGridFunctionSpace<GFS0,2,VBE,OrderingTag>;
  GFS gfs(gfs0);

  using U_BCTypeParam = Dune::PDELab::PowerConstraintsParameters<BCTypeParam<GV>, 2>;

  using namespace Dune::Indices;
  gfs.child(_0).name("u0");
  gfs.child(_1).name("u1");

  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs);

  using BCE0 = BCExtension0<GV, double>;
  using BCE1 = BCExtension1<GV, double>;
  using BCE = Dune::PDELab::CompositeGridFunction<BCE0, BCE1>;

  BCE0 bce0(gv);
  BCE1 bce1(gv);
  BCE bce(bce0, bce1);
  bce.setTime(time);
  Dune::PDELab::interpolate(bce,gfs,z);

  BCTypeParam<GV> bctype0(gv);
  bctype0.setTime(time);
  U_BCTypeParam bctype(bctype0);
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  CC cc;

  Dune::PDELab::constraints(bctype,gfs,cc);

  std::cout << "constrained dofs=" << cc.size() << " of "
            << gfs.globalSize() << std::endl;
  //set_constrained_dofs(cc,0.0,z); // set zero Dirichlet boundary conditions




  using LOP = WaveFEM<U_BCTypeParam, BCTypeParam<GV>,FEM>;
  LOP lop(bctype,bctype0,c);
  using TLOP = WaveL2<FEM>;
  TLOP tlop;
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MBE mbe((int)pow(1+2*k,dim));
  using GO0 = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC>;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);
  using GO1 = Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,CC,CC>;
  GO1 go1(gfs,cc,gfs,cc,tlop,mbe);
  using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1>;
  IGO igo(go0,go1);

  using VTKW = Dune::SubsamplingVTKWriter<GV>;
  VTKW vtkwriter(gv,Dune::RefinementIntervals{1});
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, z);  // new
  Dune::VTKSequenceWriter<GV> writer(std::make_shared<VTKW>(vtkwriter), name);
  writer.write(time);

  // Linear problem solver
  using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR;
  LS ls(5000,false);
  using SLP = Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,Z>;
  SLP slp(igo,ls,1e-8);

  Dune::PDELab::OneStepThetaParameter<RF> method1(1.0);
  Dune::PDELab::Alexander2Parameter<RF> method2;

  Dune::PDELab::OneStepMethod<RF,IGO,SLP,Z,Z> osm(method2,igo,slp);
  osm.setVerbosityLevel(2);

//vremenska petlja
  while (time<T_end-1e-8)
    {

      bce.setTime(time+dt);
      bctype0.setTime(time+dt);
      using namespace Dune::Indices;
      bctype.child(_0).setTime(time+dt);
      bctype.child(_1).setTime(time+dt);
      cc.clear();
      Dune::PDELab::constraints(bctype,gfs,cc);
      Z znew(z);
      osm.apply(time,dt,z,znew);

      z = znew;
      time+=dt;

      writer.write(time);
    }
}
#endif

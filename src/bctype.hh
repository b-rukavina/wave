#ifndef BCTYPE_HH
#define BCTYPE_HH

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

//rubni i inicijalni uvjeti za u_0 = u
template<typename GV, typename RF>
class BCExtension0
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, BCExtension0<GV,RF> > {
  const GV& gv;
  double time;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;
  BCExtension0 (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    auto x = e.geometry().global(xlocal);

  /*  if (time == 0) {
      y = (x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5);
      y = std::max(0.0, 1.0-8.0*sqrt(y));
    }
    else y = 0;*/

  if (time==0) y=sin(x[0]);
  else if(x[0]<1e-5||x[0]>3*M_PI-1e-5) y = sin(-time);

  return;
  }

  inline const GV& getGridView () {return gv;}
  void setTime (double t) {time = t;}
};

//rubni i inicijalni uvjeti za u_1 = u_t
template<typename GV, typename RF>
class BCExtension1
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, BCExtension1<GV,RF> > {
  const GV& gv;
  double time;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;
  BCExtension1 (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    auto x = e.geometry().global(xlocal);
    if (time==0) y=-cos(x[0]);
    else if(x[0]<1e-5||x[0]>3*M_PI-1e-5) y = -cos(-time);

  //  y=0;

    return;
  }
  inline const GV& getGridView () {return gv;}
  void setTime (double t) {time = t;}
};

//Klasa koja određuje Dirichletove rubne uvjete
template <typename GV>
class BCTypeParam : public Dune::PDELab::DirichletConstraintsParameters
{
  double time;
  const GV& gv;
public:
BCTypeParam(const GV& gv_) : gv(gv_) {}
public:
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
    return true;
  }
  inline const GV& getGridView () {return gv;}
  void setTime (double t) {time = t;}
};
#endif

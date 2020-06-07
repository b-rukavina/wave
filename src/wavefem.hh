#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/type.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>

//sadrži i prostorni i vremenski operator

/** a local operator for the spatial part of the wave equation as a first order system in time
 *
 * \f{align*}{
 *   -c^2\Delta u_1(x) &=& 0  x\in\Omega, \\
 *             -u_0(x) &=& 0  x\in\Omega, \\
 *              u_1(x) &=& 0  x\in\partial\Omega \\
 *  no boundary condition for u_1
 * \f}
 *
 * We assume that both components use the same finite element space.
 * \tparam FEM is a finite element map of the space used for both components
 */
template<typename U_BCType,typename BCType, typename FEM>
class WaveFEM :
  public Dune::PDELab::NumericalJacobianApplyVolume  <WaveFEM<U_BCType,BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianVolume       <WaveFEM<U_BCType,BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianApplyBoundary<WaveFEM<U_BCType,BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianBoundary     <WaveFEM<U_BCType,BCType,FEM>>,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  // types
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  using RF = typename LocalBasis::Traits::RangeFieldType;

  // private data members
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
  RF c;

public:
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };

  WaveFEM (U_BCType& ubctype_,BCType& bctype_,RF c_) : ubctype(ubctype_),bctype(bctype_),c(c_) {}

  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    assert(LFSU::CHILDREN==2);
    using namespace Dune::Indices;
    auto lfsu0 = lfsu.child(_0);
    auto lfsu1 = lfsu.child(_1);

    const int dim = EG::Entity::dimension;

    auto geo = eg.geometry();
    int order = 2*lfsu0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    for (const auto& ip : rule)
      {
        // evaluate basis functions, assume basis is the same for both components
        auto& phihat = cache.evaluateFunction(ip.position(),
                          lfsu0.finiteElement().localBasis());

        // evaluate u1
        RF u1=0.0;
        for (size_t i=0; i<lfsu1.size(); i++)
          u1 += x(lfsu1,i)*phihat[i];

        auto& gradphihat = cache.evaluateJacobian(ip.position(),
                          lfsu0.finiteElement().localBasis());

        const auto &jac =
          geo.jacobianInverseTransposed(ip.position());
        auto gradphi = makeJacobianContainer(lfsu0);
        for (std::size_t i=0; i<lfsu0.size(); i++)
          jac.mv(gradphihat[i][0],gradphi[i][0]);

        // gradient of u0
        Dune::FieldVector<RF,dim> gradu0(0.0), glob(0.0);
        for (std::size_t i=0; i<lfsu0.size(); i++)
          gradu0.axpy(x(lfsu0,i),gradphi[i][0]);

    //    glob = geo.global(ip.position());
    //    double f=5*glob[0]*glob[1];
        double f = 0;

        RF factor = ip.weight()*
          geo.integrationElement(ip.position());
        for (std::size_t i=0; i<lfsu0.size(); i++) {
          r.accumulate(lfsu0,i,c*c*(gradu0*gradphi[i][0])*factor-f*phihat[i]*factor);
          r.accumulate(lfsu1,i,-u1*phihat[i]*factor);
        }
      }
  }

    void preStep (double time, double dt, int stages) {
      bctype.setTime(time+dt);
      using namespace Dune::Indices;
      ubctype.child(_0).setTime(time+dt);
      ubctype.child(_1).setTime(time+dt);
      Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::preStep(time,dt,stages);
    }
    void setTime (double time){
      bctype.setTime(time);
      using namespace Dune::Indices;
      ubctype.child(_0).setTime(time);
      ubctype.child(_1).setTime(time);
      Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::setTime(time);
    }
private:
  BCType & bctype;
  U_BCType & ubctype;
};


/** a local operator for the temporal operator in the wave equation assuming identical components
 *
 * \f{align*}{
 \int_\Omega uv dx
 * \f}
 * \tparam FEM      Type of a finite element map
 */
template<typename FEM>
class WaveL2:
    public Dune::PDELab::NumericalJacobianApplyVolume<WaveL2<FEM>>,
    public Dune::PDELab::NumericalJacobianVolume<WaveL2<FEM>>,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  using RF = typename LocalBasis::Traits::RangeFieldType;

  Dune::PDELab::LocalBasisCache<LocalBasis> cache;

public:

  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };

  WaveL2 () : time(0.0){}

  void setTime(double t) {time = t;}

  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {

    using namespace Dune::Indices;
    auto lfsu0 = lfsu.child(_0);
    auto lfsu1 = lfsu.child(_1);

    auto geo = eg.geometry();
    const int order = 2*lfsu0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    for (const auto& ip : rule)
      {
        // evaluate basis functions at first child
        auto& phihat = cache.evaluateFunction(ip.position(),
                         lfsu0.finiteElement().localBasis());

        // evaluate u0
        RF u0=0.0, u1=0.0;
        for (std::size_t i=0; i<lfsu0.size(); i++) {
          u0 += x(lfsu0,i)*phihat[i];
          u1 += x(lfsu1,i)*phihat[i];
        }

        RF factor=ip.weight()*geo.integrationElement(ip.position());
        for (std::size_t i=0; i<lfsu0.size(); i++) {
          r.accumulate(lfsu0,i,u1*phihat[i]*factor);
          r.accumulate(lfsu1,i,u0*phihat[i]*factor);
        }
      }
  }
private:
  double time;
};

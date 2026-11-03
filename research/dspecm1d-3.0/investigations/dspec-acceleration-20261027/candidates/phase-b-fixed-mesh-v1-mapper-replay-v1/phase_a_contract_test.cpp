#include "phase_a_reconstruction.hpp"

#include <cmath>
#include <stdexcept>

namespace {
void near(double got,double expected) {
  if(!(std::isfinite(got)&&std::abs(got-expected)<1e-13))
    throw std::runtime_error("Phase A synthetic contract assertion failed");
}
}

int main() {
  using namespace phasea;

  const auto dilation=strainOf({Complex{2.0,0.0},Complex{1.0,0.0},
      Complex{},Complex{},Complex{},Complex{}},0,2.0,false);
  near(dilation[0].real(),1.0);
  near(dilation[1].real(),2.0);
  const auto translation=strainOf({Complex{1.0,0.0},Complex{},Complex{1.0,0.0},
      Complex{},Complex{},Complex{}},1,3.0,false);
  const auto rigidRotation=strainOf({Complex{},Complex{},Complex{},Complex{},
      Complex{4.0,0.0},Complex{2.0,0.0}},1,2.0,false);
  for(std::size_t c=0;c<Components;++c){
    near(std::abs(translation[c]),0.0);
    near(std::abs(rigidRotation[c]),0.0);
  }

  StrainTable l0(1,2,1);
  RadialSideDegreeRecord radial;
  radial.degree=0;
  radial.radiusNormalized=2.0;
  radial.familyMask[2]=static_cast<std::uint8_t>(Family::Radial);
  radial.radial[2*6+0]=ComplexPair::from({2.0,0.0});
  radial.radial[2*6+1]=ComplexPair::from({3.0,0.0});
  assignRadialMode(l0,0,radial,2,FieldRole::ForwardSource);
  near(l0(0,0,0,0).real(),3.0);
  near(l0(0,1,0,0).real(),2.0);

  StrainTable l1(1,2,1);
  RadialSideDegreeRecord mixed;
  mixed.degree=1;
  mixed.radiusNormalized=3.0;
  mixed.familyMask[2]=static_cast<std::uint8_t>(Family::Spheroidal) |
                      static_cast<std::uint8_t>(Family::Toroidal);
  mixed.radial[2*6+0]=ComplexPair::from({1.0,0.0});
  mixed.radial[2*6+1]=ComplexPair::from({0.5,0.0});
  mixed.radial[2*6+2]=ComplexPair::from({2.0,0.0});
  mixed.radial[2*6+3]=ComplexPair::from({0.25,0.0});
  mixed.radial[2*6+4]=ComplexPair::from({4.0,0.0});
  mixed.radial[2*6+5]=ComplexPair::from({0.75,0.0});
  assignRadialMode(l1,0,mixed,2,FieldRole::ForwardSource);
  near(l1(0,4,0,1).real(),0.0);
  near(l1(0,5,0,1).real(),0.0);

  const Strain raw{Complex{1.0,0.0},Complex{2.0,0.0},Complex{3.0,0.0},
                   Complex{4.0,0.0},Complex{5.0,0.0},Complex{6.0,0.0}};
  const auto fluid=inheritedFluidProjection(raw);
  near(fluid[0].real(),3.0);
  for(std::size_t c=1;c<Components;++c) near(fluid[c].real(),0.0);

  const auto direct= strainOf({Complex{2.0,0.0},Complex{3.0,0.0},
      Complex{},Complex{},Complex{},Complex{}},0,2.0,false);
  near(direct[0].real(),3.0);
  near(direct[1].real(),2.0);
  const auto ds=strainOfDS({Complex{2.0,0.0},Complex{3.0,0.0},
      Complex{},Complex{},Complex{},Complex{}},0,2.0,false);
  for(std::size_t c=0;c<Components;++c)near(std::abs(direct[c]-ds[c]),0.0);

  const auto contraction=[](const Strain &f,const Strain &a){
    const Complex F=-(std::conj(a[0])*f[1]+std::conj(a[1])*f[0]);
    const Complex L=-2.0*(std::conj(a[2])*f[2]+std::conj(a[3])*f[3]);
    const Complex N=-(.5*(std::conj(a[4])*f[4]+std::conj(a[5])*f[5])-
        std::conj(a[1])*f[1]);
    return -2.0*F+L+N;
  };
  const Strain f{Complex{1.0,2.0},Complex{2.0,-1.0},Complex{-.5,.25},
      Complex{.75,1.5},Complex{2.0,.5},Complex{-.25,1.0}};
  const Strain a{Complex{.5,-.5},Complex{-1.0,2.0},Complex{2.0,1.0},
      Complex{-.75,.5},Complex{1.0,-2.0},Complex{.5,.25}};
  const auto rotate=[](const Strain &v,Complex phase){
    Strain out=v;out[2]*=phase;out[3]*=std::conj(phase);
    out[4]*=phase*phase;out[5]*=std::conj(phase*phase);return out;
  };
  const Complex phase=std::polar(1.0,.713);
  near(std::abs(contraction(f,a)-contraction(rotate(f,phase),rotate(a,phase))),0.0);
  return 0;
}

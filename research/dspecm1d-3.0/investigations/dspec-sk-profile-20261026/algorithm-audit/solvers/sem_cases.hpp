#pragma once

#include "solver_backends.hpp"

#include <cmath>
#include <stdexcept>

namespace dssolver {
enum class Family { Radial, Spheroidal, Toroidal };

struct ActualOperatorCase {
  Sparse matrix;
  Dense sourceRhs;
  ElementPartitionInput elements;
  Family family=Family::Radial;
  int degree=0;
  int rawBin=-1;
  Complex z{};
};

/** Build H and the discrete source-treatment RHS directly from the prepared
 * SEM. This keeps solver comparisons on the real accepted matrix and source
 * lift/load construction rather than synthetic test matrices. */
template<class Sem,class Source>
ActualOperatorCase makeActualOperatorCase(Sem &sem,Source &source,
                                          Family family,int degree,Complex z) {
  if(!std::isfinite(z.real())||!std::isfinite(z.imag())||z==Complex{})
    throw std::invalid_argument("actual SEM case requires finite nonzero frequency");
  ActualOperatorCase out;out.family=family;out.degree=degree;out.z=z;
  switch(family){
  case Family::Radial:{
    if(degree!=0)throw std::invalid_argument("radial case degree must be zero");
    const auto stiffness=sem.hR(),mass=sem.pR();
    out.matrix=stiffness.template cast<Complex>();
    out.matrix-=z*z*mass.template cast<Complex>();
    out.sourceRhs=sem.sourceTreatedRadial(source,z).force;
    out.elements=radialElements(sem);break;
  }
  case Family::Spheroidal:{
    if(degree<1||degree>sem.maxDegree())throw std::invalid_argument("invalid spheroidal case degree");
    const auto stiffness=sem.hS(degree),mass=sem.pS(degree);
    out.matrix=stiffness.template cast<Complex>();
    out.matrix-=z*z*mass.template cast<Complex>();
    out.sourceRhs=sem.sourceTreatedSpheroidal(source,degree,z).force;
    out.elements=spheroidalElements(sem);break;
  }
  case Family::Toroidal:{
    if(degree<1||degree>sem.maxDegree())throw std::invalid_argument("invalid toroidal case degree");
    const auto stiffness=sem.hTk(degree),mass=sem.pTk(degree);
    out.matrix=stiffness.template cast<Complex>();
    out.matrix-=z*z*mass.template cast<Complex>();
    out.sourceRhs=sem.sourceTreatedToroidal(source,degree,z).force;
    out.elements=toroidalElements(sem);break;
  }
  }
  out.matrix.makeCompressed();
  bool matrixFinite=true;
  for(int col=0;col<out.matrix.outerSize();++col)
    for(Sparse::InnerIterator it(out.matrix,col);it;++it)
      matrixFinite=matrixFinite&&std::isfinite(std::abs(it.value()));
  if(out.matrix.rows()!=out.matrix.cols()||out.matrix.rows()!=out.sourceRhs.rows()||
     out.sourceRhs.cols()==0||!matrixFinite||!out.sourceRhs.allFinite())
    throw std::runtime_error("actual SEM operator or source RHS is invalid");
  return out;
}
} // namespace dssolver

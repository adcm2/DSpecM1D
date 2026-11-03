#ifndef DSPECM1D_PHASE_A_RECONSTRUCTION_HPP
#define DSPECM1D_PHASE_A_RECONSTRUCTION_HPP

#include "phase_a_contract.hpp"
#include "engine_headers/Windowed3DFields.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <span>
#include <vector>

namespace phasea {

inline constexpr std::size_t Components = 6;
inline constexpr std::array<int, Components> Spin{0, 0, 1, -1, 2, -2};
using Strain = std::array<Complex, Components>;

/** Literal six-spin radial strain formula copied from pinned YSpec StrainOf.
 * No normalization is added here: DS radialToStrain and YSpec StrainOf use the
 * identity slot map [U,U',V,V',W,W'] and the same normalized radius. */
inline Strain strainOf(const std::array<Complex, 6> &a, int l, double r,
                       bool fluid) {
  if (!(std::isfinite(r) && r > 0.0))
    throw std::invalid_argument("six-spin reconstruction requires r_nd > 0");
  const double lambda = l * (l + 1.0);
  const double zeta = std::sqrt(lambda);
  const double twist = l >= 2 ? zeta * std::sqrt((l - 1.0) * (l + 2.0)) / r : 0.0;
  const Complex U=a[0], Up=a[1], V=a[2], Vp=a[3], W=a[4], Wp=a[5];
  Strain out{};
  out[0] = Up;
  out[1] = (2.0 * U - lambda * V) / r;
  if (!fluid) {
    const Complex shear = Vp - (V - U) / r;
    const Complex turn = Wp - W / r;
    out[2] = -zeta / 2.0 * (shear + Complex(0.0, 1.0) * turn);
    out[3] = zeta / 2.0 * (shear - Complex(0.0, 1.0) * turn);
  }
  out[4] = twist * (V + Complex(0.0, 1.0) * W);
  out[5] = twist * (V - Complex(0.0, 1.0) * W);
  return out;
}

/** Independent transcription of the accepted DS radialToStrain slot algebra,
 * used only for the same-packet A2 old-DS replay. Angular factors are applied
 * later by the legacy complex-basis contraction. */
inline Strain strainOfDS(const std::array<Complex,6> &a,int l,double r,bool fluid) {
  if (!(std::isfinite(r)&&r>0.0))
    throw std::invalid_argument("DS radial replay requires r_nd > 0");
  const double lambda=l*(l+1.0),zeta=std::sqrt(lambda);
  const double eta=l>=2?std::sqrt((l-1.0)*(l+2.0)):0.0;
  const Complex U=a[0],Up=a[1],V=a[2],Vp=a[3],W=a[4],Wp=a[5];
  const Complex S=Vp-(V-U)/r,T=Wp-W/r;
  Strain out{};out[0]=Up;out[1]=(2.0*U-lambda*V)/r;
  if(!fluid){out[2]=-zeta/2.0*(S+Complex(0.0,1.0)*T);
    out[3]=zeta/2.0*(S-Complex(0.0,1.0)*T);}
  out[4]=zeta*eta/r*(V+Complex(0.0,1.0)*W);
  out[5]=zeta*eta/r*(V-Complex(0.0,1.0)*W);
  return out;
}

/** Match the accepted packed/tensor fluid observable. Raw rr/h are retained
 * elsewhere; the inherited fluid projection is rr+h in slot rr, with all
 * unavailable shear slots zero. */
inline Strain inheritedFluidProjection(const Strain &raw) noexcept {
  Strain out{};
  out[0] = raw[0] + raw[1];
  return out;
}

inline Strain rotateToGlobal(const Strain &v, double a, double b) noexcept {
  const Complex c(a,-b), c2=c*c;
  return {v[0],v[1],c*v[2],std::conj(c)*v[3],c2*v[4],std::conj(c2)*v[5]};
}

/** YSpec StrainTable layout: degree is the contiguous innermost dimension. */
class StrainTable {
public:
  StrainTable(std::size_t nodes, int orders, int lMax)
      : nodes_(nodes), orders_(orders), degrees_(static_cast<std::size_t>(lMax)+1),
        values_(nodes*Components*static_cast<std::size_t>(2*orders+1)*degrees_) {
    if (orders < 0 || orders > 2 || lMax < 0)
      throw std::invalid_argument("invalid Phase A strain table dimensions");
  }
  Complex &operator()(std::size_t node, std::size_t component, int m, int l) {
    return values_.at(index(node,component,m)+static_cast<std::size_t>(l));
  }
  const Complex *row(std::size_t node,std::size_t component,int m) const {
    return values_.data()+index(node,component,m);
  }
  int orders() const noexcept { return orders_; }
  std::size_t degrees() const noexcept { return degrees_; }
  std::size_t nodes() const noexcept { return nodes_; }
private:
  std::size_t nodes_;
  int orders_;
  std::size_t degrees_;
  std::vector<Complex> values_;
  std::size_t index(std::size_t node,std::size_t component,int m) const {
    if (node>=nodes_ || component>=Components || m < -orders_ || m > orders_)
      throw std::out_of_range("Phase A strain table index");
    return ((node*Components+component)*static_cast<std::size_t>(2*orders_+1)+
            static_cast<std::size_t>(m+orders_))*degrees_;
  }
};

/** Degree-inner spin-weighted angular table from the pinned LegendreTable
 * layout. The Wigner values use the same GSHTrans convention as the accepted
 * DSpec engine; longitude phase is applied only by evaluate(). */
class LegendreTable {
public:
  LegendreTable(int lMax,int orders,const std::vector<double> &theta)
      : lMax_(lMax), orders_(orders), thetaCount_(theta.size()),
        degrees_(static_cast<std::size_t>(lMax)+1),
        values_(theta.size()*5*static_cast<std::size_t>(2*orders+1)*degrees_) {
    if (lMax<0 || orders<0 || orders>2)
      throw std::invalid_argument("invalid Phase A Legendre table dimensions");
    for (std::size_t p=0;p<theta.size();++p) {
      if (!(std::isfinite(theta[p]) && theta[p]>=0.0 && theta[p]<=std::numbers::pi))
        throw std::invalid_argument("invalid Phase A colatitude");
      SPARSESPEC::Windowed3DFieldEngine::CachedAngularWigner w(
          lMax,std::min(2,lMax),std::min(2,lMax),std::vector<double>{theta[p]});
      for (int l=0;l<=lMax;++l)
        for (int m=-std::min(l,orders);m<=std::min(l,orders);++m)
          for (int n=-std::min(l,2);n<=std::min(l,2);++n)
            (*this)(p,n,m,l)=w[n,0][l,m];
    }
  }
  const double *row(std::size_t point,int n,int m) const {
    return values_.data()+index(point,n,m);
  }
private:
  int lMax_,orders_;
  std::size_t thetaCount_,degrees_;
  std::vector<double> values_;
  std::size_t index(std::size_t point,int n,int m) const {
    if (point>=thetaCount_ || n < -2 || n > 2 || m < -orders_ || m > orders_)
      throw std::out_of_range("Phase A Legendre table index");
    return ((point*5+static_cast<std::size_t>(n+2))*
            static_cast<std::size_t>(2*orders_+1)+static_cast<std::size_t>(m+orders_))*degrees_;
  }
  double &operator()(std::size_t point,int n,int m,int l) {
    return values_.at(index(point,n,m)+static_cast<std::size_t>(l));
  }
};

/** Literal pinned YSpec Evaluate summation order: m, component, then l. */
inline Strain evaluate(const StrainTable &table,std::size_t node,
                       const LegendreTable &legendre,std::size_t point,
                       double phi) {
  Strain out{};
  for (int m=-table.orders();m<=table.orders();++m) {
    const Complex phase=std::polar(1.0,m*phi);
    for (std::size_t c=0;c<Components;++c) {
      const Complex *f=table.row(node,c,m);
      const double *x=legendre.row(point,Spin[c],m);
      double re=0.0,im=0.0;
      for (std::size_t l=0;l<table.degrees();++l) {
        re+=f[l].real()*x[l];
        im+=f[l].imag()*x[l];
      }
      out[c]+=Complex(re,im)*phase;
    }
  }
  return out;
}

inline void assignRadialMode(StrainTable &table,std::size_t node,
                             const RadialSideDegreeRecord &record,
                             std::size_t slot,FieldRole role) {
  if (record.degree<0 || static_cast<std::size_t>(record.degree)>=table.degrees())
    throw std::out_of_range("Phase A radial record degree exceeds table");
  if (slot>=ModalOrderSlots || (role==FieldRole::ForwardSource && slot>=SourceOrderSlots) ||
      (role==FieldRole::AdjointDual && slot<SourceOrderSlots))
    throw std::out_of_range("Phase A modal slot and role do not match");
  const int m=role==FieldRole::ForwardSource?sourceOrder(slot):adjointOrder(slot);
  if (m < -table.orders() || m > table.orders() || std::abs(m)>record.degree)
    throw std::out_of_range("Phase A radial record order exceeds table");
  const std::uint8_t mask=record.familyMask[slot];
  if (mask==0) throw std::runtime_error("valid Phase A modal slot has no family provenance");
  std::array<Complex,6> radial{};
  for (std::size_t i=0;i<6;++i) radial[i]=record.radial[slot*6+i].value();
  const Complex zero{};
  const auto radialBit=static_cast<std::uint8_t>(Family::Radial);
  const auto sphBit=static_cast<std::uint8_t>(Family::Spheroidal);
  const auto torBit=static_cast<std::uint8_t>(Family::Toroidal);
  if (mask&radialBit) {
    if (mask!=radialBit || record.degree!=0 || m!=0 || radial[2]!=zero || radial[3]!=zero ||
        radial[4]!=zero || radial[5]!=zero)
      throw std::runtime_error("invalid separate radial l=0 packet row");
  } else if (mask & ~(sphBit|torBit)) {
    throw std::runtime_error("invalid family mask in Phase A modal slot");
  } else if ((mask&sphBit) && !(mask&torBit)) {
    if (radial[4]!=zero || radial[5]!=zero)
      throw std::runtime_error("spheroidal packet row contains toroidal slots");
  } else if ((mask&torBit) && !(mask&sphBit)) {
    if (radial[0]!=zero || radial[1]!=zero || radial[2]!=zero || radial[3]!=zero)
      throw std::runtime_error("toroidal packet row contains spheroidal slots");
  }
  auto strain=strainOf(radial,record.degree,record.radiusNormalized,record.fluid!=0);
  for (std::size_t c=0;c<Components;++c)
    table(node,c,m,record.degree)=strain[c];
}

inline void assignRadialModeDS(StrainTable &table,std::size_t node,
                               const RadialSideDegreeRecord &record,
                               std::size_t slot,FieldRole role) {
  if(record.degree<0||static_cast<std::size_t>(record.degree)>=table.degrees()||
     slot>=ModalOrderSlots||(role==FieldRole::ForwardSource&&slot>=SourceOrderSlots)||
     (role==FieldRole::AdjointDual&&slot<SourceOrderSlots))
    throw std::out_of_range("DS replay radial record index/role mismatch");
  const int m=role==FieldRole::ForwardSource?sourceOrder(slot):adjointOrder(slot);
  if(m < -table.orders()||m>table.orders()||std::abs(m)>record.degree)
    throw std::out_of_range("DS replay order exceeds angular table");
  std::array<Complex,6> radial{};
  for(std::size_t c=0;c<6;++c)radial[c]=record.radial[slot*6+c].value();
  const auto strain=strainOfDS(radial,record.degree,record.radiusNormalized,record.fluid!=0);
  table(node,0,m,record.degree)=strain[0];
  for(std::size_t c=1;c<Components;++c)
    table(node,c,m,record.degree)=strain[c];
}

} // namespace phasea

#endif

#include "algorithms.hpp"

#include <GSHTrans/Core>
#include <fftw3.h>

#include <Eigen/QR>
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>

namespace fs = std::filesystem;
using dspectral::Complex;
using Clock = std::chrono::steady_clock;

namespace {
constexpr int slotSpin[6] = {0, 0, 1, -1, 2, -2};
constexpr std::uint64_t radialForwardSet=std::numeric_limits<std::uint64_t>::max();
constexpr std::uint64_t radialAdjointSet=std::numeric_limits<std::uint64_t>::max()-1;
struct Point { double radius, theta, phi; int element; };
struct Label { int m, l; };
struct Pole {
  std::uint64_t set = 0;
  bool adjoint = false;
  int element = -1;
  double radiusModel = 0;
  // The three rows are local theta, east, radial unit vectors in global XYZ.
  std::array<std::array<double,3>,3> axes{};
  std::vector<Label> labels;
  std::array<std::vector<Complex>,6> values;
};
struct Record {
  std::uint32_t objective = 0, bin = 0;
  std::vector<Point> points;
  std::vector<std::array<Complex,12>> fields;
  std::vector<Pole> poles;
};
struct Packet {
  std::uint32_t schema = 0, objectives = 0, frequencies = 0;
  std::uint64_t pointCount = 0;
  std::vector<Record> records;
};

template<class T> T read(std::istream &in) {
  T value{};
  in.read(reinterpret_cast<char*>(&value), sizeof(T));
  if (!in) throw std::runtime_error("truncated solver_packet.bin");
  return value;
}
Complex readComplex(std::istream &in) {
  const double re = read<double>(in), im = read<double>(in);
  return {re, im};
}
void skipMatrix(std::istream &in) {
  const auto rows = read<std::uint32_t>(in), cols = read<std::uint32_t>(in);
  const std::uint64_t n = static_cast<std::uint64_t>(rows) * cols;
  in.seekg(static_cast<std::streamoff>(n * 2 * sizeof(double)), std::ios::cur);
  if (!in) throw std::runtime_error("truncated modal matrix in solver_packet.bin");
}
Packet readPacket(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open solver packet");
  std::array<char,8> magic{};
  in.read(magic.data(), magic.size());
  if (!in || std::string(magic.data(), magic.size()) != "DSFASTPK")
    throw std::runtime_error("invalid solver packet magic");
  Packet p;
  p.schema = read<std::uint32_t>(in);
  p.objectives = read<std::uint32_t>(in);
  p.frequencies = read<std::uint32_t>(in);
  p.pointCount = read<std::uint64_t>(in);
  if (p.schema != 1 || p.objectives == 0 || p.frequencies == 0)
    throw std::runtime_error("unsupported solver packet header");
  p.records.reserve(static_cast<std::size_t>(p.objectives) * p.frequencies);
  for (std::uint64_t ri = 0; ri < static_cast<std::uint64_t>(p.objectives)*p.frequencies; ++ri) {
    Record r;
    r.objective = read<std::uint32_t>(in); r.bin = read<std::uint32_t>(in);
    const auto np = read<std::uint64_t>(in);
    const auto nb = read<std::uint32_t>(in), npr = read<std::uint32_t>(in);
    if (np != p.pointCount) throw std::runtime_error("record point count changed");
    r.points.resize(np); r.fields.resize(np);
    for (std::size_t i=0; i<np; ++i) {
      r.points[i] = {read<double>(in), read<double>(in), read<double>(in), read<std::int32_t>(in)};
      for (auto &x : r.fields[i]) x = readComplex(in);
    }
    std::map<std::tuple<std::uint64_t,bool,int,std::uint64_t>,Pole> mergedPoles;
    for (std::uint32_t j=0; j<npr; ++j) {
      Pole q;
      q.set = read<std::uint64_t>(in);
      q.adjoint = read<std::uint8_t>(in) != 0;
      q.element = read<std::int32_t>(in); q.radiusModel = read<double>(in);
      for (auto &axis : q.axes) for (double &x : axis) x = read<double>(in);
      const auto nq = read<std::uint64_t>(in);
      q.labels.resize(nq);
      for (std::size_t k=0; k<nq; ++k) {
        q.labels[k] = {read<std::int32_t>(in), read<std::int32_t>(in)};
        for (auto &slot : q.values) slot.push_back(readComplex(in));
      }
      const auto key=std::tuple{q.set,q.adjoint,q.element,std::bit_cast<std::uint64_t>(q.radiusModel)};
      auto [it,inserted]=mergedPoles.try_emplace(key,std::move(q));
      if(!inserted) {
        const Pole &prior=it->second;
        const Pole &part=q;
        if(prior.labels.size()!=part.labels.size())
          throw std::runtime_error("duplicate pole coefficient layout differs");
        for(std::size_t k=0;k<prior.labels.size();++k)
          if(prior.labels[k].m!=part.labels[k].m||prior.labels[k].l!=part.labels[k].l)
            throw std::runtime_error("duplicate pole coefficient labels differ");
        for(int slot=0;slot<6;++slot) for(std::size_t k=0;k<prior.labels.size();++k)
          if(prior.values[slot][k]!=part.values[slot][k])
            throw std::runtime_error("repeated angular chunk coefficient records are inconsistent");
      }
    }
    r.poles.reserve(mergedPoles.size());
    for(auto &[key,pole]:mergedPoles) { (void)key; r.poles.push_back(std::move(pole)); }
    for (std::uint32_t j=0; j<nb; ++j) {
      (void)read<std::int32_t>(in); (void)read<std::int32_t>(in); (void)read<std::int32_t>(in);
      for (int m=0; m<5; ++m) skipMatrix(in);
    }
    p.records.push_back(std::move(r));
  }
  return p;
}

std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
using TraceKey = std::pair<int,std::uint64_t>;
using AngleKey = std::pair<std::uint64_t,std::uint64_t>;
TraceKey traceKey(const Point &p) { return {p.element,bits(p.radius)}; }
AngleKey angleKey(const Point &p) { return {bits(p.theta),bits(p.phi)}; }

std::array<double,3> unit(double theta,double phi) {
  return {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)};
}
double dot(const std::array<double,3> &a,const std::array<double,3> &b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
std::array<double,3> frameAngle(const Pole &pole,double theta,double phi,
                                double &localTheta,double &localPhi) {
  const auto x = unit(theta,phi);
  std::array<double,3> local{};
  for (int i=0;i<3;++i) local[i]=dot(pole.axes[i],x);
  const double horizontal=std::hypot(local[0],local[1]);
  localTheta=std::atan2(horizontal,local[2]);
  localPhi=horizontal==0.0?0.0:std::atan2(local[1],local[0]);
  const std::array<double,3> localEtComponents={
      std::cos(localTheta)*std::cos(localPhi),
      std::cos(localTheta)*std::sin(localPhi),-std::sin(localTheta)};
  std::array<double,3> localEt{};
  for (int i=0;i<3;++i) for (int j=0;j<3;++j)
    localEt[j]+=localEtComponents[i]*pole.axes[i][j];
  const std::array<double,3> globalEt={std::cos(theta)*std::cos(phi),
      std::cos(theta)*std::sin(phi),-std::sin(theta)};
  const std::array<double,3> globalEp={-std::sin(phi),std::cos(phi),0.0};
  const double a=dot(globalEt,localEt), b=dot(globalEp,localEt);
  return {a,b,0.0};
}

using Wigner = GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,GSHTrans::Multiple>;
Eigen::MatrixXcd makeBasis(const Pole &pole,const std::vector<Label> &labels,
                           int spin,const std::vector<Point> &angles,int lmax) {
  std::vector<double> localTheta(angles.size()), localPhi(angles.size());
  std::vector<Complex> rotation(angles.size());
  for (std::size_t p=0;p<angles.size();++p) {
    double t=0, ph=0; const auto ab=frameAngle(pole,angles[p].theta,angles[p].phi,t,ph);
    localTheta[p]=t; localPhi[p]=ph;
    Complex c{ab[0],-ab[1]}, cp{1.0,0.0};
    for(int k=0;k<std::abs(spin);++k) cp*=spin<0?std::conj(c):c;
    rotation[p]=cp;
  }
  Wigner w(lmax,lmax,std::min(2,lmax),localTheta);
  Eigen::MatrixXcd h(angles.size(),labels.size());
  for (std::size_t p=0;p<angles.size();++p) for(std::size_t j=0;j<labels.size();++j) {
    const auto [m,l]=labels[j];
    h(p,j) = l < std::abs(spin) ? Complex{} :
        rotation[p]*w[spin,static_cast<std::ptrdiff_t>(p)][l,m]*
        std::exp(Complex(0.0,m*localPhi[p]));
  }
  return h;
}

std::vector<Label> labelsForSpin(int lmax,int spin) {
  std::vector<Label> out;
  for(int l=std::abs(spin);l<=lmax;++l) for(int m=-l;m<=l;++m) out.push_back({m,l});
  return out;
}
std::vector<Label> globalLabels(const std::vector<Label> &local,int lmax) {
  (void)local;
  return labelsForSpin(lmax,0);
}
Eigen::MatrixXcd makeGlobalBasis(const std::vector<Point> &angles,int spin,int lmax,
                                 const std::vector<Label> &labels) {
  std::vector<double> theta(angles.size());
  for(std::size_t i=0;i<angles.size();++i) theta[i]=angles[i].theta;
  Wigner w(lmax,lmax,std::min(2,lmax),theta);
  Eigen::MatrixXcd h(angles.size(),labels.size());
  for(std::size_t p=0;p<angles.size();++p) for(std::size_t j=0;j<labels.size();++j) {
    const auto [m,l]=labels[j];
    h(p,j)=l<std::abs(spin)?Complex{}:w[spin,static_cast<std::ptrdiff_t>(p)][l,m]*
        std::exp(Complex(0,m*angles[p].phi));
  }
  return h;
}

struct Err { double l2=0, scaled=0; };
Err compare(const Eigen::VectorXcd &ref,const Eigen::VectorXcd &test) {
  const auto d=test-ref;
  Err e; e.l2=ref.norm()==0?d.norm():d.norm()/ref.norm();
  const double floor=0.01*ref.cwiseAbs().maxCoeff();
  if(floor==0.0) e.scaled=d.cwiseAbs().maxCoeff();
  else for(Eigen::Index i=0;i<ref.size();++i)
    e.scaled=std::max(e.scaled,std::abs(d[i])/std::max(std::abs(ref[i]),floor));
  return e;
}

struct TraceData {
  TraceKey key;
  std::vector<Point> angles;
  std::vector<std::vector<std::size_t>> pointRows; // frequency -> packet point index by angle
  std::vector<const Record*> records;
  std::vector<std::map<std::uint64_t,const Pole*>> poles; // frequency -> synthesisSet -> pole
};

const Pole *findPole(const Record &r,std::uint64_t set) {
  for(const auto &p:r.poles) if(p.set==set) return &p;
  return nullptr;
}
std::vector<std::uint64_t> adjointSets(const Record &r) {
  std::set<std::uint64_t> ids;
  for(const auto &p:r.poles) if(p.adjoint&&p.set!=radialAdjointSet) ids.insert(p.set);
  return {ids.begin(),ids.end()};
}

std::vector<TraceData> makeTraces(const Packet &packet,std::uint32_t objective) {
  std::vector<const Record*> records;
  for(const auto &r:packet.records) if(r.objective==objective) records.push_back(&r);
  std::sort(records.begin(),records.end(),[](auto a,auto b){return a->bin<b->bin;});
  if(records.size()!=packet.frequencies) throw std::runtime_error("objective frequency records missing");
  std::map<TraceKey,std::vector<Point>> groups;
  for(const auto &p:records.front()->points) groups[traceKey(p)].push_back(p);
  std::map<std::pair<int,std::uint64_t>,std::size_t> poleTrace;
  std::map<int,std::vector<std::pair<double,std::size_t>>> traceOrder;
  std::map<int,std::vector<std::pair<double,std::uint64_t>>> poleOrder;
  for(const auto &[key,unused]:groups) { (void)unused; traceOrder[key.first].push_back({std::bit_cast<double>(key.second),0}); }
  for(auto &[element,items]:traceOrder) {
    std::sort(items.begin(),items.end(),[](auto a,auto b){return a.first<b.first;});
    std::size_t i=0; for(auto &item:items) item.second=i++;
  }
  for(const auto &p:records.front()->poles)
    if(traceOrder.contains(p.element))
      poleOrder[p.element].push_back({p.radiusModel,bits(p.radiusModel)});
  for(auto &[element,items]:poleOrder) {
    std::sort(items.begin(),items.end(),[](auto a,auto b){return a.first<b.first;});
    items.erase(std::unique(items.begin(),items.end(),[](auto a,auto b){return a.second==b.second;}),items.end());
    auto &ti=traceOrder[element];
    if(items.size()!=ti.size()) throw std::runtime_error("pole and point radial traces differ");
    for(std::size_t i=0;i<items.size();++i) poleTrace[{element,items[i].second}]=i;
  }
  std::vector<TraceData> traces;
  for(auto &[key,points]:groups) {
    std::sort(points.begin(),points.end(),[](const auto&a,const auto&b){return std::tie(a.theta,a.phi)<std::tie(b.theta,b.phi);});
    points.erase(std::unique(points.begin(),points.end(),[](const auto&a,const auto&b){return angleKey(a)==angleKey(b);}),points.end());
    TraceData td; td.key=key; td.angles=points; td.records=records;
    td.pointRows.resize(records.size()); td.poles.resize(records.size());
    for(std::size_t k=0;k<records.size();++k) {
      std::map<AngleKey,std::size_t> rowByAngle;
      for(std::size_t i=0;i<records[k]->points.size();++i)
        if(traceKey(records[k]->points[i])==key) rowByAngle.emplace(angleKey(records[k]->points[i]),i);
      td.pointRows[k].resize(td.angles.size());
      for(std::size_t i=0;i<td.angles.size();++i) {
        const auto it=rowByAngle.find(angleKey(td.angles[i]));
        if(it==rowByAngle.end()) throw std::runtime_error("frequency has mismatched angular samples");
        td.pointRows[k][i]=it->second;
      }
      const auto rankIt=std::find_if(traceOrder[key.first].begin(),traceOrder[key.first].end(),
          [&](auto item){return bits(item.first)==key.second;});
      if(rankIt==traceOrder[key.first].end()) throw std::runtime_error("point trace rank missing");
      const std::size_t localRank=rankIt->second;
      for(const auto &p:records[k]->poles) {
        if(p.element!=key.first) continue;
        const auto pt=poleTrace.find({p.element,bits(p.radiusModel)});
        if(pt!=poleTrace.end()&&pt->second==localRank) td.poles[k].emplace(p.set,&p);
      }
    }
    traces.push_back(std::move(td));
  }
  return traces;
}

struct PairOutput {
  Eigen::MatrixXcd covariance;
  Eigen::MatrixXcd basisForward;
  Eigen::VectorXcd direct, directSynthesized, pair;
  std::vector<Eigen::MatrixXcd> basisAdjoint;
  std::vector<Eigen::MatrixXcd> receiverCovariances;
  std::vector<Eigen::VectorXcd> receiverPair;
  double directSynthesisSeconds=0, covarianceFormationSeconds=0, covarianceEvaluationSeconds=0;
};

struct PoleBasis {
  std::array<Eigen::MatrixXcd,6> local;
  std::array<Eigen::MatrixXcd,6> global;
  std::array<Eigen::MatrixXcd,6> rotation;
  std::array<std::vector<Label>,6> globalLabel;
};
PoleBasis makePoleBasis(const Pole &pole,const std::vector<Point> &angles,int lmax,
                        double &compactSeconds,double &globalFrameSeconds) {
  PoleBasis out;
  for(int slot=0;slot<6;++slot) {
    const int spin=slotSpin[slot];
    const auto compactStart=Clock::now();
    out.local[slot]=makeBasis(pole,pole.labels,spin,angles,lmax);
    compactSeconds+=std::chrono::duration<double>(Clock::now()-compactStart).count();
    out.globalLabel[slot]=labelsForSpin(lmax,spin);
    const auto globalStart=Clock::now();
    out.global[slot]=makeGlobalBasis(angles,spin,lmax,out.globalLabel[slot]);
    out.rotation[slot]=out.global[slot].colPivHouseholderQr().solve(out.local[slot]);
    globalFrameSeconds+=std::chrono::duration<double>(Clock::now()-globalStart).count();
  }
  return out;
}
std::uint64_t matrixPayloadBytes(const Eigen::MatrixXcd &m) {
  return static_cast<std::uint64_t>(m.size())*sizeof(Complex);
}
std::uint64_t compactBasisBytes(const PoleBasis &b) {
  std::uint64_t n=0;
  for(const auto &m:b.local) n+=matrixPayloadBytes(m);
  return n;
}
std::uint64_t globalBasisBytes(const PoleBasis &b) {
  std::uint64_t n=0;
  for(int i=0;i<6;++i) {
    n+=matrixPayloadBytes(b.global[i])+matrixPayloadBytes(b.rotation[i]);
    n+=static_cast<std::uint64_t>(b.globalLabel[i].capacity())*sizeof(Label);
  }
  return n;
}
Eigen::VectorXcd poleCoefficients(const Pole *pole,int slot) {
  if(!pole) throw std::runtime_error("missing pole coefficient row");
  Eigen::VectorXcd v(pole->values[slot].size());
  for(std::size_t i=0;i<pole->values[slot].size();++i) v[i]=pole->values[slot][i];
  return v;
}
Eigen::MatrixXcd poleFrequencyCoefficients(const TraceData &td,std::uint64_t set,
    int slot,const PoleBasis &basis,bool rotateGlobal) {
  const Pole *first=td.poles.front().at(set);
  const std::size_t q=rotateGlobal?basis.globalLabel[slot].size():first->labels.size();
  Eigen::MatrixXcd c(td.records.size(),q);
  for(std::size_t k=0;k<td.records.size();++k) {
    auto v=poleCoefficients(td.poles[k].at(set),slot);
    if(rotateGlobal) v=basis.rotation[slot]*v;
    if(v.size()!=q) throw std::runtime_error("coefficient basis dimension changed");
    c.row(k)=v.transpose();
  }
  return c;
}

PairOutput evaluatePair(const TraceData &td,std::uint64_t objective,
    int aslot,int fslot,const std::vector<Eigen::MatrixXcd> &basisA,
    const Eigen::MatrixXcd &basisF,const std::vector<std::uint64_t> &receivers) {
  const int nf=td.records.size();
  (void)objective;
  const Pole *firstF=td.poles.front().at(0);
  if(!firstF) throw std::runtime_error("source pole record missing");
  const std::size_t qf=firstF->labels.size()+(fslot<2?1:0);
  Eigen::MatrixXcd cf(nf,qf);
  for(int k=0;k<nf;++k) {
    const Pole *p=td.poles[k].at(0);
    if(!p||p->labels.size()!=firstF->labels.size()) throw std::runtime_error("source coefficient layout changed");
    const std::size_t qCompact=p->labels.size();
    for(std::size_t j=0;j<qCompact;++j) cf(k,j)=p->values[fslot][j];
    if(fslot<2) cf(k,qCompact)=td.poles[k].at(radialForwardSet)->values[fslot][0];
  }
  PairOutput out; out.direct=Eigen::VectorXcd::Zero(td.angles.size());
  out.directSynthesized=Eigen::VectorXcd::Zero(td.angles.size());
  out.basisForward=basisF;out.basisAdjoint=basisA;
  for(int k=0;k<nf;++k) for(std::size_t p=0;p<td.angles.size();++p) {
    const auto &field=td.records[k]->fields[td.pointRows[k][p]];
    out.direct[p]+=std::conj(field[6+aslot])*field[fslot];
  }
  out.pair=Eigen::VectorXcd::Zero(td.angles.size());
  Eigen::MatrixXcd fFields;
  {
    const auto t0=Clock::now();
    fFields=basisF*cf.transpose();
    out.directSynthesisSeconds+=std::chrono::duration<double>(Clock::now()-t0).count();
  }
  for(std::size_t ri=0;ri<receivers.size();++ri) {
    const Pole *firstA=td.poles.front().at(receivers[ri]);
    if(!firstA) throw std::runtime_error("adjoint pole record missing");
    const std::size_t qa=firstA->labels.size()+(aslot<2?1:0);
    Eigen::MatrixXcd ca(nf,qa);
    for(int k=0;k<nf;++k) {
      const Pole *p=td.poles[k].at(receivers[ri]);
      if(!p||p->labels.size()!=firstA->labels.size()) throw std::runtime_error("adjoint coefficient layout changed");
      const std::size_t qCompact=p->labels.size();
      for(std::size_t j=0;j<qCompact;++j) ca(k,j)=p->values[aslot][j];
      if(aslot<2) ca(k,qCompact)=td.poles[k].at(radialAdjointSet)->values[aslot][0];
    }
    Eigen::MatrixXcd aFields;
    {
      const auto t0=Clock::now();
      aFields=basisA[ri]*ca.transpose();
      out.directSynthesized.array()+=(aFields.conjugate().array()*fFields.array()).rowwise().sum();
      out.directSynthesisSeconds+=std::chrono::duration<double>(Clock::now()-t0).count();
    }
    Eigen::MatrixXcd c;
    {
      const auto t0=Clock::now();
      c=dspectral::frequencyCovariance(ca,cf);
      out.covarianceFormationSeconds+=std::chrono::duration<double>(Clock::now()-t0).count();
    }
    Eigen::VectorXcd v;
    {
      const auto t0=Clock::now();
      v=dspectral::evaluateCovariance(basisA[ri],basisF,c);
      out.covarianceEvaluationSeconds+=std::chrono::duration<double>(Clock::now()-t0).count();
    }
    out.pair+=v;
    out.receiverCovariances.push_back(std::move(c));
    out.receiverPair.push_back(std::move(v));
  }
  return out;
}

std::pair<std::vector<double>,std::vector<double>> gaussLegendre(int n) {
  std::vector<double> x(n),w(n);
  const double eps=2e-15;
  for(int i=0;i<(n+1)/2;++i) {
    double z=std::cos(std::numbers::pi*(i+0.75)/(n+0.5)),z1=0,pp=0;
    do {
      double p1=1,p2=0;
      for(int j=1;j<=n;++j){double p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/j;}
      pp=n*(z*p1-p2)/(z*z-1); z1=z; z=z1-p1/pp;
    } while(std::abs(z-z1)>eps);
    x[i]=-z; x[n-1-i]=z;
    const double wi=2/((1-z*z)*pp*pp); w[i]=wi;w[n-1-i]=wi;
  }
  return {x,w};
}

struct GauntTerm { Eigen::Index a,b,out; double value; };
std::vector<GauntTerm> makeGauntTerms(const std::vector<Label> &a,
    const std::vector<Label> &b,int spin,int lmax) {
  const int maxJ=2*lmax, nlat=2*lmax+1, nmodes=(maxJ+1)*(maxJ+1);
  const auto [x,wg]=gaussLegendre(nlat);
  std::vector<double> theta(nlat);
  for(int i=0;i<nlat;++i) theta[i]=std::acos(x[i]);
  Wigner w(maxJ,maxJ,std::min(2,maxJ),theta);
  std::vector<GauntTerm> terms;
  const double twoPi=2*std::numbers::pi;
  for(std::size_t ia=0;ia<a.size();++ia) for(std::size_t ib=0;ib<b.size();++ib) {
    const int m=a[ia].m, n=b[ib].m, M=n-m;
    const int l=a[ia].l, lp=b[ib].l;
    if(std::abs(M)>maxJ) continue;
    const int jlo=std::max(std::abs(l-lp),std::abs(M)), jhi=l+lp;
    for(int J=jlo;J<=jhi;++J) {
      double integral=0;
      for(int it=0;it<nlat;++it)
        integral+=wg[it]*w[spin,it][l,m]*w[spin,it][lp,n]*w[0,it][J,M];
      // GSHTrans stores sqrt((2J+1)/(4 pi))*d, so each scalar harmonic has
      // unit norm and the quadrature integral is already the projection
      // coefficient. The azimuthal phase contributes the 2 pi factor.
      integral*=twoPi;
      const int idx=J*J+J+M;
      if(std::abs(integral)>2e-15) terms.push_back({static_cast<Eigen::Index>(ia),
          static_cast<Eigen::Index>(ib),idx,integral});
    }
  }
  if(terms.empty()) throw std::runtime_error("no Gaunt terms generated");
  (void)nmodes;
  return terms;
}
Eigen::VectorXcd gauntRecouple(const Eigen::MatrixXcd &c,
    const std::vector<GauntTerm> &terms,int lmax) {
  Eigen::VectorXcd out=Eigen::VectorXcd::Zero((2*lmax+1)*(2*lmax+1));
  for(const auto &t:terms) out[t.out]+=t.value*c(t.a,t.b);
  return out;
}

struct FFTProjection { Eigen::VectorXcd coeff; Eigen::MatrixXd reconstructed; Err error; };
FFTProjection fftProject(const Eigen::MatrixXd &grid,const std::vector<double> &theta,
                         const std::vector<double> &phi,int lmax) {
  const int ntheta=theta.size(), nphi=phi.size(), maxJ=2*lmax;
  if(ntheta<maxJ+1||nphi<2*maxJ+1) throw std::runtime_error("global grid undersamples scalar kernel bandlimit");
  const double dphi=2*std::numbers::pi/nphi;
  for(int j=1;j<nphi;++j) if(std::abs((phi[j]-phi[j-1])-dphi)>1e-10)
    throw std::runtime_error("global longitude samples are not equally spaced");
  for(int i=0;i<ntheta;++i) for(int j=0;j<nphi;++j)
    (void)grid(i,j);
  fftw_complex *in=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nphi);
  fftw_complex *out=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nphi);
  auto fwd=fftw_plan_dft_1d(nphi,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  auto back=fftw_plan_dft_1d(nphi,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
  Eigen::MatrixXcd fourier(ntheta,maxJ*2+1);
  const double phi0=phi.front();
  for(int i=0;i<ntheta;++i) {
    for(int j=0;j<nphi;++j){in[j][0]=grid(i,j);in[j][1]=0;}
    fftw_execute(fwd);
    for(int M=-maxJ;M<=maxJ;++M) {
      const int k=M>=0?M:nphi+M;
      fourier(i,M+maxJ)=Complex(out[k][0],out[k][1])/double(nphi)*std::exp(Complex(0,-M*phi0));
    }
  }
  Wigner w(maxJ,maxJ,std::min(2,maxJ),theta);
  Eigen::VectorXcd coeff=Eigen::VectorXcd::Zero((maxJ+1)*(maxJ+1));
  Eigen::MatrixXcd lat(ntheta,maxJ*2+1);
  for(int M=-maxJ;M<=maxJ;++M) {
    const int j0=std::abs(M), ncol=maxJ-j0+1;
    Eigen::MatrixXcd a(ntheta,ncol);
    for(int i=0;i<ntheta;++i) for(int k=0;k<ncol;++k)
      a(i,k)=w[0,i][j0+k,M];
    const Eigen::VectorXcd c=a.colPivHouseholderQr().solve(fourier.col(M+maxJ));
    for(int k=0;k<ncol;++k) coeff[(j0+k)*(j0+k)+(j0+k)+M]=c[k];
    lat.col(M+maxJ)=a*c;
  }
  Eigen::MatrixXd reconstructed(ntheta,nphi);
  for(int i=0;i<ntheta;++i) {
    for(int k=0;k<nphi;++k){in[k][0]=0;in[k][1]=0;}
    for(int M=-maxJ;M<=maxJ;++M) {
      const int k=M>=0?M:nphi+M;
      const Complex z=double(nphi)*lat(i,M+maxJ)*std::exp(Complex(0,M*phi0));
      in[k][0]=z.real();in[k][1]=z.imag();
    }
    fftw_execute(back);
    for(int k=0;k<nphi;++k) reconstructed(i,k)=out[k][0]/nphi;
  }
  fftw_destroy_plan(fwd);fftw_destroy_plan(back);fftw_free(in);fftw_free(out);
  Eigen::VectorXcd ref=Eigen::Map<const Eigen::VectorXd>(grid.data(),grid.size()).cast<Complex>();
  Eigen::VectorXcd test=Eigen::Map<const Eigen::VectorXd>(reconstructed.data(),reconstructed.size()).cast<Complex>();
  return {std::move(coeff),std::move(reconstructed),compare(ref,test)};
}

Eigen::VectorXcd realHarmonics(const Eigen::VectorXcd &c,int maxJ,double factor=1.0) {
  Eigen::VectorXcd out(c.size());
  for(int J=0;J<=maxJ;++J) for(int M=-J;M<=J;++M) {
    const int i=J*J+J+M, j=J*J+J-M;
    const double phase=(std::abs(M)%2)?-1.0:1.0;
    out[i]=factor*0.5*(c[i]+phase*std::conj(c[j]));
  }
  return out;
}
Eigen::VectorXd synthesizeHarmonics(const Eigen::VectorXcd &c,
    const std::vector<double> &theta,const std::vector<double> &phi,int maxJ) {
  std::vector<double> thetaRepeated=theta;
  Wigner w(maxJ,maxJ,std::min(2,maxJ),thetaRepeated);
  Eigen::VectorXd out(theta.size()*phi.size());
  std::size_t p=0;
  for(std::size_t i=0;i<theta.size();++i) for(double ph:phi) {
    Complex value{};
    for(int J=0;J<=maxJ;++J) for(int M=-J;M<=J;++M)
      value+=c[J*J+J+M]*w[0,i][J,M]*std::exp(Complex(0,M*ph));
    out[p++]=value.real();
  }
  return out;
}

Eigen::VectorXd fftSynthesizeHarmonics(const Eigen::VectorXcd &c,
    const std::vector<double> &theta,const std::vector<double> &phi,int maxJ) {
  const int ntheta=theta.size(),nphi=phi.size();
  if(ntheta<maxJ+1||nphi<2*maxJ+1) throw std::runtime_error("global grid undersamples scalar kernel bandlimit");
  const double dphi=2*std::numbers::pi/nphi,phi0=phi.front();
  for(int j=1;j<nphi;++j) if(std::abs((phi[j]-phi[j-1])-dphi)>1e-10)
    throw std::runtime_error("global longitude samples are not equally spaced");
  fftw_complex *in=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nphi);
  fftw_complex *out=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nphi);
  auto plan=fftw_plan_dft_1d(nphi,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
  Wigner w(maxJ,maxJ,std::min(2,maxJ),theta);
  Eigen::VectorXd values(ntheta*nphi);
  for(int i=0;i<ntheta;++i) {
    for(int k=0;k<nphi;++k){in[k][0]=0;in[k][1]=0;}
    for(int M=-maxJ;M<=maxJ;++M) {
      Complex lat{};
      for(int J=std::abs(M);J<=maxJ;++J)
        lat+=c[J*J+J+M]*w[0,i][J,M];
      const int k=M>=0?M:nphi+M;
      const Complex z=double(nphi)*lat*std::exp(Complex(0,M*phi0));
      in[k][0]=z.real();in[k][1]=z.imag();
    }
    fftw_execute(plan);
    for(int j=0;j<nphi;++j) values[i*nphi+j]=out[j][0]/nphi;
  }
  fftw_destroy_plan(plan);fftw_free(in);fftw_free(out);
  return values;
}

Eigen::MatrixXcd globalFrequencyCoefficients(const TraceData &td,std::uint64_t set,
    int slot,const PoleBasis &basis,std::uint64_t radialSet) {
  auto c=poleFrequencyCoefficients(td,set,slot,basis,true);
  // Reserved l=0 rows are angularly synthesized radial-only samples, not
  // coefficient vectors. They are isotropic and already include normalized
  // Y00, so convert their sampled value back to the normalized harmonic
  // coefficient before adding it to the common-frame covariance.
  if(slot<2) {
    const double y00=1.0/std::sqrt(4.0*std::numbers::pi);
    const auto it=std::find_if(basis.globalLabel[slot].begin(),basis.globalLabel[slot].end(),
        [](const Label &x){return x.m==0&&x.l==0;});
    if(it==basis.globalLabel[slot].end()) throw std::runtime_error("global scalar basis has no l=0 mode");
    const Eigen::Index col=it-basis.globalLabel[slot].begin();
    for(std::size_t k=0;k<td.records.size();++k)
      c(k,col)+=td.poles[k].at(radialSet)->values[slot][0]/y00;
  }
  return c;
}

std::array<Eigen::VectorXcd,6> makeKernels(const std::array<PairOutput,36> &pairs,
                                           bool direct) {
  const auto get=[&](int a,int f)->const Eigen::VectorXcd& {
    return direct?pairs[6*a+f].direct:pairs[6*a+f].pair;
  };
  const Eigen::Index n=get(0,0).size();
  std::array<Eigen::VectorXcd,6> k;
  for(auto &v:k)v=Eigen::VectorXcd::Zero(n);
  for(Eigen::Index p=0;p<n;++p) {
    k[0][p]=-std::real(get(1,1)[p]);
    k[1][p]=-std::real(get(0,0)[p]);
    k[2][p]=-std::real(get(0,1)[p]+get(1,0)[p]);
    k[3][p]=-2*std::real(get(2,2)[p]+get(3,3)[p]);
    k[4][p]=-std::real(0.5*(get(4,4)[p]+get(5,5)[p])-get(1,1)[p]);
    k[5][p]=-std::real(get(0,0)[p]+get(0,1)[p]+get(1,0)[p]+get(1,1)[p]);
  }
  return k;
}

std::array<Eigen::VectorXcd,6> gauntKernels(const std::array<Eigen::VectorXcd,36> &pair,
                                           int lmax) {
  auto real=[&](int a,int f,double scale=1.0){return realHarmonics(pair[6*a+f],2*lmax,scale);};
  std::array<Eigen::VectorXcd,6> k;
  k[0]=-real(1,1); k[1]=-real(0,0); k[2]=-real(0,1)-real(1,0);
  k[3]=-2.0*(real(2,2)+real(3,3));
  k[4]=-(0.5*(real(4,4)+real(5,5))-real(1,1));
  k[5]=-(real(0,0)+real(0,1)+real(1,0)+real(1,1));
  return k;
}

} // namespace

int main(int argc,char **argv) try {
  if(argc<3||argc>4) throw std::runtime_error("usage: spectral_driver SOLVER_PACKET LMAX [OBJECTIVE]");
  const auto packet=readPacket(argv[1]); const int lmax=std::stoi(argv[2]);
  const std::uint32_t objective=argc==4?static_cast<std::uint32_t>(std::stoi(argv[3])):0;
  if(lmax!=3&&lmax!=6&&lmax!=12) throw std::runtime_error("LMAX must be 3, 6, or 12");
  if(objective>=packet.objectives) throw std::runtime_error("objective index outside packet");
  auto traces=makeTraces(packet,objective);
  if(traces.empty()) throw std::runtime_error("packet contains no radial traces");
  // The three-trace packet shares angles and frames; construct geometric
  // bases and Gaunt tables once, then evaluate every exported radial trace.
  const auto &basisTrace=traces.front();
  const auto receiverIds=adjointSets(*basisTrace.records.front());
  if(receiverIds.empty()) throw std::runtime_error("packet has no receiver adjoint coefficients");
  double compactBasisSeconds=0,globalFrameSeconds=0,gauntPlanSeconds=0,fftGridPlanSeconds=0;
  const Pole *sourcePole=basisTrace.poles.front().at(0);
  PoleBasis sourceBasis=makePoleBasis(*sourcePole,basisTrace.angles,lmax,compactBasisSeconds,globalFrameSeconds);
  PoleBasis radialBasis=makePoleBasis(*basisTrace.poles.front().at(radialForwardSet),basisTrace.angles,lmax,compactBasisSeconds,globalFrameSeconds);
  // The packet's reserved radial records may already carry the scalar Y00
  // normalization. Keep their l=0 basis as unity for this exported convention.
  for(int slot=0;slot<2;++slot) radialBasis.local[slot].setOnes();
  std::map<std::uint64_t,PoleBasis> receiverBasis;
  for(auto id:receiverIds) receiverBasis.emplace(id,makePoleBasis(*basisTrace.poles.front().at(id),basisTrace.angles,lmax,compactBasisSeconds,globalFrameSeconds));
  std::map<int,std::vector<GauntTerm>> gauntTerms;
  const auto gauntStart=Clock::now();
  for(int slot=0;slot<6;++slot) {
    const int spin=slotSpin[slot];
    if(!gauntTerms.count(spin))
      gauntTerms.emplace(spin,makeGauntTerms(labelsForSpin(lmax,spin),labelsForSpin(lmax,spin),spin,lmax));
  }
  gauntPlanSeconds=std::chrono::duration<double>(Clock::now()-gauntStart).count();
  const auto fftGridStart=Clock::now();
  std::set<double> thetaSet,phiSet;
  for(const auto &p:basisTrace.angles){thetaSet.insert(p.theta);phiSet.insert(p.phi);}
  std::vector<double> theta(thetaSet.begin(),thetaSet.end()),phi(phiSet.begin(),phiSet.end());
  std::map<std::pair<std::uint64_t,std::uint64_t>,std::pair<int,int>> gridIndex;
  for(int i=0;i<(int)theta.size();++i)for(int j=0;j<(int)phi.size();++j)
    gridIndex[{bits(theta[i]),bits(phi[j])}]={i,j};
  const bool fftAvailable=theta.size()*phi.size()==basisTrace.angles.size()&&
                          theta.size()>=static_cast<std::size_t>(2*lmax+1)&&
                          phi.size()>=static_cast<std::size_t>(4*lmax+1);
  fftGridPlanSeconds=std::chrono::duration<double>(Clock::now()-fftGridStart).count();
  const double planSeconds=compactBasisSeconds+globalFrameSeconds+gauntPlanSeconds+fftGridPlanSeconds;
  std::uint64_t compactBasisPayload=compactBasisBytes(sourceBasis)+compactBasisBytes(radialBasis);
  std::uint64_t globalFramePayload=globalBasisBytes(sourceBasis)+globalBasisBytes(radialBasis);
  for(const auto &[id,b]:receiverBasis) { (void)id; compactBasisPayload+=compactBasisBytes(b); globalFramePayload+=globalBasisBytes(b); }
  std::uint64_t gauntPayload=0;
  for(const auto &[spin,terms]:gauntTerms) { (void)spin; gauntPayload+=static_cast<std::uint64_t>(terms.capacity())*sizeof(GauntTerm); }
  const std::uint64_t fftGridPayload=static_cast<std::uint64_t>(basisTrace.angles.capacity())*sizeof(Point)+
      static_cast<std::uint64_t>(theta.capacity()+phi.capacity())*sizeof(double)+
      static_cast<std::uint64_t>(gridIndex.size())*(sizeof(std::pair<const std::pair<std::uint64_t,std::uint64_t>,std::pair<int,int>>)+48);
  double aggregateFieldL2=0,aggregateFieldPoint=0,aggregateDirectL2=0,aggregateDirectPoint=0;
  double aggregateKernelL2=0,aggregateKernelPoint=0,aggregateGauntL2=0,aggregateGauntPoint=0;
  double aggregateFFTPoint=0,aggregateSynthesisSeconds=0,aggregateCovarianceFormationSeconds=0;
  double aggregateCovarianceEvaluationSeconds=0,aggregateSvdSeconds=0;
  double aggregateRecouplingSeconds=0,aggregateKernelAssemblySeconds=0,aggregateFFTProjectionSeconds=0,aggregateFFTSynthesisSeconds=0;
  std::uint64_t peakCovarianceBytes=0,peakSvdFactorBytes=0;
  Eigen::Index aggregateMinRank=std::numeric_limits<Eigen::Index>::max(),aggregateMaxRank=0;
  for(const auto &td:traces) {
    if(td.angles.size()!=basisTrace.angles.size()) throw std::runtime_error("radial traces have different angular counts");
    for(std::size_t i=0;i<td.angles.size();++i)
      if(angleKey(td.angles[i])!=angleKey(basisTrace.angles[i])) throw std::runtime_error("radial traces have different angular grids");

  double fieldL2=0,fieldPoint=0;
  std::array<double,6> fieldL2BySlot{},fieldPointBySlot{};
  std::array<double,6> forwardL2BySlot{},adjointL2BySlot{};
  double radialForwardResidualL2=0,radialAdjointResidualL2=0;
  std::array<Complex,6> firstForwardRef{},firstForwardPred{};
  std::array<Complex,3> firstAdjointScalar{};
  for(std::size_t k=0;k<td.records.size();++k) for(int slot=0;slot<6;++slot) {
    Eigen::VectorXcd fCompact=sourceBasis.local[slot]*poleCoefficients(td.poles[k].at(0),slot);
    Eigen::VectorXcd fPred=fCompact;
    Eigen::VectorXcd fRef(td.angles.size()),aRef(td.angles.size());
    Eigen::VectorXcd aPred=Eigen::VectorXcd::Zero(td.angles.size());
    for(std::size_t p=0;p<td.angles.size();++p) {
      const auto &field=td.records[k]->fields[td.pointRows[k][p]];
      fRef[p]=field[slot]; aRef[p]=field[6+slot];
    }
    Eigen::VectorXcd fRadial=Eigen::VectorXcd::Zero(td.angles.size());
    if(slot<2){fRadial=radialBasis.local[slot]*poleCoefficients(td.poles[k].at(radialForwardSet),slot);fPred+=fRadial;}
    for(auto id:receiverIds)
      aPred+=receiverBasis.at(id).local[slot]*poleCoefficients(td.poles[k].at(id),slot);
    Eigen::VectorXcd aCompact=aPred;
    if(slot<2){
      const Eigen::VectorXcd aRadial=radialBasis.local[slot]*poleCoefficients(td.poles[k].at(radialAdjointSet),slot);
      if(k==0&&slot==0) firstAdjointScalar={aRef[0],aCompact[0],aRadial[0]};
      if(slot<2) radialAdjointResidualL2=std::max(radialAdjointResidualL2,
          compare(aRef-aCompact,aRadial).l2);
      aPred+=aRadial;
    }
    const auto ef=compare(fRef,fPred),ea=compare(aRef,aPred);
    if(slot<2) radialForwardResidualL2=std::max(radialForwardResidualL2,
        compare(fRef-fCompact,fRadial).l2);
    if(k==0){firstForwardRef[slot]=fRef[0];firstForwardPred[slot]=fPred[0];}
    fieldL2=std::max({fieldL2,ef.l2,ea.l2}); fieldPoint=std::max({fieldPoint,ef.scaled,ea.scaled});
    fieldL2BySlot[slot]=std::max({fieldL2BySlot[slot],ef.l2,ea.l2});
    fieldPointBySlot[slot]=std::max({fieldPointBySlot[slot],ef.scaled,ea.scaled});
    forwardL2BySlot[slot]=std::max(forwardL2BySlot[slot],ef.l2);
    adjointL2BySlot[slot]=std::max(adjointL2BySlot[slot],ea.l2);
  }

  std::array<PairOutput,36> pairs;
  double svdSeconds=0;
  double directSynthesisSeconds=0,covarianceFormationSeconds=0,covarianceEvaluationSeconds=0;
  double directSynthesisL2=0,directSynthesisPoint=0;
  std::uint64_t covarianceMatrixBytes=0,maxCovarianceMatrixBytes=0,maxSvdFactorBytes=0;
  Eigen::Index minRank=std::numeric_limits<Eigen::Index>::max(),maxRank=0;
  for(int aslot=0;aslot<6;++aslot) for(int fslot=0;fslot<6;++fslot) {
    if(slotSpin[aslot]!=slotSpin[fslot]) continue;
    std::vector<Eigen::MatrixXcd> ba;
    for(auto id:receiverIds) {
      auto b=receiverBasis.at(id).local[aslot];
      if(aslot<2){b.conservativeResize(Eigen::NoChange,b.cols()+1);b.col(b.cols()-1)=radialBasis.local[aslot].col(0);}
      ba.push_back(std::move(b));
    }
    auto bf=sourceBasis.local[fslot];
    if(fslot<2){bf.conservativeResize(Eigen::NoChange,bf.cols()+1);bf.col(bf.cols()-1)=radialBasis.local[fslot].col(0);}
    pairs[6*aslot+fslot]=evaluatePair(td,objective,aslot,fslot,ba,
        bf,receiverIds);
    directSynthesisSeconds+=pairs[6*aslot+fslot].directSynthesisSeconds;
    covarianceFormationSeconds+=pairs[6*aslot+fslot].covarianceFormationSeconds;
    covarianceEvaluationSeconds+=pairs[6*aslot+fslot].covarianceEvaluationSeconds;
  }
  const auto directK=makeKernels(pairs,true),covK=makeKernels(pairs,false);
  double pairL2=0,pairPoint=0,kernelL2=0,kernelPoint=0;
  for(int a=0;a<6;++a) for(int f=0;f<6;++f) if(slotSpin[a]==slotSpin[f]) {
    const auto e=compare(pairs[6*a+f].direct,pairs[6*a+f].pair);
    pairL2=std::max(pairL2,e.l2); pairPoint=std::max(pairPoint,e.scaled);
    const auto ed=compare(pairs[6*a+f].direct,pairs[6*a+f].directSynthesized);
    directSynthesisL2=std::max(directSynthesisL2,ed.l2);directSynthesisPoint=std::max(directSynthesisPoint,ed.scaled);
    for(std::size_t ri=0;ri<pairs[6*a+f].receiverPair.size();++ri) {
      const auto &receiverPair=pairs[6*a+f].receiverPair[ri];
      const auto t0=Clock::now();
      const auto svd=dspectral::factorCovariance(pairs[6*a+f].receiverCovariances[ri]);
      const auto matrixBytes=static_cast<std::uint64_t>(pairs[6*a+f].receiverCovariances[ri].size())*sizeof(Complex);
      covarianceMatrixBytes+=matrixBytes;maxCovarianceMatrixBytes=std::max(maxCovarianceMatrixBytes,matrixBytes);
      maxSvdFactorBytes=std::max(maxSvdFactorBytes,static_cast<std::uint64_t>(svd.factorBytes));
      const auto rank=dspectral::smallestVerifiedRank(
          pairs[6*a+f].basisAdjoint[ri],pairs[6*a+f].basisForward,svd,receiverPair);
      minRank=std::min(minRank,rank.first); maxRank=std::max(maxRank,rank.first);
      svdSeconds+=std::chrono::duration<double>(Clock::now()-t0).count();
    }
  }
  for(int j=0;j<6;++j) {
    const auto e=compare(directK[j],covK[j]);
    kernelL2=std::max(kernelL2,e.l2); kernelPoint=std::max(kernelPoint,e.scaled);
  }

  // Build the compact-to-global spin-harmonic coefficient rotations and
  // convention-safe quadrature Gaunt tables, then compare to FFT/latitude
  // projection on the actual tensor grid.
  std::array<Eigen::VectorXcd,36> gauntPairs;
  double recouplingSeconds=0,kernelAssemblySeconds=0,fftProjectionSeconds=0,fftSynthesisSeconds=0;
  const auto recoupleStart=Clock::now();
  std::array<Eigen::MatrixXcd,6> sourceGlobal;
  std::map<std::uint64_t,std::array<Eigen::MatrixXcd,6>> receiverGlobal;
  for(int slot=0;slot<6;++slot)
    sourceGlobal[slot]=globalFrequencyCoefficients(td,0,slot,sourceBasis,radialForwardSet);
  for(auto id:receiverIds) for(int slot=0;slot<6;++slot)
    receiverGlobal[id][slot]=globalFrequencyCoefficients(td,id,slot,receiverBasis.at(id),radialAdjointSet);
  for(int a=0;a<6;++a) for(int f=0;f<6;++f) if(slotSpin[a]==slotSpin[f]) {
    const int spin=slotSpin[a];
    const auto labels=labelsForSpin(lmax,spin);
    Eigen::MatrixXcd c=Eigen::MatrixXcd::Zero(labels.size(),labels.size());
    for(auto id:receiverIds)c+=receiverGlobal.at(id)[a].adjoint()*sourceGlobal[f];
    gauntPairs[6*a+f]=gauntRecouple(c,gauntTerms.at(spin),lmax);
  }
  recouplingSeconds=std::chrono::duration<double>(Clock::now()-recoupleStart).count();
  const auto kernelStart=Clock::now();
  const auto gauntK=gauntKernels(gauntPairs,lmax);
  kernelAssemblySeconds=std::chrono::duration<double>(Clock::now()-kernelStart).count();
  double fftKernelPoint=0,gauntFftL2=0,gauntFftPoint=0,fftVsDirectHarmonicPoint=0;
  std::array<FFTProjection,6> projections;
  std::array<Eigen::VectorXcd,6> candidateGrid;
  if(fftAvailable) for(int k=0;k<6;++k) {
    Eigen::MatrixXd grid(theta.size(),phi.size());
    for(std::size_t p=0;p<td.angles.size();++p) {
      const auto ij=gridIndex.at(angleKey(td.angles[p])); grid(ij.first,ij.second)=std::real(covK[k][p]);
    }
    const auto fftProjStart=Clock::now();
    projections[k]=fftProject(grid,theta,phi,lmax);
    fftProjectionSeconds+=std::chrono::duration<double>(Clock::now()-fftProjStart).count();
    fftKernelPoint=std::max(fftKernelPoint,projections[k].error.scaled);
    const auto e=compare(projections[k].coeff,gauntK[k]);
    gauntFftL2=std::max(gauntFftL2,e.l2);gauntFftPoint=std::max(gauntFftPoint,e.scaled);
    const auto fftSynthStart=Clock::now();
    candidateGrid[k]=fftSynthesizeHarmonics(gauntK[k],theta,phi,2*lmax);
    fftSynthesisSeconds+=std::chrono::duration<double>(Clock::now()-fftSynthStart).count();
    const Eigen::VectorXcd directHarmonicGrid=synthesizeHarmonics(gauntK[k],theta,phi,2*lmax).cast<Complex>();
    fftVsDirectHarmonicPoint=std::max(fftVsDirectHarmonicPoint,
        compare(directHarmonicGrid,candidateGrid[k]).scaled);
  }
  if(fftAvailable) for(int k=0;k<6;++k) {
    Eigen::VectorXcd covGrid=covK[k];
    const auto e=compare(covGrid,candidateGrid[k]);
    fftKernelPoint=std::max(fftKernelPoint,e.scaled);
  }
  aggregateFieldL2=std::max(aggregateFieldL2,fieldL2);aggregateFieldPoint=std::max(aggregateFieldPoint,fieldPoint);
  aggregateDirectL2=std::max(aggregateDirectL2,directSynthesisL2);aggregateDirectPoint=std::max(aggregateDirectPoint,directSynthesisPoint);
  aggregateKernelL2=std::max(aggregateKernelL2,kernelL2);aggregateKernelPoint=std::max(aggregateKernelPoint,kernelPoint);
  aggregateGauntL2=std::max(aggregateGauntL2,gauntFftL2);aggregateGauntPoint=std::max(aggregateGauntPoint,gauntFftPoint);
  aggregateFFTPoint=std::max(aggregateFFTPoint,fftKernelPoint);
  aggregateSynthesisSeconds+=directSynthesisSeconds;
  aggregateCovarianceFormationSeconds+=covarianceFormationSeconds;
  aggregateCovarianceEvaluationSeconds+=covarianceEvaluationSeconds;
  aggregateSvdSeconds+=svdSeconds;
  aggregateRecouplingSeconds+=recouplingSeconds;aggregateKernelAssemblySeconds+=kernelAssemblySeconds;
  aggregateFFTProjectionSeconds+=fftProjectionSeconds;aggregateFFTSynthesisSeconds+=fftSynthesisSeconds;
  peakCovarianceBytes=std::max(peakCovarianceBytes,covarianceMatrixBytes);
  peakSvdFactorBytes=std::max(peakSvdFactorBytes,maxSvdFactorBytes);
  aggregateMinRank=std::min(aggregateMinRank,minRank);aggregateMaxRank=std::max(aggregateMaxRank,maxRank);
  std::cout<<std::setprecision(17)
    <<"{\"degree\":"<<lmax<<",\"objective\":"<<objective
    <<",\"frequency_count\":"<<td.records.size()<<",\"trace_element\":"<<td.key.first
    <<",\"trace_radius_bits\":"<<td.key.second<<",\"angles\":"<<td.angles.size()
    <<",\"receiver_count\":"<<receiverIds.size()
    <<",\"field_reconstruction_l2\":"<<fieldL2<<",\"field_reconstruction_pointwise\":"<<fieldPoint
    <<",\"radial_forward_residual_l2\":"<<radialForwardResidualL2
    <<",\"radial_adjoint_residual_l2\":"<<radialAdjointResidualL2
    <<",\"field_l2_by_slot\":["<<fieldL2BySlot[0]<<","<<fieldL2BySlot[1]<<","<<fieldL2BySlot[2]<<","<<fieldL2BySlot[3]<<","<<fieldL2BySlot[4]<<","<<fieldL2BySlot[5]<<"]"
    <<",\"forward_l2_by_slot\":["<<forwardL2BySlot[0]<<","<<forwardL2BySlot[1]<<","<<forwardL2BySlot[2]<<","<<forwardL2BySlot[3]<<","<<forwardL2BySlot[4]<<","<<forwardL2BySlot[5]<<"]"
    <<",\"adjoint_l2_by_slot\":["<<adjointL2BySlot[0]<<","<<adjointL2BySlot[1]<<","<<adjointL2BySlot[2]<<","<<adjointL2BySlot[3]<<","<<adjointL2BySlot[4]<<","<<adjointL2BySlot[5]<<"]"
    <<",\"field_pointwise_by_slot\":["<<fieldPointBySlot[0]<<","<<fieldPointBySlot[1]<<","<<fieldPointBySlot[2]<<","<<fieldPointBySlot[3]<<","<<fieldPointBySlot[4]<<","<<fieldPointBySlot[5]<<"]"
    <<",\"first_forward_ref\":["<<firstForwardRef[0].real()<<","<<firstForwardRef[0].imag()<<","<<firstForwardRef[1].real()<<","<<firstForwardRef[1].imag()<<"]"
    <<",\"first_forward_pred\":["<<firstForwardPred[0].real()<<","<<firstForwardPred[0].imag()<<","<<firstForwardPred[1].real()<<","<<firstForwardPred[1].imag()<<"]"
    <<",\"first_adjoint_scalar_ref_compact_radial\":["<<firstAdjointScalar[0].real()<<","<<firstAdjointScalar[0].imag()<<","<<firstAdjointScalar[1].real()<<","<<firstAdjointScalar[1].imag()<<","<<firstAdjointScalar[2].real()<<","<<firstAdjointScalar[2].imag()<<"]"
    <<",\"direct_synthesis_l2\":"<<directSynthesisL2<<",\"direct_synthesis_pointwise\":"<<directSynthesisPoint
    <<",\"pair_covariance_l2\":"<<pairL2<<",\"pair_covariance_pointwise\":"<<pairPoint
    <<",\"kernel_l2\":"<<kernelL2<<",\"kernel_pointwise\":"<<kernelPoint
    <<",\"fft_reconstruction_pointwise\":"<<fftKernelPoint
    <<",\"gaunt_fft_l2\":"<<gauntFftL2<<",\"gaunt_fft_pointwise\":"<<gauntFftPoint
    <<",\"fft_vs_direct_harmonic_pointwise\":"<<fftVsDirectHarmonicPoint
    <<",\"fft_grid_sufficient\":"<<(fftAvailable?"true":"false")
    <<",\"verified_svd_rank_min\":"<<minRank<<",\"verified_svd_rank_max\":"<<maxRank
    <<",\"direct_synthesis_seconds\":"<<directSynthesisSeconds
    <<",\"covariance_formation_seconds\":"<<covarianceFormationSeconds
    <<",\"covariance_evaluation_seconds\":"<<covarianceEvaluationSeconds<<",\"svd_seconds\":"<<svdSeconds
    <<",\"recoupling_seconds\":"<<recouplingSeconds<<",\"kernel_assembly_seconds\":"<<kernelAssemblySeconds
    <<",\"fft_projection_seconds\":"<<fftProjectionSeconds<<",\"fft_synthesis_seconds\":"<<fftSynthesisSeconds
    <<",\"covariance_matrix_bytes_sum\":"<<covarianceMatrixBytes
    <<",\"covariance_matrix_bytes_max\":"<<maxCovarianceMatrixBytes
    <<",\"svd_factor_bytes_max\":"<<maxSvdFactorBytes
    <<"}\n";
  }
  std::cout<<std::setprecision(17)
    <<"{\"aggregate\":true,\"degree\":"<<lmax<<",\"objective\":"<<objective
    <<",\"trace_count\":"<<traces.size()<<",\"frequency_count\":"<<basisTrace.records.size()
    <<",\"global_plan_prepare_seconds\":"<<planSeconds
    <<",\"compact_basis_prepare_seconds\":"<<compactBasisSeconds
    <<",\"global_frame_basis_prepare_seconds\":"<<globalFrameSeconds
    <<",\"gaunt_plan_prepare_seconds\":"<<gauntPlanSeconds
    <<",\"fft_grid_plan_prepare_seconds\":"<<fftGridPlanSeconds
    <<",\"compact_basis_payload_bytes\":"<<compactBasisPayload
    <<",\"global_frame_basis_payload_bytes\":"<<globalFramePayload
    <<",\"gaunt_plan_payload_bytes\":"<<gauntPayload
    <<",\"fft_grid_plan_payload_bytes_rough\":"<<fftGridPayload
    <<",\"total_plan_payload_bytes_rough\":"<<(compactBasisPayload+globalFramePayload+gauntPayload+fftGridPayload)
    <<",\"field_reconstruction_l2_max\":"<<aggregateFieldL2<<",\"field_reconstruction_pointwise_max\":"<<aggregateFieldPoint
    <<",\"direct_synthesis_l2_max\":"<<aggregateDirectL2<<",\"direct_synthesis_pointwise_max\":"<<aggregateDirectPoint
    <<",\"frequency_covariance_kernel_l2_max\":"<<aggregateKernelL2
    <<",\"frequency_covariance_kernel_pointwise_max\":"<<aggregateKernelPoint
    <<",\"gaunt_fft_l2_max\":"<<aggregateGauntL2<<",\"gaunt_fft_pointwise_max\":"<<aggregateGauntPoint
    <<",\"fft_reconstruction_pointwise_max\":"<<aggregateFFTPoint
    <<",\"direct_synthesis_seconds_sum\":"<<aggregateSynthesisSeconds
    <<",\"covariance_formation_seconds_sum\":"<<aggregateCovarianceFormationSeconds
    <<",\"covariance_evaluation_seconds_sum\":"<<aggregateCovarianceEvaluationSeconds
    <<",\"svd_seconds_sum\":"<<aggregateSvdSeconds
    <<",\"recoupling_seconds_sum\":"<<aggregateRecouplingSeconds
    <<",\"kernel_assembly_seconds_sum\":"<<aggregateKernelAssemblySeconds
    <<",\"fft_projection_seconds_sum\":"<<aggregateFFTProjectionSeconds
    <<",\"fft_synthesis_seconds_sum\":"<<aggregateFFTSynthesisSeconds
    <<",\"covariance_matrix_bytes_peak_trace\":"<<peakCovarianceBytes
    <<",\"svd_factor_bytes_peak\":"<<peakSvdFactorBytes
    <<",\"verified_svd_rank_min\":"<<aggregateMinRank<<",\"verified_svd_rank_max\":"<<aggregateMaxRank
    <<"}\n";
  return 0;
} catch(const std::exception &e) { std::cerr<<"spectral prototype: "<<e.what()<<'\n';return 1; }

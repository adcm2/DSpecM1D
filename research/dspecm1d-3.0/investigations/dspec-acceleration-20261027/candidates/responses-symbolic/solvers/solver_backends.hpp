#pragma once
#include "../profile.hpp"

// Diagnostic solver alternatives.  These classes are deliberately isolated
// from the accepted production headers and expose the same matrix/RHS shape as
// dsfast::SharedSymmetricFactor.
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/LU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace dssolver {
using Complex = std::complex<double>;
using Sparse = Eigen::SparseMatrix<Complex, Eigen::ColMajor, int>;
using Dense = Eigen::MatrixXcd;
using RealVector = Eigen::VectorXd;

// Fortran LAPACK symbols supplied by the audited OpenBLAS/LAPACK link.
extern "C" {
void zgbtrf_(const int*,const int*,const int*,const int*,Complex*,const int*,int*,int*);
void zgbtrs_(const char*,const int*,const int*,const int*,const int*,const Complex*,
             const int*,const int*,Complex*,const int*,int*);
}
static_assert(sizeof(Complex)==2*sizeof(double),"complex ABI is incompatible with ZGBTRF");

struct Bandwidth {
  int lower = 0;
  int upper = 0;
};

inline Bandwidth measureBandwidth(const Sparse &matrix) {
  if (matrix.rows() <= 0 || matrix.rows() != matrix.cols())
    throw std::invalid_argument("bandwidth requires a nonempty square matrix");
  Bandwidth bw;
  for (int col = 0; col < matrix.outerSize(); ++col)
    for (Sparse::InnerIterator it(matrix, col); it; ++it) {
      if (!std::isfinite(std::abs(it.value())))
        throw std::invalid_argument("matrix contains nonfinite values");
      bw.lower = std::max(bw.lower, static_cast<int>(it.row() - it.col()));
      bw.upper = std::max(bw.upper, static_cast<int>(it.col() - it.row()));
    }
  return bw;
}

inline bool samePattern(const Sparse &a, const Sparse &b) {
  if (a.rows() != b.rows() || a.cols() != b.cols() ||
      a.nonZeros() != b.nonZeros()) return false;
  for (int col = 0; col <= a.outerSize(); ++col)
    if (a.outerIndexPtr()[col] != b.outerIndexPtr()[col]) return false;
  for (int k = 0; k < a.nonZeros(); ++k)
    if (a.innerIndexPtr()[k] != b.innerIndexPtr()[k]) return false;
  return true;
}

// Shared positive real diagonal equilibration: A = D H D, rhs_A = D rhs,
// x = D y. The same D is correct for A^H because D is real and symmetric.
inline std::pair<Sparse, RealVector> equilibrate(const Sparse &matrix,
                                                  double *elapsedSeconds=nullptr) {
  const auto start=std::chrono::steady_clock::now();
  if (matrix.rows() <= 0 || matrix.rows() != matrix.cols())
    throw std::invalid_argument("equilibration requires a nonempty square matrix");
  RealVector d = RealVector::Zero(matrix.rows());
  for (int col = 0; col < matrix.outerSize(); ++col)
    for (Sparse::InnerIterator it(matrix, col); it; ++it) {
      const double a = std::abs(it.value());
      if (!std::isfinite(a)) throw std::invalid_argument("matrix contains nonfinite values");
      d[it.row()] = std::max(d[it.row()], a);
    }
  for (Eigen::Index i = 0; i < d.size(); ++i) {
    if (!(d[i] > 0.0)) throw std::invalid_argument("matrix contains a zero row");
    d[i] = 1.0 / std::sqrt(d[i]);
    if (!std::isfinite(d[i])) throw std::invalid_argument("invalid row scaling");
  }
  Sparse a = matrix;
  for (int col = 0; col < a.outerSize(); ++col)
    for (Sparse::InnerIterator it(a, col); it; ++it) {
      it.valueRef() *= d[it.row()] * d[it.col()];
      if (!std::isfinite(std::abs(it.valueRef())))
        throw std::invalid_argument("equilibrated matrix contains nonfinite values");
    }
  if(elapsedSeconds)*elapsedSeconds=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-start).count();
  return {std::move(a), std::move(d)};
}

inline void checkRhs(const Dense &rhs, Eigen::Index n) {
  if (rhs.rows() != n || rhs.cols() <= 0 || !rhs.allFinite())
    throw std::invalid_argument("solver RHS dimensions or values are invalid");
}

inline double relativeResidual(const Sparse &matrix,const Dense &solution,
                               const Dense &rhs) {
  if(matrix.rows()!=solution.rows()||matrix.rows()!=rhs.rows()||
     solution.cols()!=rhs.cols())
    throw std::invalid_argument("residual dimensions do not match");
  const double numerator=(matrix*solution-rhs).norm();
  const double denominator=matrix.norm()*solution.norm()+rhs.norm();
  if(denominator==0.0)return numerator==0.0?0.0:std::numeric_limits<double>::infinity();
  return numerator/denominator;
}
inline double relativeFieldDifference(const Dense &candidate,const Dense &reference) {
  if(candidate.rows()!=reference.rows()||candidate.cols()!=reference.cols())
    throw std::invalid_argument("field comparison dimensions do not match");
  const double referenceNorm=reference.norm();
  if(referenceNorm==0.0)
    return candidate.norm()==0.0?0.0:std::numeric_limits<double>::infinity();
  return (candidate-reference).norm()/referenceNorm;
}

/** LAPACK general-complex band LU, with the LP64 Fortran ABI used by the
 * linked OpenBLAS/LAPACK in the audit build. ABI and numerical validation must
 * be completed in the parent-assigned compiler slot before use. */
class EquilibratedGeneralBandLU {
public:
  struct Stats {
    int n=0, kl=0, ku=0, ldab=0; std::size_t packedBytes=0;
    double scalingSeconds=0.0, packSeconds=0.0, analysisSeconds=0.0,
           factorSeconds=0.0, substitutionSeconds=0.0, recoverySeconds=0.0;
  };

  void compute(const Sparse &matrix) { factorize(matrix); }
  void factorize(const Sparse &matrix) {
    ready_ = false;
    stats_={};
    auto scaled = equilibrate(matrix,&stats_.scalingSeconds);
    a_ = std::move(scaled.first);
    d_ = std::move(scaled.second);
    const auto packStart=std::chrono::steady_clock::now();
    const auto bw = measureBandwidth(a_);
    const int n = static_cast<int>(a_.rows());
    const int kl = bw.lower, ku = bw.upper;
    const int ldab = 2 * kl + ku + 1;
    if (ldab < kl + ku + 1) throw std::logic_error("invalid LAPACK band leading dimension");
    ab_.assign(static_cast<std::size_t>(ldab) * n, Complex{});
    for (int col=0; col<a_.outerSize(); ++col)
      for (Sparse::InnerIterator it(a_,col); it; ++it) {
        const int row=static_cast<int>(it.row());
        const int br=kl+ku+row-col;
        if (br < kl || br > 2*kl+ku)
          throw std::logic_error("measured bandwidth failed to contain a matrix entry");
        ab_[static_cast<std::size_t>(col)*ldab+br]=it.value();
      }
    ipiv_.resize(n);
    stats_.packSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-packStart).count();
    int info=0;
    const auto factorStart=std::chrono::steady_clock::now();
    zgbtrf_(&n,&n,&kl,&ku,ab_.data(),&ldab,ipiv_.data(),&info);
    stats_.factorSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-factorStart).count();
    if (info < 0) throw std::runtime_error("ZGBTRF rejected an argument");
    if (info > 0) throw std::runtime_error("ZGBTRF found a singular band factor");
    stats_.n=n;stats_.kl=kl;stats_.ku=ku;stats_.ldab=ldab;
    stats_.packedBytes=ab_.size()*sizeof(Complex);
    ready_=true;
  }
  Dense solve(const Dense &rhs) const { return solveImpl(rhs,'N'); }
  Dense solveAdjoint(const Dense &rhs) const { return solveImpl(rhs,'C'); }
  const Stats &stats() const { return stats_; }

private:
  Dense solveImpl(const Dense &rhs, char trans) const {
    if (!ready_) throw std::logic_error("band factor is not ready");
    checkRhs(rhs,d_.size());
    Dense x=rhs;
    for (Eigen::Index i=0;i<x.rows();++i) x.row(i)*=d_[i];
    const int n=stats_.n, kl=stats_.kl, ku=stats_.ku, nrhs=static_cast<int>(x.cols());
    const int ldb=std::max(1,n);
    int info=0; const auto solveStart=std::chrono::steady_clock::now();
    zgbtrs_(&trans,&n,&kl,&ku,&nrhs,ab_.data(),&stats_.ldab,ipiv_.data(),
            x.data(),&ldb,&info);
    stats_.substitutionSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-solveStart).count();
    if (info != 0 || !x.allFinite()) throw std::runtime_error("ZGBTRS solve failed");
    for (Eigen::Index i=0;i<x.rows();++i) x.row(i)*=d_[i];
    if (!x.allFinite()) throw std::runtime_error("band solve returned nonfinite values");
    return x;
  }
  Sparse a_;
  RealVector d_;
  std::vector<Complex> ab_;
  std::vector<int> ipiv_;
  mutable Stats stats_;
  bool ready_=false;
};

/** Natural physical ordering control with verified pattern-only symbolic
 * reuse. Numeric factors and diagonal equilibration are recomputed each call. */
class EquilibratedNaturalSparseLU {
public:
  struct Stats {
    bool symbolicReused=false; int n=0; Eigen::Index nnz=0;
    double scalingSeconds=0.0, analysisSeconds=0.0, factorSeconds=0.0,
           substitutionSeconds=0.0;
  };
  using Solver=Eigen::SparseLU<Sparse,Eigen::NaturalOrdering<int>>;
  void compute(const Sparse &matrix) {
    ready_=false;
    stats_={}; auto scaled=equilibrate(matrix,&stats_.scalingSeconds);
    a_=std::move(scaled.first); d_=std::move(scaled.second);
    auto tick=std::chrono::steady_clock::now();
    solver_.analyzePattern(a_);
    stats_.analysisSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
    tick=std::chrono::steady_clock::now();
    solver_.factorize(a_);
    stats_.factorSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
    checkFactor("NaturalOrdering SparseLU compute");
    previousPattern_=a_; ready_=true;
    stats_.symbolicReused=false;stats_.n=static_cast<int>(a_.rows());stats_.nnz=a_.nonZeros();
  }
  void factorize(const Sparse &matrix) {
    stats_.scalingSeconds=0.0;stats_.analysisSeconds=0.0;
    auto scaled=equilibrate(matrix,&stats_.scalingSeconds);
    Sparse next=std::move(scaled.first);
    const bool reuse=ready_ && samePattern(previousPattern_,next);
    d_=std::move(scaled.second); a_=std::move(next);
    auto tick=std::chrono::steady_clock::now();
    if (reuse) solver_.factorize(a_);
    else { solver_.analyzePattern(a_); stats_.analysisSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count(); tick=std::chrono::steady_clock::now(); solver_.factorize(a_); }
    stats_.factorSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
    checkFactor("NaturalOrdering SparseLU factorize");
    previousPattern_=a_; ready_=true;
    stats_.symbolicReused=reuse;stats_.n=static_cast<int>(a_.rows());stats_.nnz=a_.nonZeros();
  }
  Dense solve(const Dense &rhs) { return solveImpl(rhs,false); }
  Dense solveAdjoint(const Dense &rhs) { return solveImpl(rhs,true); }
  const Stats &stats() const { return stats_; }
private:
  void checkFactor(const char *where) {
    if (solver_.info()!=Eigen::Success) throw std::runtime_error(where);
  }
  Dense solveImpl(const Dense &rhs,bool adjoint) {
    if (!ready_) throw std::logic_error("natural sparse factor is not ready");
    checkRhs(rhs,d_.size()); Dense b=rhs;
    for(Eigen::Index i=0;i<b.rows();++i)b.row(i)*=d_[i];
    const auto solveStart=std::chrono::steady_clock::now();
    Dense x;
    if(adjoint)x=solver_.adjoint().solve(b);
    else x=solver_.solve(b);
    stats_.substitutionSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-solveStart).count();
    if(solver_.info()!=Eigen::Success||!x.allFinite())
      throw std::runtime_error("NaturalOrdering SparseLU solve failed");
    for(Eigen::Index i=0;i<x.rows();++i)x.row(i)*=d_[i];
    if(!x.allFinite())throw std::runtime_error("NaturalOrdering solution is nonfinite");
    return x;
  }
  Sparse a_,previousPattern_;
  RealVector d_;
  Solver solver_;
  mutable Stats stats_;
  bool ready_=false;
};

/** The caller supplies actual SEM local-to-global maps. Each element map must
 * list all local fields in a consistent order; its endpoint set marks every
 * constrained/interface/surface/centre trace DOF. All other local DOFs are
 * eliminated only if they are exclusive to exactly one element. */
struct ElementPartitionInput {
  std::vector<std::vector<int>> elementDofs;
  std::vector<std::vector<int>> elementBoundaryDofs;
};

struct CondensationPartition {
  std::vector<int> traceDofs;
  std::vector<std::vector<int>> interiorDofs;
  std::vector<int> traceIndex; // global -> condensed index, -1 for interior
};

// These builders deliberately derive local maps from the SEM's own LTG
// accessors. Endpoint rows stay in the Schur trace, including duplicated
// fluid-solid tangential IDs and surface/centre constraints.
template<class Sem> ElementPartitionInput radialElements(const Sem &sem) {
  const int ne=sem.mesh().NE(), nq=sem.mesh().NN();
  ElementPartitionInput out; out.elementDofs.resize(static_cast<std::size_t>(ne));
  out.elementBoundaryDofs.resize(static_cast<std::size_t>(ne));
  for(int e=0;e<ne;++e)for(int field=0;field<2;++field)for(int q=0;q<nq;++q){
    const int g=static_cast<int>(sem.ltgR(field,e,q));
    out.elementDofs[static_cast<std::size_t>(e)].push_back(g);
    if(q==0||q==nq-1)out.elementBoundaryDofs[static_cast<std::size_t>(e)].push_back(g);
  }
  return out;
}
template<class Sem> ElementPartitionInput spheroidalElements(const Sem &sem) {
  const int ne=sem.mesh().NE(), nq=sem.mesh().NN();
  ElementPartitionInput out; out.elementDofs.resize(static_cast<std::size_t>(ne));
  out.elementBoundaryDofs.resize(static_cast<std::size_t>(ne));
  for(int e=0;e<ne;++e)for(int field=0;field<3;++field)for(int q=0;q<nq;++q){
    const int g=static_cast<int>(sem.ltgS(field,e,q));
    out.elementDofs[static_cast<std::size_t>(e)].push_back(g);
    if(q==0||q==nq-1)out.elementBoundaryDofs[static_cast<std::size_t>(e)].push_back(g);
  }
  return out;
}
template<class Sem> ElementPartitionInput toroidalElements(const Sem &sem) {
  const int nq=sem.mesh().NN(), ne=sem.eu()-sem.el();
  ElementPartitionInput out; out.elementDofs.resize(static_cast<std::size_t>(ne));
  out.elementBoundaryDofs.resize(static_cast<std::size_t>(ne));
  for(int e=sem.el();e<sem.eu();++e)for(int q=0;q<nq;++q){
    const int g=static_cast<int>(sem.ltgT(e,q));
    auto &local=out.elementDofs[static_cast<std::size_t>(e-sem.el())];
    auto &boundary=out.elementBoundaryDofs[static_cast<std::size_t>(e-sem.el())];
    local.push_back(g); if(q==0||q==nq-1)boundary.push_back(g);
  }
  return out;
}

inline CondensationPartition makeCondensationPartition(
    int dimension,const ElementPartitionInput &input) {
  if(dimension<=0||input.elementDofs.empty()||
     input.elementDofs.size()!=input.elementBoundaryDofs.size())
    throw std::invalid_argument("invalid element partition input");
  std::vector<int> elementUses(static_cast<std::size_t>(dimension),0);
  std::vector<char> boundary(static_cast<std::size_t>(dimension),0);
  for(std::size_t e=0;e<input.elementDofs.size();++e){
    if(input.elementDofs[e].empty())throw std::invalid_argument("empty element map");
    for(int g:input.elementDofs[e]){
      if(g<0||g>=dimension)throw std::out_of_range("element map index outside matrix");
      ++elementUses[static_cast<std::size_t>(g)];
    }
    for(int g:input.elementBoundaryDofs[e]){
      if(std::find(input.elementDofs[e].begin(),input.elementDofs[e].end(),g)==input.elementDofs[e].end())
        throw std::invalid_argument("element boundary index is not in its local map");
      boundary[static_cast<std::size_t>(g)]=1;
    }
  }
  CondensationPartition result;
  result.traceIndex.assign(static_cast<std::size_t>(dimension),-1);
  result.interiorDofs.resize(input.elementDofs.size());
  std::vector<char> selected(static_cast<std::size_t>(dimension),0);
  for(std::size_t e=0;e<input.elementDofs.size();++e){
    for(int g:input.elementDofs[e]){
      if(!boundary[static_cast<std::size_t>(g)]&&elementUses[static_cast<std::size_t>(g)]==1){
        if(selected[static_cast<std::size_t>(g)])throw std::logic_error("interior DOF selected twice");
        selected[static_cast<std::size_t>(g)]=1;
        result.interiorDofs[e].push_back(g);
      }
    }
  }
  for(int g=0;g<dimension;++g){
    if(elementUses[static_cast<std::size_t>(g)]==0)
      throw std::invalid_argument("matrix contains a DOF absent from all element maps");
    if(!selected[static_cast<std::size_t>(g)]){
      result.traceIndex[static_cast<std::size_t>(g)]=static_cast<int>(result.traceDofs.size());
      result.traceDofs.push_back(g);
    }
  }
  if(result.traceDofs.empty())throw std::invalid_argument("condensation has no trace DOFs");
  return result;
}

/** Exact element-interior Schur elimination. Local blocks use pivoted complete
 * LU and are rejected if singular, poorly conditioned, or solve-residual poor;
 * callers should fall back to a whole-system pivoted solver on any exception. */
template<class SchurSolver=EquilibratedNaturalSparseLU>
class ExactStaticCondensation {
public:
  struct Stats {
    int dimension=0,traceSize=0,interiorSize=0,elementCount=0;
    bool symbolicReused=false;
    double minLocalRcond=std::numeric_limits<double>::infinity();
    double maxLocalSolveResidual=0.0;
    double scalingSeconds=0.0, localFactorSeconds=0.0,
           schurAssemblySeconds=0.0, analysisSeconds=0.0,
           schurFactorSeconds=0.0, substitutionSeconds=0.0,
           recoverySeconds=0.0;
  };
  explicit ExactStaticCondensation(double minimumLocalRcond=1e-12,
                                   double maximumLocalResidual=1e-11)
      :minimumLocalRcond_(minimumLocalRcond),maximumLocalResidual_(maximumLocalResidual) {
    if(!(minimumLocalRcond_>0.0)||!(maximumLocalResidual_>0.0))
      throw std::invalid_argument("invalid local conditioning gates");
  }
  void compute(const Sparse &matrix,const ElementPartitionInput &elements) {
    factorizeImpl(matrix,elements,false);
  }
  void factorize(const Sparse &matrix,const ElementPartitionInput &elements) {
    factorizeImpl(matrix,elements,true);
  }
  Dense solve(const Dense &rhs) const { return solveImpl(rhs,false); }
  Dense solveAdjoint(const Dense &rhs) const { return solveImpl(rhs,true); }
  const CondensationPartition &partition() const { return partition_; }
  const Stats &stats() const { return stats_; }

private:
  struct Local {
    std::vector<int> interior;
    std::vector<int> trace;
    Eigen::MatrixXcd aii, ait, ati;
    Eigen::FullPivLU<Eigen::MatrixXcd> lu;
    Eigen::FullPivLU<Eigen::MatrixXcd> luAdjoint;
    double rcond=1.0;
  };
  static double relativeResidual(const Eigen::MatrixXcd &a,const Dense &x,const Dense &b) {
    const double numerator=(a*x-b).norm();
    const double denominator=a.norm()*x.norm()+b.norm();
    if(denominator==0.0)return numerator==0.0?0.0:std::numeric_limits<double>::infinity();
    return numerator/denominator;
  }
  void factorizeImpl(const Sparse &matrix,const ElementPartitionInput &elements,bool allowReuse) {
    const bool wasReady=ready_;
    const Sparse priorSchur=schur_;
    ready_=false;
    stats_={}; auto scaled=equilibrate(matrix,&stats_.scalingSeconds);
    a_=std::move(scaled.first); d_=std::move(scaled.second);
    partition_=makeCondensationPartition(static_cast<int>(a_.rows()),elements);
    const auto nt=static_cast<Eigen::Index>(partition_.traceDofs.size());
    std::vector<Eigen::Triplet<Complex,int>> schurEntries;
    schurEntries.reserve(static_cast<std::size_t>(a_.nonZeros()));
    for(int col=0;col<a_.outerSize();++col)
      for(Sparse::InnerIterator it(a_,col);it;++it){
        const int ri=partition_.traceIndex[static_cast<std::size_t>(it.row())];
        const int ci=partition_.traceIndex[static_cast<std::size_t>(it.col())];
        if(ri>=0&&ci>=0)schurEntries.emplace_back(ri,ci,it.value());
      }
    locals_.clear(); locals_.reserve(elements.elementDofs.size());
    std::vector<int> interiorOwner(static_cast<std::size_t>(a_.rows()),-1);
    for(std::size_t e=0;e<partition_.interiorDofs.size();++e)
      for(int g:partition_.interiorDofs[e]) interiorOwner[static_cast<std::size_t>(g)]=static_cast<int>(e);
    for(int col=0;col<a_.outerSize();++col)
      for(Sparse::InnerIterator it(a_,col);it;++it){
        const int ri=interiorOwner[static_cast<std::size_t>(it.row())];
        const int ci=interiorOwner[static_cast<std::size_t>(it.col())];
        if(ri>=0&&ci>=0&&ri!=ci&&it.value()!=Complex{})
          throw std::invalid_argument("interior blocks couple across elements; condensation partition is invalid");
      }
    const double scalingSeconds=stats_.scalingSeconds;
    stats_={}; stats_.scalingSeconds=scalingSeconds;
    stats_.dimension=static_cast<int>(a_.rows());
    stats_.traceSize=static_cast<int>(nt); stats_.elementCount=static_cast<int>(elements.elementDofs.size());
    const auto localStart=std::chrono::steady_clock::now();
    for(std::size_t e=0;e<elements.elementDofs.size();++e){
      Local local;
      local.interior=partition_.interiorDofs[e];
      for(int g:elements.elementDofs[e])
        if(partition_.traceIndex[static_cast<std::size_t>(g)]>=0)local.trace.push_back(g);
      if(local.interior.empty()) { locals_.push_back(std::move(local)); continue; }
      const Eigen::Index ni=static_cast<Eigen::Index>(local.interior.size());
      const Eigen::Index ne=static_cast<Eigen::Index>(local.trace.size());
      local.aii.resize(ni,ni); local.ait.resize(ni,ne); local.ati.resize(ne,ni);
      for(Eigen::Index i=0;i<ni;++i){
        for(Eigen::Index j=0;j<ni;++j)local.aii(i,j)=a_.coeff(local.interior[static_cast<std::size_t>(i)],local.interior[static_cast<std::size_t>(j)]);
        for(Eigen::Index j=0;j<ne;++j)local.ait(i,j)=a_.coeff(local.interior[static_cast<std::size_t>(i)],local.trace[static_cast<std::size_t>(j)]);
      }
      for(Eigen::Index i=0;i<ne;++i)
        for(Eigen::Index j=0;j<ni;++j)local.ati(i,j)=a_.coeff(local.trace[static_cast<std::size_t>(i)],local.interior[static_cast<std::size_t>(j)]);
      local.lu.compute(local.aii);
      local.luAdjoint.compute(local.aii.adjoint());
      local.rcond=local.lu.rcond();
      stats_.minLocalRcond=std::min(stats_.minLocalRcond,local.rcond);
      if(local.lu.rank()!=ni||local.luAdjoint.rank()!=ni||
         !std::isfinite(local.rcond)||local.rcond<minimumLocalRcond_)
        throw std::runtime_error("element interior is singular or poorly conditioned; use global fallback");
      Eigen::MatrixXcd solvedAit=local.lu.solve(local.ait);
      const double res=relativeResidual(local.aii,solvedAit,local.ait);
      stats_.maxLocalSolveResidual=std::max(stats_.maxLocalSolveResidual,res);
      if(!solvedAit.allFinite()||res>maximumLocalResidual_)
        throw std::runtime_error("element interior elimination residual failed; use global fallback");
      for(Eigen::Index i=0;i<ne;++i){
        const int ti=partition_.traceIndex[static_cast<std::size_t>(local.trace[static_cast<std::size_t>(i)])];
        for(Eigen::Index j=0;j<ne;++j){
          const int tj=partition_.traceIndex[static_cast<std::size_t>(local.trace[static_cast<std::size_t>(j)])];
          const Complex correction=(local.ati.row(i)*solvedAit.col(j))(0,0);
          if(correction!=Complex{})schurEntries.emplace_back(ti,tj,-correction);
        }
      }
      stats_.interiorSize+=static_cast<int>(ni);
      locals_.push_back(std::move(local));
    }
    stats_.localFactorSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-localStart).count();
    const auto schurStart=std::chrono::steady_clock::now();
    Sparse nextSchur(nt,nt); nextSchur.setFromTriplets(schurEntries.begin(),schurEntries.end());
    nextSchur.makeCompressed();
    stats_.schurAssemblySeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-schurStart).count();
    const bool reuse=allowReuse&&wasReady&&samePattern(priorSchur,nextSchur);
    schur_=std::move(nextSchur);
    const auto factorStart=std::chrono::steady_clock::now();
    if(reuse)schurSolver_.factorize(schur_);
    else schurSolver_.compute(schur_);
    stats_.schurFactorSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-factorStart).count();
    stats_.analysisSeconds=schurSolver_.stats().analysisSeconds;
    stats_.symbolicReused=reuse;
    ready_=true;
  }
  Dense solveImpl(const Dense &rhs,bool adjoint) const {
    if(!ready_)throw std::logic_error("condensed factor is not ready");
    checkRhs(rhs,d_.size());
    Dense b=rhs; for(Eigen::Index i=0;i<b.rows();++i)b.row(i)*=d_[i];
    const Eigen::Index nt=static_cast<Eigen::Index>(partition_.traceDofs.size());
    Dense traceRhs=Eigen::MatrixXcd::Zero(nt,b.cols());
    for(Eigen::Index i=0;i<nt;++i)traceRhs.row(i)=b.row(partition_.traceDofs[static_cast<std::size_t>(i)]);
    const auto substitutionStart=std::chrono::steady_clock::now();
    std::vector<Dense> localBase(locals_.size());
    for(std::size_t e=0;e<locals_.size();++e){
      const auto &l=locals_[e]; if(l.interior.empty())continue;
      Dense bi(l.interior.size(),b.cols());
      for(std::size_t i=0;i<l.interior.size();++i)bi.row(static_cast<Eigen::Index>(i))=b.row(l.interior[i]);
      Dense yi=adjoint?l.luAdjoint.solve(bi):l.lu.solve(bi);
      if(!yi.allFinite())throw std::runtime_error("local interior substitution failed");
      localBase[e]=yi;
      if(adjoint){
        for(Eigen::Index i=0;i<static_cast<Eigen::Index>(l.trace.size());++i){
          const int ti=partition_.traceIndex[static_cast<std::size_t>(l.trace[static_cast<std::size_t>(i)])];
          traceRhs.row(ti)-=l.ait.adjoint().row(i)*yi;
        }
      }else{
        for(Eigen::Index i=0;i<static_cast<Eigen::Index>(l.trace.size());++i){
          const int ti=partition_.traceIndex[static_cast<std::size_t>(l.trace[static_cast<std::size_t>(i)])];
          traceRhs.row(ti)-=l.ati.row(i)*yi;
        }
      }
    }
    Dense xt=adjoint?schurSolver_.solveAdjoint(traceRhs):schurSolver_.solve(traceRhs);
    stats_.substitutionSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-substitutionStart).count();
    const auto recoveryStart=std::chrono::steady_clock::now();
    Dense x= Dense::Zero(b.rows(),b.cols());
    for(Eigen::Index i=0;i<nt;++i)x.row(partition_.traceDofs[static_cast<std::size_t>(i)])=xt.row(i);
    for(std::size_t e=0;e<locals_.size();++e){
      const auto &l=locals_[e]; if(l.interior.empty())continue;
      Dense traceSolution(l.trace.size(),b.cols());
      if(adjoint){
        for(Eigen::Index i=0;i<static_cast<Eigen::Index>(l.trace.size());++i){
          const int ti=partition_.traceIndex[static_cast<std::size_t>(l.trace[static_cast<std::size_t>(i)])];
          traceSolution.row(i)=xt.row(ti);
        }
        Dense correction=l.luAdjoint.solve(l.ati.adjoint()*traceSolution);
        if(!correction.allFinite())throw std::runtime_error("interior recovery failed");
        for(std::size_t i=0;i<l.interior.size();++i)x.row(l.interior[i])=(localBase[e]-correction).row(static_cast<Eigen::Index>(i));
      }else{
        for(Eigen::Index i=0;i<static_cast<Eigen::Index>(l.trace.size());++i){
          const int ti=partition_.traceIndex[static_cast<std::size_t>(l.trace[static_cast<std::size_t>(i)])];
          traceSolution.row(i)=xt.row(ti);
        }
        Dense correction=l.lu.solve(l.ait*traceSolution);
        if(!correction.allFinite())throw std::runtime_error("interior recovery failed");
        for(std::size_t i=0;i<l.interior.size();++i)x.row(l.interior[i])=(localBase[e]-correction).row(static_cast<Eigen::Index>(i));
      }
    }
    for(Eigen::Index i=0;i<x.rows();++i)x.row(i)*=d_[i];
    stats_.recoverySeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-recoveryStart).count();
    if(!x.allFinite())throw std::runtime_error("condensed solution is nonfinite");
    return x;
  }

  double minimumLocalRcond_,maximumLocalResidual_;
  Sparse a_,schur_;
  RealVector d_;
  CondensationPartition partition_;
  std::vector<Local> locals_;
  mutable SchurSolver schurSolver_;
  mutable Stats stats_;
  bool ready_=false;
};

} // namespace dssolver

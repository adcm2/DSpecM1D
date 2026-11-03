#pragma once

#include <cstddef>
#include <cstdint>
#include <array>
#include <string>

namespace YSpec::PhaseEProfile {

enum class ScopeId : std::uint8_t {
  OuterApi,
  ModelAndDepthSetup,
  ForwardFrequency,
  ForwardSynthesis,
  KernelFrequency,
  RadialSolveR,
  RadialSolveT,
  RadialSolveS,
  StrainTable,
  AngularWignerTable,
  AngularTablePacking,
  AngularEvaluation,
  MaterialContraction,
  MaterialTransform,
  Reduction,
  TimeTransform,
  DriverSetup,
  Measurement,
  OutputIo,
  Count
};

enum class WorkKind : std::uint8_t { Forward, Kernel };

#ifdef YSPEC_PHASE_E_PROFILE_DISABLED
class Scope {
 public:
  explicit Scope(ScopeId) noexcept {}
  ~Scope() noexcept = default;
  void Stop() noexcept {}
  Scope(const Scope&) = delete;
  Scope& operator=(const Scope&) = delete;
};

class WorkItem {
 public:
  WorkItem(int, WorkKind) noexcept {}
  ~WorkItem() noexcept = default;
  void Complete(bool = false) noexcept {}
  WorkItem(const WorkItem&) = delete;
  WorkItem& operator=(const WorkItem&) = delete;
};
#else
class Scope {
 public:
  explicit Scope(ScopeId id) noexcept;
  ~Scope() noexcept;
  void Stop() noexcept;
  Scope(const Scope&) = delete;
  Scope& operator=(const Scope&) = delete;

 private:
  bool enabled_ = false;
  bool finished_ = false;
  std::size_t index_ = 0;
};

class WorkItem {
 public:
  WorkItem(int bin, WorkKind kind) noexcept;
  ~WorkItem() noexcept;
  void Complete(bool skipped = false) noexcept { complete_ = true; skipped_ = skipped; }
  WorkItem(const WorkItem&) = delete;
  WorkItem& operator=(const WorkItem&) = delete;

 private:
  int bin_ = -1;
  WorkKind kind_ = WorkKind::Forward;
  std::size_t baseDepth_ = 0;
  std::array<std::uint64_t, static_cast<std::size_t>(ScopeId::Count)> baselineInclusive_{};
  std::array<std::uint64_t, static_cast<std::size_t>(ScopeId::Count)> baselineExclusive_{};
  std::array<std::uint64_t, static_cast<std::size_t>(ScopeId::Count)> baselineCalls_{};
  bool hasCallerParent_ = false;
  std::uint64_t parentChildrenBaseline_ = 0;
  int threadId_ = -1;
  bool enabled_ = false;
  bool complete_ = false;
  bool skipped_ = false;
};
#endif

void Reset() noexcept;
std::string FinalizeJson();
const char* Name(ScopeId id) noexcept;

}  // namespace YSpec::PhaseEProfile

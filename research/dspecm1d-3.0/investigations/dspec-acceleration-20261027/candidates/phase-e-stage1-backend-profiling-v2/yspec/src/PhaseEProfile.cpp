#include <YSpec/PhaseEProfile.hpp>

#include <array>
#include <chrono>
#include <exception>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace YSpec::PhaseEProfile {
namespace {

using Clock = std::chrono::steady_clock;
constexpr auto ScopeCount = static_cast<std::size_t>(ScopeId::Count);
constexpr std::size_t MaxBin = 64;

struct Entry {
  std::uint64_t inclusive = 0;
  std::uint64_t exclusive = 0;
  std::uint64_t calls = 0;
};

struct Frame {
  ScopeId id{};
  Clock::time_point start{};
  std::uint64_t children = 0;
};

struct Local {
  std::array<Entry, ScopeCount> entries{};
  std::array<Frame, 64> stack{};
  std::size_t depth = 0;
};

struct Aggregate {
  std::array<Entry, ScopeCount> workerEntries{};
  std::array<Entry, ScopeCount> callerEntries{};
  std::array<std::array<bool, MaxBin>, 2> attempted{};
  std::array<std::array<bool, MaxBin>, 2> completed{};
  std::array<bool, MaxBin> kernelSkipped{};
  std::array<std::uint64_t, 2> itemCount{};
  std::uint64_t masterThreadWorkItemNanoseconds = 0;
  std::uint64_t masterThreadWorkItemCount = 0;
  struct ThreadCoverage {
    std::array<std::array<bool, MaxBin>, 2> attempted{};
    std::array<std::array<bool, MaxBin>, 2> completed{};
    std::array<bool, MaxBin> skipped{};
    std::array<std::uint64_t, 2> itemCount{};
  };
  std::array<ThreadCoverage, 64> threadCoverage{};
};

thread_local Local local;
std::mutex aggregateMutex;
Aggregate aggregate;

constexpr std::array<const char*, ScopeCount> ScopeNames{
    "outer_api", "model_and_depth_setup", "forward_frequency",
    "forward_degree_synthesis", "kernel_frequency", "radial_solve_R", "radial_solve_T",
    "radial_solve_S", "strain_table", "angular_wigner_table",
    "angular_table_packing", "angular_evaluation", "material_contraction", "material_transform",
    "reduction", "time_transform", "driver_setup", "measurement", "output_io"};

std::size_t Index(ScopeId id) noexcept {
  return static_cast<std::size_t>(id);
}

void MergeLocal() noexcept {
  std::scoped_lock lock(aggregateMutex);
  for (std::size_t i = 0; i < ScopeCount; ++i) {
    aggregate.callerEntries[i].inclusive += local.entries[i].inclusive;
    aggregate.callerEntries[i].exclusive += local.entries[i].exclusive;
    aggregate.callerEntries[i].calls += local.entries[i].calls;
    local.entries[i] = {};
  }
}

void MergeWorkItem(
    int bin, WorkKind kind, bool complete, bool skipped,
    const std::array<std::uint64_t, ScopeCount>& baselineInclusive,
    const std::array<std::uint64_t, ScopeCount>& baselineExclusive,
    const std::array<std::uint64_t, ScopeCount>& baselineCalls,
    bool masterThreadItem, std::uint64_t masterElapsed, int threadId) noexcept {
  const auto index = static_cast<std::size_t>(kind);
  std::scoped_lock lock(aggregateMutex);
  for (std::size_t i = 0; i < ScopeCount; ++i) {
    const auto inclusive = local.entries[i].inclusive - baselineInclusive[i];
    const auto exclusive = local.entries[i].exclusive - baselineExclusive[i];
    const auto calls = local.entries[i].calls - baselineCalls[i];
    aggregate.workerEntries[i].inclusive += inclusive;
    aggregate.workerEntries[i].exclusive += exclusive;
    aggregate.workerEntries[i].calls += calls;
    local.entries[i].inclusive = baselineInclusive[i];
    local.entries[i].exclusive = baselineExclusive[i];
    local.entries[i].calls = baselineCalls[i];
  }
  ++aggregate.itemCount[index];
  if (bin >= 0 && static_cast<std::size_t>(bin) < MaxBin) {
    aggregate.attempted[index][static_cast<std::size_t>(bin)] = true;
    if (complete) aggregate.completed[index][static_cast<std::size_t>(bin)] = true;
    if (kind == WorkKind::Kernel && skipped)
      aggregate.kernelSkipped[static_cast<std::size_t>(bin)] = true;
  }
  if (masterThreadItem) {
    aggregate.masterThreadWorkItemNanoseconds += masterElapsed;
    ++aggregate.masterThreadWorkItemCount;
  }
  if (threadId >= 0 && static_cast<std::size_t>(threadId) < aggregate.threadCoverage.size()) {
    auto& coverage = aggregate.threadCoverage[static_cast<std::size_t>(threadId)];
    ++coverage.itemCount[index];
    if (bin >= 0 && static_cast<std::size_t>(bin) < MaxBin) {
      coverage.attempted[index][static_cast<std::size_t>(bin)] = true;
      if (complete) coverage.completed[index][static_cast<std::size_t>(bin)] = true;
      if (kind == WorkKind::Kernel && skipped)
        coverage.skipped[static_cast<std::size_t>(bin)] = true;
    }
  }
}

void WriteBins(std::ostringstream& out,
               const std::array<bool, MaxBin>& bins) {
  bool first = true;
  out << '[';
  for (std::size_t bin = 0; bin < MaxBin; ++bin) {
    if (!bins[bin]) continue;
    out << (first ? "" : ",") << bin;
    first = false;
  }
  out << ']';
}

}  // namespace

const char* Name(ScopeId id) noexcept {
  const auto index = Index(id);
  return index < ScopeNames.size() ? ScopeNames[index] : "invalid";
}

#ifdef YSPEC_PHASE_E_PROFILE_DISABLED

void Reset() noexcept {}
std::string FinalizeJson() { return R"({"enabled":false})"; }

#else

Scope::Scope(ScopeId id) noexcept {
  if (local.depth >= local.stack.size()) std::terminate();
  enabled_ = true;
  index_ = local.depth++;
  local.stack[index_] = {id, Clock::now(), 0};
}

Scope::~Scope() noexcept { Stop(); }

void Scope::Stop() noexcept {
  if (!enabled_ || finished_) return;
  if (local.depth != index_ + 1) std::terminate();
  const auto end = Clock::now();
  auto& frame = local.stack[index_];
  const auto elapsed = static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - frame.start)
          .count());
  --local.depth;
  auto& entry = local.entries[Index(frame.id)];
  entry.inclusive += elapsed;
  entry.exclusive += elapsed - frame.children;
  ++entry.calls;
  if (local.depth) local.stack[local.depth - 1].children += elapsed;
  finished_ = true;
}

WorkItem::WorkItem(int bin, WorkKind kind) noexcept
    : bin_{bin}, kind_{kind}, baseDepth_{local.depth},
      hasCallerParent_{local.depth > 0}, enabled_{true} {
#ifdef _OPENMP
  threadId_ = omp_get_thread_num();
#else
  threadId_ = 0;
#endif
  if (hasCallerParent_) parentChildrenBaseline_ = local.stack[local.depth - 1].children;
  for (std::size_t i = 0; i < ScopeCount; ++i) {
    baselineInclusive_[i] = local.entries[i].inclusive;
    baselineExclusive_[i] = local.entries[i].exclusive;
    baselineCalls_[i] = local.entries[i].calls;
  }
}

WorkItem::~WorkItem() noexcept {
  if (!enabled_) return;
  if (local.depth != baseDepth_) std::terminate();
  const auto masterElapsed = hasCallerParent_
      ? local.stack[baseDepth_ - 1].children - parentChildrenBaseline_
      : 0;
  MergeWorkItem(bin_, kind_, complete_, skipped_, baselineInclusive_,
                baselineExclusive_, baselineCalls_, hasCallerParent_, masterElapsed,
                threadId_);
}

void Reset() noexcept {
  if (local.depth != 0) std::terminate();
  local.entries = {};
  std::scoped_lock lock(aggregateMutex);
  aggregate = {};
}

std::string FinalizeJson() {
  if (local.depth != 0) throw std::logic_error("profile finalized with open scope");
  MergeLocal();
  std::scoped_lock lock(aggregateMutex);
  std::ostringstream out;
  out << std::setprecision(17) << "{\"enabled\":true,\"clock\":\"steady_clock elapsed per thread\",\"merge_policy\":\"thread-local records merged once per completed frequency work item and once for main thread after joins\",\"work_items\":{\"forward_count\":" << aggregate.itemCount[static_cast<std::size_t>(WorkKind::Forward)] << ",\"kernel_count\":" << aggregate.itemCount[static_cast<std::size_t>(WorkKind::Kernel)] << ",\"forward_attempted_bins\":";
  WriteBins(out, aggregate.attempted[static_cast<std::size_t>(WorkKind::Forward)]);
  out << ",\"forward_completed_bins\":";
  WriteBins(out, aggregate.completed[static_cast<std::size_t>(WorkKind::Forward)]);
  out << ",\"kernel_attempted_bins\":";
  WriteBins(out, aggregate.attempted[static_cast<std::size_t>(WorkKind::Kernel)]);
  out << ",\"kernel_completed_bins\":";
  WriteBins(out, aggregate.completed[static_cast<std::size_t>(WorkKind::Kernel)]);
  out << ",\"kernel_exact_zero_skipped_bins\":";
  WriteBins(out, aggregate.kernelSkipped);
  out << "},\"thread_coverage\":[";
  bool firstThread = true;
  for (std::size_t thread = 0; thread < aggregate.threadCoverage.size(); ++thread) {
    const auto& coverage = aggregate.threadCoverage[thread];
    if (coverage.itemCount[0] == 0 && coverage.itemCount[1] == 0) continue;
    if (!firstThread) out << ',';
    firstThread = false;
    out << "{\"thread_id\":" << thread
        << ",\"forward_work_items\":" << coverage.itemCount[0]
        << ",\"kernel_work_items\":" << coverage.itemCount[1]
        << ",\"forward_attempted_bins\":";
    WriteBins(out, coverage.attempted[0]);
    out << ",\"forward_completed_bins\":";
    WriteBins(out, coverage.completed[0]);
    out << ",\"kernel_attempted_bins\":";
    WriteBins(out, coverage.attempted[1]);
    out << ",\"kernel_completed_bins\":";
    WriteBins(out, coverage.completed[1]);
    out << ",\"kernel_exact_zero_skipped_bins\":";
    WriteBins(out, coverage.skipped);
    out << '}';
  }
  out << "],\"worker_scope_sums\":{";
  for (std::size_t i = 0; i < ScopeCount; ++i) {
    if (i) out << ',';
    const auto& entry = aggregate.workerEntries[i];
    out << '\"' << ScopeNames[i] << "\":{\"inclusive_s\":"
        << entry.inclusive * 1e-9 << ",\"exclusive_s\":"
        << entry.exclusive * 1e-9 << ",\"calls\":" << entry.calls << '}';
  }
  out << "},\"caller_thread_scope_sums\":{";
  for (std::size_t i = 0; i < ScopeCount; ++i) {
    if (i) out << ',';
    const auto& entry = aggregate.callerEntries[i];
    out << '\"' << ScopeNames[i] << "\":{\"inclusive_s\":"
        << entry.inclusive * 1e-9 << ",\"exclusive_s\":"
        << entry.exclusive * 1e-9 << ",\"calls\":" << entry.calls << '}';
  }
  out << "},\"master_thread_work_item_root_scope_bridge\":{\"elapsed_s\":"
      << aggregate.masterThreadWorkItemNanoseconds * 1e-9
      << ",\"work_item_count\":" << aggregate.masterThreadWorkItemCount
      << ",\"meaning\":\"exact direct root-scope inclusive time of frequency items executed by the OpenMP master thread; included in worker scope sums and inside the outer API caller wall scope\"}}";
  return out.str();
}

#endif

}  // namespace YSpec::PhaseEProfile

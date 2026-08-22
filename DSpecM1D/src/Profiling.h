#ifndef DSPECM1D_PROFILING_H
#define DSPECM1D_PROFILING_H

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#ifdef DSPECM1D_ENABLE_PROFILING
#include <omp.h>
#endif

namespace SPARSESPEC::detail::profiling {

enum class Mode : std::size_t { all = 0, radial, toroidal, spheroidal, count };
enum class Category : std::size_t {
  total = 0,
  sem_mesh,
  base_operator,
  truncation,
  dynamic_matrix,
  compression,
  band_pack,
  factorization,
  solve,
  projection,
  unclassified,
  count
};

inline constexpr std::size_t modeCount = static_cast<std::size_t>(Mode::count);
inline constexpr std::size_t categoryCount =
    static_cast<std::size_t>(Category::count);

inline constexpr const char *modeName(Mode mode) {
  switch (mode) {
  case Mode::all: return "all";
  case Mode::radial: return "radial";
  case Mode::toroidal: return "toroidal";
  case Mode::spheroidal: return "spheroidal";
  case Mode::count: break;
  }
  return "unknown";
}

inline constexpr const char *categoryName(Category category) {
  switch (category) {
  case Category::total: return "total_spectra";
  case Category::sem_mesh: return "adaptive_sem_mesh";
  case Category::base_operator: return "base_operator_preparation";
  case Category::truncation: return "start_truncation_extraction";
  case Category::dynamic_matrix: return "frequency_matrix_construction";
  case Category::compression: return "sparse_compression";
  case Category::band_pack: return "lapack_band_packing";
  case Category::factorization: return "factorization";
  case Category::solve: return "solve";
  case Category::projection: return "source_receiver_projection";
  case Category::unclassified: return "unclassified";
  case Category::count: break;
  }
  return "unknown";
}

struct Counts {
  std::size_t sems = 0, degrees = 0, frequencySystems = 0;
  std::size_t eigenCompute = 0, eigenFactorize = 0;
  std::size_t lapackFactorize = 0, solves = 0, bandPacks = 0, rhs = 0;
  std::size_t dimensions = 0, nonzeros = 0;
  long dimensionMin = std::numeric_limits<long>::max();
  long dimensionMax = 0;
  long nonzeroMin = std::numeric_limits<long>::max();
  long nonzeroMax = 0;
  long klMin = std::numeric_limits<long>::max();
  long klMax = 0;
  long kuMin = std::numeric_limits<long>::max();
  long kuMax = 0;
};

struct Data {
  std::array<std::array<double, categoryCount>, modeCount> seconds{};
  std::array<double, modeCount> workerTotal{};
  std::array<double, modeCount> workerCategorized{};
  std::array<double, modeCount> workerUnclassified{};
  Counts counts;
};

#ifdef DSPECM1D_ENABLE_PROFILING

class Context {
  using Clock = std::chrono::steady_clock;
  struct Slot {
    std::array<std::array<double, categoryCount>, modeCount> seconds{};
    std::array<std::array<double, categoryCount>, modeCount> workerSeconds{};
    std::array<double, modeCount> workerTotal{};
    Counts counts;
  };

public:
  Context() : slots_(static_cast<std::size_t>(std::max(1, omp_get_max_threads()))) {}

  void activate() { active_ = this; }
  void deactivate() { if (active_ == this) active_ = nullptr; }
  void setMode(Mode mode) { currentMode_ = mode; }
  static Mode mode() { return currentMode_; }

  void addTime(Category category, Mode mode, double seconds) {
    auto &current = slot();
    auto &values = current.seconds;
    values[static_cast<std::size_t>(mode)]
           [static_cast<std::size_t>(category)] += seconds;
    if (mode != Mode::all)
      values[static_cast<std::size_t>(Mode::all)]
            [static_cast<std::size_t>(category)] += seconds;
    if (omp_in_parallel()) {
      auto &worker = current.workerSeconds;
      worker[static_cast<std::size_t>(mode)]
            [static_cast<std::size_t>(category)] += seconds;
      if (mode != Mode::all)
        worker[static_cast<std::size_t>(Mode::all)]
              [static_cast<std::size_t>(category)] += seconds;
    }
  }

  void addWorkerTotal(Mode mode, double seconds) {
    auto &values = slot().workerTotal;
    values[static_cast<std::size_t>(mode)] += seconds;
    if (mode != Mode::all)
      values[static_cast<std::size_t>(Mode::all)] += seconds;
  }

  void countSem() { ++slot().counts.sems; }
  void countDegree() { ++slot().counts.degrees; }
  void countFrequencySystem(long dimension, long nonzeros, long kl, long ku) {
    auto &counts = slot().counts;
    ++counts.frequencySystems;
    ++counts.dimensions;
    ++counts.nonzeros;
    counts.dimensionMin = std::min(counts.dimensionMin, dimension);
    counts.dimensionMax = std::max(counts.dimensionMax, dimension);
    counts.nonzeroMin = std::min(counts.nonzeroMin, nonzeros);
    counts.nonzeroMax = std::max(counts.nonzeroMax, nonzeros);
    counts.klMin = std::min(counts.klMin, kl);
    counts.klMax = std::max(counts.klMax, kl);
    counts.kuMin = std::min(counts.kuMin, ku);
    counts.kuMax = std::max(counts.kuMax, ku);
  }
  void countCompute() { ++slot().counts.eigenCompute; }
  void countFactorize(bool lapack) {
    if (lapack) ++slot().counts.lapackFactorize;
    else ++slot().counts.eigenFactorize;
  }
  void countSolve(long rhs) {
    auto &counts = slot().counts;
    ++counts.solves;
    counts.rhs += static_cast<std::size_t>(std::max(0L, rhs));
  }
  void countBandPack() { ++slot().counts.bandPacks; }

  Data finish() const {
    Data result;
    for (const auto &slot : slots_) {
      for (std::size_t mode = 0; mode < modeCount; ++mode)
        for (std::size_t category = 0; category < categoryCount; ++category) {
          result.seconds[mode][category] += slot.seconds[mode][category];
          result.workerCategorized[mode] +=
              slot.workerSeconds[mode][category];
        }
      for (std::size_t mode = 0; mode < modeCount; ++mode)
        result.workerTotal[mode] += slot.workerTotal[mode];
      const auto &source = slot.counts;
      auto &target = result.counts;
      target.sems += source.sems;
      target.degrees += source.degrees;
      target.frequencySystems += source.frequencySystems;
      target.eigenCompute += source.eigenCompute;
      target.eigenFactorize += source.eigenFactorize;
      target.lapackFactorize += source.lapackFactorize;
      target.solves += source.solves;
      target.bandPacks += source.bandPacks;
      target.rhs += source.rhs;
      target.dimensions += source.dimensions;
      target.nonzeros += source.nonzeros;
      target.dimensionMin = std::min(target.dimensionMin, source.dimensionMin);
      target.dimensionMax = std::max(target.dimensionMax, source.dimensionMax);
      target.nonzeroMin = std::min(target.nonzeroMin, source.nonzeroMin);
      target.nonzeroMax = std::max(target.nonzeroMax, source.nonzeroMax);
      target.klMin = std::min(target.klMin, source.klMin);
      target.klMax = std::max(target.klMax, source.klMax);
      target.kuMin = std::min(target.kuMin, source.kuMin);
      target.kuMax = std::max(target.kuMax, source.kuMax);
    }
    for (std::size_t mode = 0; mode < modeCount; ++mode) {
      result.workerUnclassified[mode] = std::max(
          0.0, result.workerTotal[mode] - result.workerCategorized[mode]);
      result.seconds[mode][static_cast<std::size_t>(Category::unclassified)] +=
          result.workerUnclassified[mode];
    }
    return result;
  }

  static Context *active() { return active_; }

private:
  Slot &slot() {
    const int index = omp_in_parallel() ? omp_get_thread_num() : 0;
    return slots_[static_cast<std::size_t>(index) % slots_.size()];
  }

  std::vector<Slot> slots_;
  static inline Context *active_ = nullptr;
  static inline thread_local Mode currentMode_ = Mode::all;
};

class Scope {
  using Clock = std::chrono::steady_clock;
public:
  Scope(Context &context, Category category, Mode mode = Mode::all)
      : context_(&context), category_(category), mode_(mode), start_(Clock::now()) {}
  Scope(Context *context, Category category, Mode mode = Mode::all)
      : context_(context), category_(category), mode_(mode), start_(Clock::now()) {}
  ~Scope() {
    if (context_)
      context_->addTime(category_, mode_,
                        std::chrono::duration<double>(Clock::now() - start_).count());
  }
private:
  Context *context_;
  Category category_;
  Mode mode_;
  Clock::time_point start_;
};

class WorkerScope {
  using Clock = std::chrono::steady_clock;
public:
  WorkerScope(Context &context, Mode mode)
      : context_(&context), mode_(mode), start_(Clock::now()) {}
  ~WorkerScope() {
    context_->addWorkerTotal(
        mode_, std::chrono::duration<double>(Clock::now() - start_).count());
  }
private:
  Context *context_;
  Mode mode_;
  Clock::time_point start_;
};

inline thread_local Data lastData;
inline void publish(const Data &data) { lastData = data; }
inline const Data &last() { return lastData; }

#else

class Context {
public:
  void activate() {}
  void deactivate() {}
  void setMode(Mode) {}
  void addTime(Category, Mode, double) {}
  void addWorkerTotal(Mode, double) {}
  void countSem() {}
  void countDegree() {}
  void countFrequencySystem(long, long, long, long) {}
  void countCompute() {}
  void countFactorize(bool) {}
  void countSolve(long) {}
  void countBandPack() {}
  Data finish() const { return {}; }
  static Context *active() { return nullptr; }
  static Mode mode() { return Mode::all; }
};

class Scope {
public:
  Scope(Context &, Category, Mode = Mode::all) {}
  Scope(Context *, Category, Mode = Mode::all) {}
};
class WorkerScope {
public:
  WorkerScope(Context &, Mode) {}
};
inline void publish(const Data &) {}
inline const Data &last() { static const Data empty{}; return empty; }

#endif

} // namespace SPARSESPEC::detail::profiling

#endif

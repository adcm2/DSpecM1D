#pragma once
#include <cstddef>
#include <exception>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <omp.h>

namespace dspec::frequency {
class TaskError final : public std::runtime_error {
public:
  TaskError(std::string pass, std::size_t bin, std::string detail)
      : std::runtime_error(message(pass, bin, detail)), pass_(std::move(pass)), bin_(bin) {}
  const std::string &pass() const noexcept { return pass_; }
  std::size_t bin() const noexcept { return bin_; }
private:
  static std::string message(const std::string &pass, std::size_t bin,
                             const std::string &detail) {
    std::ostringstream out;
    out << "frequency task failed: pass=" << pass << " bin=" << bin << ": " << detail;
    return out.str();
  }
  std::string pass_;
  std::size_t bin_;
};

// Indexed ownership, static scheduling, join-before-error, lowest-bin rethrow.
template<class ContextFactory, class Function>
void runIndexedWithContext(std::size_t count, int workers, std::size_t firstBin,
                const std::string &pass, ContextFactory &&contextFactory,
                Function &&function) {
  if (workers < 1 || workers > 8)
    throw std::invalid_argument("frequency worker count must be in [1,8]");
  if (count > static_cast<std::size_t>(std::numeric_limits<long long>::max()))
    throw std::overflow_error("frequency task count exceeds OpenMP loop range");
  if (firstBin > std::numeric_limits<std::size_t>::max() - count)
    throw std::overflow_error("frequency bin range overflows size_t");
  if (count == 0) return;
  static_assert(noexcept(contextFactory(0)), "frequency worker context factory must not throw");
  std::vector<std::exception_ptr> failures(count);
  int actualWorkers = 0;
  const auto task = std::forward<Function>(function);
#pragma omp parallel num_threads(workers) shared(actualWorkers, failures, task, contextFactory)
  {
    auto context = contextFactory(omp_get_thread_num());
    static_assert(std::is_nothrow_destructible_v<decltype(context)>, "frequency worker context destructor must not throw");
#pragma omp single
    actualWorkers = omp_get_num_threads();
#pragma omp for schedule(static)
    for (long long index = 0; index < static_cast<long long>(count); ++index) {
      try {
        task(static_cast<std::size_t>(index), context);
      } catch (...) {
        failures[static_cast<std::size_t>(index)] = std::current_exception();
      }
    }
  }
  if (actualWorkers != workers)
    throw std::runtime_error("OpenMP frequency team size differs from requested workers");
  for (std::size_t index = 0; index < failures.size(); ++index) {
    if (!failures[index]) continue;
    const std::size_t bin = firstBin + index;
    try {
      std::rethrow_exception(failures[index]);
    } catch (const std::exception &error) {
      throw TaskError(pass, bin, error.what());
    } catch (...) {
      throw TaskError(pass, bin, "non-standard exception");
    }
  }
}

template<class Function>
void runIndexed(std::size_t count, int workers, std::size_t firstBin,
                const std::string &pass, Function &&function) {
  runIndexedWithContext(count, workers, firstBin, pass,
      [](int) noexcept { return 0; },
      [&](std::size_t index, int &) { function(index); });
}
} // namespace dspec::frequency

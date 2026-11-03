#pragma once

// Progress reported from inside a threaded loop: one call at a time, and an
// exception from the callback, which may be how a caller stops the
// calculation, caught inside the critical section, which it may not leave,
// and thrown again outside it.

#include <YSpec/Calculation.hpp>
#include <cstddef>
#include <exception>

namespace YSpec {

inline void Report(const Progress& progress, std::size_t done,
                   std::size_t total) {
  if (!progress) return;
  auto failure = std::exception_ptr{};
#pragma omp critical(YSpecProgress)
  try {
    progress(done, total);
  } catch (...) {
    failure = std::current_exception();
  }
  if (failure) std::rethrow_exception(failure);
}

}  // namespace YSpec

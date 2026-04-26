#ifndef TIMER_CLASS_GUARD_H
#define TIMER_CLASS_GUARD_H
#include <chrono>
#include <iostream>
#include <string>

class Timer {
private:
  std::chrono::time_point<std::chrono::high_resolution_clock> _start, _stop;
  std::chrono::microseconds _duration;

public:
  Timer()
      : _start{std::chrono::high_resolution_clock::now()}, _stop{_start},
        _duration{0} {}

  void start() { _start = std::chrono::high_resolution_clock::now(); }

  void stop() {
    _stop = std::chrono::high_resolution_clock::now();
    _duration =
        std::chrono::duration_cast<std::chrono::microseconds>(_stop - _start);
  }

  /// Returns elapsed time in seconds.
  double duration_seconds() const {
    return static_cast<double>(_duration.count()) / 1'000'000.0;
  }

  void stop(const std::string &message) {
    stop();
    std::cout << message << ": " << duration_seconds() << " seconds\n";
  }
};

#endif
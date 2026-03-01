#ifndef TIMER_CLASS_GUARD_H
#define TIMER_CLASS_GUARD_H
#include <chrono>
#include <iostream>

using namespace std::chrono;
class Timer {
 private:
   std::chrono::time_point<std::chrono::_V2::system_clock,
                           std::chrono::_V2::system_clock::duration>
       _start, _stop;
   std::chrono::__enable_if_is_duration<std::chrono::microseconds> _duration;

 public:
   Timer()
       : _start{high_resolution_clock::now()},
         _stop{high_resolution_clock::now()},
         _duration{duration_cast<microseconds>(_stop - _start)} {};

   void start() { _start = high_resolution_clock::now(); };

   void stop() { _stop = high_resolution_clock::now(); };

   double duration() { return static_cast<double>(_duration.count()); };

   void stop(const std::string &strmessage) {
      _stop = high_resolution_clock::now();
      _duration = duration_cast<microseconds>(_stop - _start);
      std::cout << strmessage << ": " << _duration.count() / 1000000.0
                << " seconds\n";
   };
};

#endif
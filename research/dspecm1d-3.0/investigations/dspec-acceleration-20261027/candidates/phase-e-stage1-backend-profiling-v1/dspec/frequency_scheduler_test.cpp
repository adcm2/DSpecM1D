#include "frequency_scheduler.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
}
int main() {
  using namespace std::chrono_literals;
  std::vector<int> slots(8, -1);
  std::vector<std::size_t> completionOrder;
  std::mutex completionMutex;
  dspec::frequency::runIndexed(slots.size(), 4, 20, "ordering", [&](std::size_t i) {
    if (i == 0) std::this_thread::sleep_for(60ms);
    slots[i] = static_cast<int>(i * 11);
    std::lock_guard<std::mutex> lock(completionMutex);
    completionOrder.push_back(i);
  });
  require(completionOrder.size() == slots.size(), "not every frequency task completed");
  require(completionOrder.front() != 0, "test did not induce out-of-order completion");
  for (std::size_t i = 0; i < slots.size(); ++i)
    require(slots[i] == static_cast<int>(i * 11), "indexed slot assignment changed");
  std::vector<int> collected;
  for (int value : slots) collected.push_back(value);
  require(collected.front() == 0 && collected.back() == 77,
          "parent indexed collection was not deterministic");

  std::atomic<int> completed{0};
  bool returned = false;
  try {
    dspec::frequency::runIndexed(8, 4, 100, "multiple-errors", [&](std::size_t i) {
      if (i == 1) std::this_thread::sleep_for(25ms);
      ++completed;
      if (i == 1 || i == 5) throw std::runtime_error("injected bin failure");
    });
    returned = true;
  } catch (const dspec::frequency::TaskError &error) {
    require(error.pass() == "multiple-errors", "wrong pass in deterministic failure");
    require(error.bin() == 101, "did not select lowest failing frequency bin");
    require(std::string(error.what()).find("bin=101") != std::string::npos,
            "failure message omitted bin context");
  }
  require(!returned, "a failed pass returned partial success");
  require(completed == 8, "scheduler stopped dispatched work before join");

  struct PassContext {
    std::atomic<int> &completed;
    std::atomic<int> &joined;
    explicit PassContext(std::atomic<int> &done,std::atomic<int> &closed) noexcept:completed(done),joined(closed){}
    ~PassContext() noexcept {++joined;}
    PassContext(const PassContext&)=delete;
  };
  std::atomic<int> contexts{0},joined{0},contextTasks{0};
  dspec::frequency::runIndexedWithContext(20,4,0,"context-lifetime",
      [&](int) noexcept {++contexts;return PassContext(contextTasks,joined);},
      [&](std::size_t,PassContext &context){++context.completed;});
  require(contexts==4&&joined==4,"worker context did not span exactly one outer pass lifetime per worker");
  require(contextTasks==20,"worker context did not cover every task in its pass");

  completed = 0;
  returned = false;
  auto injectedCacheCallback = [](std::size_t i) {
    if (i == 3) throw std::overflow_error("injected response-cache capacity failure");
  };
  try {
    dspec::frequency::runIndexed(8, 4, 200, "cache-callback", [&](std::size_t i) {
      injectedCacheCallback(i);
      ++completed;
    });
    returned = true;
  } catch (const dspec::frequency::TaskError &error) {
    require(error.pass() == "cache-callback" && error.bin() == 203,
            "cache callback error lost frequency context");
    require(std::string(error.what()).find("capacity") != std::string::npos,
            "cache callback detail was not propagated");
  }
  require(!returned, "cache failure returned partial success");
  require(completed == 7, "cache failure did not join the remaining workers");
  return 0;
}

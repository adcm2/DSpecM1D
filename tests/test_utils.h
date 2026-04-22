#ifndef DSPECM1D_TEST_UTILS_H
#define DSPECM1D_TEST_UTILS_H

#include <cstdlib>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>
#include <unistd.h>
#include "config.h"

namespace DSpecMTest {

inline std::filesystem::path
repoRoot() {
  return std::filesystem::path(DSPECM1D_TEST_SOURCE_DIR);
}

inline std::filesystem::path
modelPath() {
  return repoRoot() / "data" / "models" / "prem.200.noatten.txt";
}

inline std::string
makeParameterText(const std::string &earthModelPath, int outputType = 0,
                  double relativeError = 1e-5) {
  return "\"./output/test.out\"\n"
         "\"" + earthModelPath + "\"\n" +
         "4\n"
         "0\n"
         "2\n" +
         std::to_string(outputType) + "\n"
         "0\n" + std::to_string(relativeError) + "\n"
         "0\n"
         "6\n"
         "0.1\n"
         "1.0\n"
         "5.0\n"
         "1.0\n"
         "0.1\n"
         "0.2\n"
         "0.8\n"
         "1.0\n"
         "33.0\n"
         "10.0\n"
         "20.0\n"
         "1.0\n"
         "2.0\n"
         "3.0\n"
         "4.0\n"
         "5.0\n"
         "6.0\n"
         "0.0\n"
         "2\n"
         "45.0 90.0\n"
         "-10.0 120.0\n";
}

inline std::filesystem::path
writeFile(const std::filesystem::path &path, const std::string &contents) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  out << contents;
  return path;
}

class TempDir {
public:
  TempDir() {
    namespace fs = std::filesystem;
    std::mt19937_64 rng(
        static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count()) ^
        static_cast<std::uint64_t>(::getpid()) ^
        static_cast<std::uint64_t>(
            std::hash<std::thread::id>{}(std::this_thread::get_id())));
    std::uniform_int_distribution<std::uint64_t> dist;

    const fs::path base = fs::temp_directory_path();
    for (int attempt = 0; attempt < 32; ++attempt) {
      std::ostringstream name;
      name << "dspecm1d-" << ::getpid() << "-" << std::hex << dist(rng);
      fs::path candidate = base / name.str();
      std::error_code ec;
      if (fs::create_directories(candidate, ec)) {
        m_path = std::move(candidate);
        return;
      }
    }

    throw std::runtime_error("Failed to create unique temporary directory");
  }

  ~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(m_path, ec);
  }

  const std::filesystem::path &path() const { return m_path; }

private:
  std::filesystem::path m_path;
};

}   // namespace DSpecMTest

#endif

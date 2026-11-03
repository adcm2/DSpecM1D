#include <YSpec/Instrument.hpp>
#include <algorithm>
#include <cctype>
#include <format>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

namespace YSpec {

namespace {

std::string Upper(std::string text) {
  std::ranges::transform(text, text.begin(), [](unsigned char c) {
    return static_cast<char>(std::toupper(c));
  });
  return text;
}

[[noreturn]] void Fail(const std::filesystem::path& path, int line,
                       const std::string& message) {
  throw std::runtime_error(
      std::format("{}:{}: {}", path.string(), line, message));
}

std::string Trim(const std::string& text) {
  const auto first = text.find_first_not_of(" \t\r");
  if (first == std::string::npos) return {};
  const auto last = text.find_last_not_of(" \t\r");
  return text.substr(first, last - first + 1);
}

}  // namespace

Instrument ReadSacPoleZero(const std::filesystem::path& path) {
  auto file = std::ifstream(path);
  if (!file) throw std::runtime_error("cannot read " + path.string());
  auto number = 0;

  auto instrument = Instrument{};
  auto unit = std::optional<Quantity>{};
  // The list being read, the number it may hold, and how many it holds.
  std::vector<Complex>* list = nullptr;
  auto count = std::size_t{0};
  auto listed = std::size_t{0};
  auto seen = std::vector<std::string>{};

  for (auto line = std::string{}; std::getline(file, line);) {
    ++number;
    const auto text = Trim(line);
    if (text.empty()) continue;
    if (text.front() == '*') {
      // A comment, which may give the input unit.
      const auto upper = Upper(text);
      const auto colon = upper.find(':');
      if (upper.find("INPUT UNIT") != std::string::npos &&
          colon != std::string::npos) {
        const auto value = Trim(upper.substr(colon + 1));
        if (value == "M") {
          unit = Quantity::Displacement;
        } else if (value == "M/S") {
          unit = Quantity::Velocity;
        } else if (value == "M/S**2" || value == "M/S/S" || value == "M/S^2") {
          unit = Quantity::Acceleration;
        } else {
          Fail(path, number,
               std::format("the input unit is {}, which is not M, M/S or "
                           "M/S**2",
                           value));
        }
      }
      continue;
    }

    auto stream = std::istringstream(text);
    auto first = std::string{};
    stream >> first;
    const auto keyword = Upper(first);
    if (keyword == "ZEROS" || keyword == "POLES" || keyword == "CONSTANT") {
      if (std::ranges::find(seen, keyword) != seen.end()) {
        Fail(path, number,
             "a second " + keyword + ": the file holds more than one response");
      }
      seen.push_back(keyword);
      list = nullptr;
      if (keyword == "CONSTANT") {
        if (!(stream >> instrument.constant))
          Fail(path, number, "CONSTANT needs a value");
        continue;
      }
      auto size = 0L;
      if (!(stream >> size) || size < 0) {
        Fail(path, number, keyword + " needs a number of entries");
      }
      list = keyword == "ZEROS" ? &instrument.zeros : &instrument.poles;
      count = static_cast<std::size_t>(size);
      listed = 0;
      list->assign(count, Complex{});
      continue;
    }

    // The real and imaginary parts of a zero or a pole.
    auto values = std::istringstream(text);
    auto real = 0.0;
    auto imaginary = 0.0;
    if (!(values >> real >> imaginary)) {
      Fail(path, number, "expected a keyword or two numbers");
    }
    if (!list) Fail(path, number, "numbers outside a list of zeros or poles");
    if (listed == count)
      Fail(path, number, "more entries than the count given");
    (*list)[listed++] = {real, imaginary};
  }
  if (seen.empty()) Fail(path, number, "holds no ZEROS, POLES or CONSTANT");
  if (unit) instrument.input = *unit;
  return instrument;
}

}  // namespace YSpec

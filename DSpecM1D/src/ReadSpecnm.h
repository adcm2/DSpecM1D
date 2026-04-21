#ifndef READ_SPECNM_H
#define READ_SPECNM_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace SPECNMREADER {

class DataColumns {
  std::vector<double> column1;
  std::vector<double> column2;
  std::vector<double> column3;
  std::vector<double> column4;

public:
  DataColumns() = default;
  explicit DataColumns(const std::string &filename);

  const std::vector<double> &getColumn1() const;
  const std::vector<double> &getColumn2() const;
  const std::vector<double> &getColumn3() const;
  const std::vector<double> &getColumn4() const;

  void clear();
};

inline DataColumns::DataColumns(const std::string &filename) {
  std::ifstream infile(filename);
  std::string line;
  int line_number = 0;

  if (!infile.is_open()) {
    throw std::runtime_error("Error: Could not open file " + filename);
  }

  while (std::getline(infile, line)) {
    ++line_number;

    if (line.find_first_not_of(" \t\n\v\f\r") == std::string::npos) {
      continue;
    }

    std::stringstream ss(line);
    std::string token;
    std::vector<double> values;
    values.reserve(4);

    bool parseFailed = false;
    while (std::getline(ss, token, ';')) {
      try {
        values.push_back(std::stod(token));
      } catch (const std::exception &) {
        parseFailed = true;
        break;
      }
    }

    if (!parseFailed && values.size() == 4) {
      column1.push_back(values[0]);
      column2.push_back(values[1]);
      column3.push_back(values[2]);
      column4.push_back(values[3]);
    } else {
      std::cerr << "Warning: Could not parse semicolon-delimited line "
                << line_number << ": '" << line << "'. Skipping line."
                << std::endl;
    }
  }

  infile.close();

  if (column1.empty()) {
    std::cerr << "Warning: No valid data lines read from file " << filename
              << std::endl;
  }
}

inline const std::vector<double> &
DataColumns::getColumn1() const {
  return column1;
}

inline const std::vector<double> &
DataColumns::getColumn2() const {
  return column2;
}

inline const std::vector<double> &
DataColumns::getColumn3() const {
  return column3;
}

inline const std::vector<double> &
DataColumns::getColumn4() const {
  return column4;
}

inline void
DataColumns::clear() {
  column1.clear();
  column2.clear();
  column3.clear();
  column4.clear();
}

}   // namespace SPECNMREADER

namespace SpecnmReader = SPECNMREADER;

#endif

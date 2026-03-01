#ifndef READ_YSPEC_H
#define READ_YSPEC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>     // Required for stringstream
#include <stdexcept>   // Required for runtime_error

namespace YSPECREADER {
// The DataColumns class with constructors
class DataColumns {
  // Member variables are now private by default in a class
  std::vector<double> column1;
  std::vector<double> column2;
  std::vector<double> column3;
  std::vector<double> column4;

public:
  // --- Constructors ---
  // 1. Default constructor: Creates an empty DataColumns object
  DataColumns() = default;

  // 2. Constructor that reads data from a file (declaration only)
  explicit DataColumns(const std::string &filename);

  // --- Getter methods (declarations only) ---
  const std::vector<double> &getColumn1() const;
  const std::vector<double> &getColumn2() const;
  const std::vector<double> &getColumn3() const;
  const std::vector<double> &getColumn4() const;

  // --- Other methods (declaration only) ---
  void clear();
};

// --- Member Function Definitions ---

// 2. Definition of the file-reading constructor
inline DataColumns::DataColumns(const std::string &filename) {
  std::ifstream infile(filename);
  std::string line;
  int line_number = 0;

  // Check if the file was opened successfully
  if (!infile.is_open()) {
    // Constructors should throw exceptions on failure, not just print errors
    throw std::runtime_error("Error: Could not open file " + filename);
  }

  // Read the file line by line
  while (std::getline(infile, line)) {
    line_number++;
    std::stringstream ss(line);   // Create a stringstream from the line
    double val1, val2, val3, val4;

    // Attempt to extract four double values from the stringstream
    if (ss >> val1 >> val2 >> val3 >> val4) {
      // Check if there's any extra non-whitespace data left on the line
      std::string remaining;
      if (ss >> remaining) {
        std::cerr << "Warning: Extra data found on line " << line_number
                  << ": '" << line << "'. Skipping line." << std::endl;
      } else {
        // Successfully read four values, add them directly to member vectors
        column1.push_back(val1);
        column2.push_back(val2);
        column3.push_back(val3);
        column4.push_back(val4);
      }
    } else {
      // Handle empty lines or lines that don't contain exactly four numbers
      // Only issue a warning if the line wasn't just whitespace
      if (!line.empty() &&
          line.find_first_not_of(" \t\n\v\f\r") != std::string::npos) {
        std::cerr << "Warning: Could not parse line " << line_number << ": '"
                  << line << "'. Skipping line." << std::endl;
      }
    }
  }

  infile.close();   // Close the file

  // Check if any data was read
  if (column1.empty()) {
    std::cerr << "Warning: No valid data lines read from file " << filename
              << std::endl;
  }
}

// --- Getter method definitions ---

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

// --- Other method definitions ---
inline void
DataColumns::clear() {
  column1.clear();
  column2.clear();
  column3.clear();
  column4.clear();
}

}   // namespace YSPECREADER
#endif
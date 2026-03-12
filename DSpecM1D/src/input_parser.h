#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>   // For std::pair
#include <vector>

// Helper function to read the next non-comment, non-empty line from a file
inline std::string
get_next_value_line(std::ifstream &file) {
  std::string line;
  while (std::getline(file, line)) {
    // Trim leading and trailing whitespace
    size_t first = line.find_first_not_of(" \t\n\r");
    if (std::string::npos == first) {
      continue;
    }
    size_t last = line.find_last_not_of(" \t\n\r");
    line = line.substr(first, (last - first + 1));

    if (!line.empty() && line[0] != '#') {
      // If the line is quoted, remove the quotes
      if (line.length() >= 2 && line.front() == '"' && line.back() == '"') {
        return line.substr(1, line.length() - 2);
      }
      return line;
    }
  }
  return "";   // Return empty if nothing is found
}

// Class to store all parameters from the input file
class InputParameters {
private:
  std::string m_output_prefix;
  std::string m_earth_model;
  int m_attenuation;
  int m_gravitation;
  int m_output_type;
  int m_corrections;
  int m_lmin;
  int m_lmax;
  int m_type;
  double m_f1;
  double m_f2;
  double m_t_out;
  double m_time_step_sec;
  double m_f11, m_f12, m_f21, m_f22;   // Filter params
  double m_source_depth_km;
  double m_source_lat_deg;
  double m_source_lon_deg;
  double m_m_rr, m_m_rt, m_m_rp, m_m_tt, m_m_tp, m_m_pp;   // Moment tensor
  double m_receiver_depth;
  double relerr;
  int m_num_receivers;
  std::vector<std::pair<double, double>> m_receivers;   // lat, lon

public:
  // Constructor that reads the file
  explicit InputParameters(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open input file: " + filename);
    }

    m_output_prefix = get_next_value_line(file);
    m_earth_model = get_next_value_line(file);
    std::stringstream(get_next_value_line(file)) >> m_type;
    std::stringstream(get_next_value_line(file)) >> m_attenuation;
    std::stringstream(get_next_value_line(file)) >> m_gravitation;
    std::stringstream(get_next_value_line(file)) >> m_output_type;
    std::stringstream(get_next_value_line(file)) >> m_corrections;
    std::stringstream(get_next_value_line(file)) >> relerr;
    std::stringstream(get_next_value_line(file)) >> m_lmin;
    std::stringstream(get_next_value_line(file)) >> m_lmax;
    std::stringstream(get_next_value_line(file)) >> m_f1;
    std::stringstream(get_next_value_line(file)) >> m_f2;
    std::stringstream(get_next_value_line(file)) >> m_t_out;
    std::stringstream(get_next_value_line(file)) >> m_time_step_sec;
    std::stringstream(get_next_value_line(file)) >> m_f11;
    std::stringstream(get_next_value_line(file)) >> m_f12;
    std::stringstream(get_next_value_line(file)) >> m_f21;
    std::stringstream(get_next_value_line(file)) >> m_f22;
    std::stringstream(get_next_value_line(file)) >> m_source_depth_km;
    std::stringstream(get_next_value_line(file)) >> m_source_lat_deg;
    std::stringstream(get_next_value_line(file)) >> m_source_lon_deg;
    std::stringstream(get_next_value_line(file)) >> m_m_rr;
    std::stringstream(get_next_value_line(file)) >> m_m_rt;
    std::stringstream(get_next_value_line(file)) >> m_m_rp;
    std::stringstream(get_next_value_line(file)) >> m_m_tt;
    std::stringstream(get_next_value_line(file)) >> m_m_tp;
    std::stringstream(get_next_value_line(file)) >> m_m_pp;
    std::stringstream(get_next_value_line(file)) >> m_receiver_depth;
    std::stringstream(get_next_value_line(file)) >> m_num_receivers;

    // Read receiver locations
    for (int i = 0; i < m_num_receivers; ++i) {
      std::string line = get_next_value_line(file);
      if (line.empty()) {
        throw std::runtime_error(
            "Input file ended unexpectedly while reading receivers.");
      }
      std::stringstream ss(line);
      double lat, lon;
      ss >> lat >> lon;
      m_receivers.push_back({lat, lon});
    }
  }

  // --- Public Getter Functions ---
  const std::string &output_prefix() const { return m_output_prefix; }
  const std::string &earth_model() const { return m_earth_model; }
  int type() const { return m_type; }
  int attenuation() const { return m_attenuation; }
  int gravitation() const { return m_gravitation; }
  int output_type() const { return m_output_type; }
  int corrections() const { return m_corrections; }
  double relative_error() const { return relerr; }
  int lmin() const { return m_lmin; }
  int lmax() const { return m_lmax; }
  double f1() const { return m_f1; }
  double f2() const { return m_f2; }
  double t_out() const { return m_t_out; }
  double time_step_sec() const { return m_time_step_sec; }
  double f11() const { return m_f11; }
  double f12() const { return m_f12; }
  double f21() const { return m_f21; }
  double f22() const { return m_f22; }
  double source_depth_km() const { return m_source_depth_km; }
  double source_lat_deg() const { return m_source_lat_deg; }
  double source_lon_deg() const { return m_source_lon_deg; }
  double m_rr() const { return m_m_rr; }
  double m_rt() const { return m_m_rt; }
  double m_rp() const { return m_m_rp; }
  double m_tt() const { return m_m_tt; }
  double m_tp() const { return m_m_tp; }
  double m_pp() const { return m_m_pp; }
  double receiver_depth() const { return m_receiver_depth; }
  int num_receivers() const { return m_num_receivers; }
  const std::vector<std::pair<double, double>> &receivers() const {
    return m_receivers;
  }
};

#endif   // INPUT_PARSER_H
#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>   // For std::pair
#include <vector>
#include <stdexcept>
#include <type_traits>

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

// Reads next value line and parses exactly one scalar of type T.
// Throws with a field name if missing/invalid/trailing tokens exist.
template <typename T>
inline T
read_required_scalar(std::ifstream &file, const char *field_name) {
  const std::string line = get_next_value_line(file);
  if (line.empty()) {
    throw std::runtime_error(std::string("Missing value for field: ") +
                             field_name);
  }

  std::stringstream ss(line);
  T value{};
  if (!(ss >> value)) {
    throw std::runtime_error(std::string("Invalid value for field: ") +
                             field_name + " (line: \"" + line + "\")");
  }

  // Reject trailing junk: "10 abc"
  ss >> std::ws;
  if (!ss.eof()) {
    throw std::runtime_error(std::string("Trailing tokens for field: ") +
                             field_name + " (line: \"" + line + "\")");
  }

  return value;
}

// Helper function to read a latitude and longitude from a file
inline std::pair<double, double>
read_required_lat_lon(std::ifstream &file, const char *field_name) {
  const std::string line = get_next_value_line(file);
  if (line.empty()) {
    throw std::runtime_error(std::string("Missing value for field: ") +
                             field_name);
  }

  std::stringstream ss(line);
  double lat = 0.0, lon = 0.0;
  if (!(ss >> lat >> lon)) {
    throw std::runtime_error(std::string("Invalid lat/lon for field: ") +
                             field_name + " (line: \"" + line + "\")");
  }

  ss >> std::ws;
  if (!ss.eof()) {
    throw std::runtime_error(std::string("Trailing tokens for field: ") +
                             field_name + " (line: \"" + line + "\")");
  }

  return {lat, lon};
}

inline void
require_condition(bool ok, const char *message) {
  if (!ok) {
    throw std::runtime_error(message);
  }
}

template <typename T>
inline void
require_in_range(T value, T lo, T hi, const char *field_name) {
  if (value < lo || value > hi) {
    throw std::runtime_error(std::string(field_name) + " out of range [" +
                             std::to_string(lo) + ", " + std::to_string(hi) +
                             "], got " + std::to_string(value));
  }
}

template <typename T>
inline void
require_positive(T value, const char *field_name) {
  if (!(value > static_cast<T>(0))) {
    throw std::runtime_error(std::string(field_name) + " must be > 0, got " +
                             std::to_string(value));
  }
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
  double m_relative_error;
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
    if (m_output_prefix.empty())
      throw std::runtime_error("Missing value for field: output_prefix");

    m_earth_model = get_next_value_line(file);
    if (m_earth_model.empty())
      throw std::runtime_error("Missing value for field: earth_model");

    m_type = read_required_scalar<int>(file, "type");
    m_attenuation = read_required_scalar<int>(file, "attenuation");
    m_gravitation = read_required_scalar<int>(file, "gravitation");
    m_output_type = read_required_scalar<int>(file, "output_type");
    m_corrections = read_required_scalar<int>(file, "corrections");
    m_relative_error = read_required_scalar<double>(file, "relative_error");
    m_lmin = read_required_scalar<int>(file, "lmin");
    m_lmax = read_required_scalar<int>(file, "lmax");
    m_f1 = read_required_scalar<double>(file, "f1");
    m_f2 = read_required_scalar<double>(file, "f2");
    m_t_out = read_required_scalar<double>(file, "t_out");
    m_time_step_sec = read_required_scalar<double>(file, "time_step_sec");

    m_f11 = read_required_scalar<double>(file, "f11");
    m_f12 = read_required_scalar<double>(file, "f12");
    m_f21 = read_required_scalar<double>(file, "f21");
    m_f22 = read_required_scalar<double>(file, "f22");

    m_source_depth_km = read_required_scalar<double>(file, "source_depth_km");
    m_source_lat_deg = read_required_scalar<double>(file, "source_lat_deg");
    m_source_lon_deg = read_required_scalar<double>(file, "source_lon_deg");

    m_m_rr = read_required_scalar<double>(file, "m_rr");
    m_m_rt = read_required_scalar<double>(file, "m_rt");
    m_m_rp = read_required_scalar<double>(file, "m_rp");
    m_m_tt = read_required_scalar<double>(file, "m_tt");
    m_m_tp = read_required_scalar<double>(file, "m_tp");
    m_m_pp = read_required_scalar<double>(file, "m_pp");

    m_receiver_depth = read_required_scalar<double>(file, "receiver_depth");
    m_num_receivers = read_required_scalar<int>(file, "num_receivers");

    // Read receiver locations
    for (int i = 0; i < m_num_receivers; ++i) {
      m_receivers.push_back(read_required_lat_lon(file, "receiver_lat_lon"));
    }

    require_condition(m_num_receivers > 0, "num_receivers must be > 0");
    require_condition(m_lmin >= 0 && m_lmax >= m_lmin, "Invalid lmin/lmax");
    require_positive(m_time_step_sec, "time_step_sec");
    require_positive(m_t_out, "t_out");
    require_positive(m_relative_error, "relative_error");
    require_condition(m_receiver_depth >= 0.0, "receiver_depth must be >= 0");
    require_condition(m_source_depth_km >= 0.0, "source_depth_km must be >= 0");

    require_in_range(m_source_lat_deg, -90.0, 90.0, "source_lat_deg");
    require_in_range(m_source_lon_deg, -180.0, 360.0, "source_lon_deg");

    require_condition(m_f11 <= m_f12 && m_f1 < m_f2 && m_f21 <= m_f22,
                      "Invalid frequency/taper ordering");

    for (int i = 0; i < m_num_receivers; ++i) {
      const auto [lat, lon] = m_receivers[i];
      require_in_range(lat, -90.0, 90.0, "receiver latitude");
      require_in_range(lon, -180.0, 360.0, "receiver longitude");
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
  double relative_error() const { return m_relative_error; }
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
#ifndef READ_STATION_H
#define READ_STATION_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @class SiteChanEntry
 * @brief Represents all data from a single line of a sitechan file.
 *
 * This class encapsulates the properties of a station channel entry.
 * Data is private and accessed through public getter methods.
 */
class SiteChanEntry {
private:
  std::string m_station;
  std::string m_channel;
  int m_jdate;
  int m_ondate;
  int m_offdate;
  char m_chanflag;
  double m_elevation_depth;
  double m_azimuth;
  double m_zenith;
  char m_description;
  std::string m_timestamp;

public:
  // Constructor to initialize all members
  SiteChanEntry(std::string station, std::string channel, int jdate, int ondate,
                int offdate, char chanflag, double elevation_depth,
                double azimuth, double zenith, char description,
                std::string timestamp)
      : m_station(std::move(station)), m_channel(std::move(channel)),
        m_jdate(jdate), m_ondate(ondate), m_offdate(offdate),
        m_chanflag(chanflag), m_elevation_depth(elevation_depth),
        m_azimuth(azimuth), m_zenith(zenith), m_description(description),
        m_timestamp(std::move(timestamp)) {}

  // Public getter functions to access private members
  const std::string &station() const { return m_station; }
  const std::string &channel() const { return m_channel; }
  int jdate() const { return m_jdate; }
  int ondate() const { return m_ondate; }
  int offdate() const { return m_offdate; }
  char chanflag() const { return m_chanflag; }
  double elevation_depth() const { return m_elevation_depth; }
  double azimuth() const { return m_azimuth; }
  double zenith() const { return m_zenith; }
  char description() const { return m_description; }
  const std::string &timestamp() const { return m_timestamp; }
};

/**
 * @brief Reads all data from a sitechan file.
 * @param filename The path to the input file.
 * @return A vector of SiteChanEntry objects, each containing the data from one
 * line.
 */
std::vector<SiteChanEntry>
read_full_sitechan_file(const std::string &filename) {
  std::vector<SiteChanEntry> entries;
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return entries;   // Return empty vector on error
  }

  std::string line;
  while (std::getline(file, line)) {
    std::stringstream ss(line);

    // Temporary variables to hold parsed data
    std::string station, channel, timestamp;
    int jdate, ondate, offdate;
    char chanflag, description;
    double elevation_depth, azimuth, zenith;

    // Read the fixed-format columns
    ss >> station >> channel >> jdate >> ondate >> offdate >> chanflag >>
        elevation_depth >> azimuth >> zenith >> description;

    if (ss.fail()) {
      std::cerr << "Warning: Skipping malformed line: " << line << std::endl;
      continue;
    }

    // Consume whitespace and read the rest of the line
    ss >> std::ws;
    std::getline(ss, timestamp);

    // Use the constructor to create and add the new entry
    entries.emplace_back(station, channel, jdate, ondate, offdate, chanflag,
                         elevation_depth, azimuth, zenith, description,
                         timestamp);
  }

  return entries;
}

#endif   // READ_STATION_H
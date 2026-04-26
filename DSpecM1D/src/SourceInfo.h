#ifndef SOURCEINFO_GUARD_H
#define SOURCEINFO_GUARD_H

#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "InputParser.h"

namespace SourceInfo {

// --- Class Declaration ---
class EarthquakeCMT {
public:
  // Constructors
  explicit EarthquakeCMT(const std::string &pathtocmt);
  EarthquakeCMT(double lat, double lon, double depth, double mrr, double mrt,
                double mrp, double mtt, double mtp, double mpp);
  explicit EarthquakeCMT(const InputParameters &params);   // New constructor

  // Function to print CMT information
  void PrintCMTInfo() const;

  // Getter functions
  const std::string &ID() const;
  int Year() const;
  int Day() const;
  int Hour() const;
  int Minute() const;
  double Second() const;
  double Latitude() const;
  double Longitude() const;
  double Depth() const;
  double Step() const;
  double HalfDuration() const;
  double M0() const;
  double Mrr() const;
  double Mtt() const;
  double Mpp() const;
  double Mrt() const;
  double Mrp() const;
  double Mtp() const;
  double MomentTensorNorm() const;
  double Strike1() const;
  double Dip1() const;
  double Slip1() const;
  double Strike2() const;
  double Dip2() const;
  double Slip2() const;

  // Canonical moment tensor functions
  std::complex<double> MC00() const;
  std::complex<double> MC0p() const;
  std::complex<double> MC0m() const;
  std::complex<double> MCpp() const;
  std::complex<double> MCmm() const;
  std::complex<double> MCmp() const;

private:
  std::string cmt_id;
  int cmt_yr = 0, cmt_day = 0, cmt_hr = 0, cmt_min = 0;
  double cmt_sec = 0.0, cmt_lat = 0.0, cmt_lon = 0.0, cmt_depth = 0.0,
         cmt_step = 0.0, cmt_halfd = 0.0, cmt_m0 = 0.0, cmt_mrr = 0.0,
         cmt_mtt = 0.0, cmt_mpp = 0.0, cmt_mrt = 0.0, cmt_mrp = 0.0,
         cmt_mtp = 0.0, crt_mn = 0.0;
  double cmt_s1 = 0.0, cmt_d1 = 0.0, cmt_sl1 = 0.0, cmt_s2 = 0.0, cmt_d2 = 0.0,
         cmt_sl2 = 0.0;
};

// --- Method Definitions ---

// Constructor (reading from file)
inline EarthquakeCMT::EarthquakeCMT(const std::string &pathtocmt) {
  std::ifstream modelfile(pathtocmt);
  if (modelfile.is_open()) {
    // extract information
    modelfile >> cmt_id >> cmt_yr >> cmt_day >> cmt_hr >> cmt_min >> cmt_sec >>
        cmt_lat >> cmt_lon >> cmt_depth >> cmt_step >> cmt_halfd >> cmt_m0 >>
        cmt_mrr >> cmt_mtt >> cmt_mpp >> cmt_mrt >> cmt_mrp >> cmt_mtp >>
        crt_mn >> cmt_s1 >> cmt_d1 >> cmt_sl1 >> cmt_s2 >> cmt_d2 >> cmt_sl2;
  } else {
    std::cerr << "Error: Could not open CMT file " << pathtocmt << std::endl;
  }
}

// Overloaded constructor for direct initialization
inline EarthquakeCMT::EarthquakeCMT(double lat, double lon, double depth,
                                    double mrr, double mrt, double mrp,
                                    double mtt, double mtp, double mpp)
    : cmt_lat(lat), cmt_lon(lon), cmt_depth(depth), cmt_mrr(mrr), cmt_mtt(mtt),
      cmt_mpp(mpp), cmt_mrt(mrt), cmt_mrp(mrp), cmt_mtp(mtp) {
  cmt_id = "DirectInput";
}

// New constructor from InputParameters object
inline EarthquakeCMT::EarthquakeCMT(const InputParameters &params)
    : EarthquakeCMT(params.source_lat_deg(), params.source_lon_deg(),
                    params.source_depth_km(), params.m_rr(), params.m_rt(),
                    params.m_rp(), params.m_tt(), params.m_tp(),
                    params.m_pp()) {
  // Body is empty due to constructor delegation
}

// Function to print CMT information
inline void
EarthquakeCMT::PrintCMTInfo() const {
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "--- CMT Source Information ---\n";
  std::cout << "ID: " << cmt_id << " (" << cmt_yr << "/" << cmt_day << " "
            << cmt_hr << ":" << cmt_min << ":" << cmt_sec << ")\n";
  std::cout << "Location (Lat/Lon/Depth): " << cmt_lat << " / " << cmt_lon
            << " / " << cmt_depth << " km\n";
  std::cout << "M0: " << cmt_m0 << "\n";
  std::cout << "Mrr: " << cmt_mrr << ", Mtt: " << cmt_mtt
            << ", Mpp: " << cmt_mpp << "\n";
  std::cout << "Mrt: " << cmt_mrt << ", Mrp: " << cmt_mrp
            << ", Mtp: " << cmt_mtp << "\n";
  std::cout << "------------------------------\n";
}

// --- Getter Definitions ---
inline const std::string &
EarthquakeCMT::ID() const {
  return cmt_id;
}
inline int
EarthquakeCMT::Year() const {
  return cmt_yr;
}
inline int
EarthquakeCMT::Day() const {
  return cmt_day;
}
inline int
EarthquakeCMT::Hour() const {
  return cmt_hr;
}
inline int
EarthquakeCMT::Minute() const {
  return cmt_min;
}
inline double
EarthquakeCMT::Second() const {
  return cmt_sec;
}
inline double
EarthquakeCMT::Latitude() const {
  return cmt_lat;
}
inline double
EarthquakeCMT::Longitude() const {
  return cmt_lon;
}
inline double
EarthquakeCMT::Depth() const {
  return cmt_depth;
}
inline double
EarthquakeCMT::Step() const {
  return cmt_step;
}
inline double
EarthquakeCMT::HalfDuration() const {
  return cmt_halfd;
}
inline double
EarthquakeCMT::M0() const {
  return cmt_m0;
}
inline double
EarthquakeCMT::Mrr() const {
  return cmt_mrr;
}
inline double
EarthquakeCMT::Mtt() const {
  return cmt_mtt;
}
inline double
EarthquakeCMT::Mpp() const {
  return cmt_mpp;
}
inline double
EarthquakeCMT::Mrt() const {
  return cmt_mrt;
}
inline double
EarthquakeCMT::Mrp() const {
  return cmt_mrp;
}
inline double
EarthquakeCMT::Mtp() const {
  return cmt_mtp;
}
inline double
EarthquakeCMT::MomentTensorNorm() const {
  return crt_mn;
}
inline double
EarthquakeCMT::Strike1() const {
  return cmt_s1;
}
inline double
EarthquakeCMT::Dip1() const {
  return cmt_d1;
}
inline double
EarthquakeCMT::Slip1() const {
  return cmt_sl1;
}
inline double
EarthquakeCMT::Strike2() const {
  return cmt_s2;
}
inline double
EarthquakeCMT::Dip2() const {
  return cmt_d2;
}
inline double
EarthquakeCMT::Slip2() const {
  return cmt_sl2;
}

// --- Canonical Moment Tensor Definitions ---
inline std::complex<double>
EarthquakeCMT::MC00() const {
  return cmt_mrr;
}
inline std::complex<double>
EarthquakeCMT::MC0p() const {
  return std::complex<double>(-cmt_mrt, cmt_mrp) / std::sqrt(2.0);
}
inline std::complex<double>
EarthquakeCMT::MC0m() const {
  return std::complex<double>(cmt_mrt, cmt_mrp) / std::sqrt(2.0);
}
inline std::complex<double>
EarthquakeCMT::MCpp() const {
  return std::complex<double>(0.5 * (cmt_mtt - cmt_mpp), -cmt_mtp);
}
inline std::complex<double>
EarthquakeCMT::MCmm() const {
  return std::complex<double>(0.5 * (cmt_mtt - cmt_mpp), cmt_mtp);
}
inline std::complex<double>
EarthquakeCMT::MCmp() const {
  return std::complex<double>(-0.5 * (cmt_mtt + cmt_mpp), 0.0);
}

}   // namespace SourceInfo

#endif   // SOURCEINFO_GUARD_H
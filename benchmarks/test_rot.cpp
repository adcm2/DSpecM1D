#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#include <iostream>
#include <new_coupling/Timer>
#include "input_parser.h"   // Use the new input parser
#include "SR_Info.h"
// #include <cstdlib>

template <typename FLOAT> class prem_norm {
public:
  FLOAT LengthNorm() const { return _length_norm; };
  FLOAT MassNorm() const { return _mass_norm; };
  FLOAT TimeNorm() const { return _time_norm; };

  prem_norm() {};

private:
  FLOAT _length_norm = 1000.0;
  FLOAT _mass_norm = 5515 * std::pow(_length_norm, 3.0);
  FLOAT _time_norm = 1 / std::sqrt(3.1415926535897932 * 6.67230e-11 * 5515);
};

int
main() {

  // Read all parameters from the input file_w
  InputParameters params("../YSpec/input_bench_mf.txt");
  SRInfo sr_info(params);
  std::cout << "Successful\n";
  return 0;
}
#endif
// #endif
#include <iostream>
#include <PlanetaryModel/All>
#include <new_coupling/Timer>
#include "sem_full.h"
#include "start_element.h"
// #include "sem_spheroidal_debug.h"
// #include "../SpectraSolver/SpectraSolver/FF"
// #include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/filter_base.h"
// #include
// "../SpectraSolver/SpectraSolver/src/ODE_Spectra/postprocessfunctions.h"
// #include "read_station.h"
// #include "input_parser.h"   // Use the new input parser
// #include "read_yspec.h"
// #include "read_mineos.h"
// #include "full_spec.h"

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
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // parameters of sem
  int lmin = 0;
  int lmax = 10;
  int NQ = 5;
  double maxstep = 0.05;
  double omega = 0.005;

  std::cout << "omega: \n";
  std::cin >> omega;

  std::cout << "lmin: \n";
  std::cin >> lmin;

  std::cout << "lmax: \n";
  std::cin >> lmax;

  prem_norm<double> norm_class;
  std::string earth_model_path = "../YSpec/examples/prem.200.no.noatten.txt";
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");
  Timer timer1;

  // sem class
  timer1.start();
  Full1D::sem sem(prem, maxstep, NQ, lmax);
  timer1.stop("Time to construct SEM");

  auto mesh = sem.mesh();
  //   for (int idxe = 0; idxe < mesh.NE(); ++idxe) {
  //     std::cout << "idxe: " << idxe << ", solid: " << mesh.IsSolid(idxe)
  //               << ", low: " << mesh.ELR(idxe) << ", up: " << mesh.EUR(idxe)
  //               << "\n";
  //   }

  timer1.start();
  for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
    auto tmpval = SpectralTools::StartElement_Sph(sem, prem, idxl, omega);
    auto tmp2 = SpectralTools::StartElement_S(sem, prem, idxl, omega);
    // std::cout << "\n";
  }
  timer1.stop("Time for StartElement function");

  return 0;
}
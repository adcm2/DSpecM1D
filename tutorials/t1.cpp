#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

/**
 * @file t1.cpp
 * @brief Tutorial 1: Basic seismogram synthesis using the Spectral Element
 * Method (SEM).
 *
 * @details
 * This tutorial walks through the complete workflow for computing synthetic
 * seismograms using the DSpecM1D library:
 *
 *  1. Read simulation parameters from an input file.
 *  2. Initialise the 1D Earth model (PREM).
 *  3. Set up the frequency axis and taper.
 *  4. Define the earthquake source (CMT).
 *  5. Compute the raw frequency-domain spectrum via the sparse SEM solver.
 *  6. Apply physical normalisation (displacement / velocity / acceleration).
 *  7. Transform back to the time domain and apply Hann-window filtering.
 *
 * @note
 * The input parameter file is expected at:
 *   `<build_dir>/data/params/ex1.txt`
 *
 * @author  Alex Myhill
 * @date    2026-03-06
 */

// --- Standard library ---
#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// --- Project headers ---
#include "config.h"
#include <DSpecM1D/All>
#include <DSpecM1D/Timer>
#include <PlanetaryModel/All>
#include <SpectraSolver/FF>

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// @brief Value of π to machine precision.
// constexpr double PI = 3.1415926535897932;

// ---------------------------------------------------------------------------
// Normalisation helper
// ---------------------------------------------------------------------------

/**
 * @brief Dimensional normalisation scales for PREM.
 *
 * @details
 * Stores the three fundamental scales used to non-dimensionalise the
 * governing equations:
 *
 * | Scale  | Value                                          |
 * |--------|------------------------------------------------|
 * | Length | 1 km = 1000 m                                  |
 * | Mass   | `ρ_ref × L³`  where `ρ_ref = 5515 kg m⁻³`     |
 * | Time   | `1 / √(π G ρ_ref)`                             |
 *
 * @tparam FLOAT Floating-point type (typically `double`).
 */
// template <typename FLOAT> class prem_norm {
// public:
//   prem_norm() = default;

//   /// @brief Returns the length normalisation scale (m).
//   FLOAT LengthNorm() const { return _length_norm; }

//   /// @brief Returns the mass normalisation scale (kg).
//   FLOAT MassNorm() const { return _mass_norm; }

//   /// @brief Returns the time normalisation scale (s).
//   FLOAT TimeNorm() const { return _time_norm; }

// private:
//   /// Length scale: 1 km in metres.
//   FLOAT _length_norm = 1000.0;

//   /// Mass scale derived from reference density and length scale.
//   FLOAT _mass_norm = 5515.0 * std::pow(_length_norm, 3.0);

//   /// Time scale derived from reference density and gravitational constant.
//   FLOAT _time_norm = 1.0 / std::sqrt(PI * 6.67230e-11 * 5515.0);
// };

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

/**
 * @brief Entry point for Tutorial 1.
 *
 * @details
 * Runs the full SEM seismogram synthesis pipeline and writes the filtered
 * time-domain seismograms to stdout. See the file-level documentation for
 * a step-by-step description.
 *
 * @return 0 on success.
 */
int
main() {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;

  // -------------------------------------------------------------------------
  // Step 1: Load input parameters
  // -------------------------------------------------------------------------

  /// Path to the parameter file, resolved relative to the CMake build tree.
  std::string param_path =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex1.txt";

  InputParameters params(param_path);
  SRInfo sr_info(params);

  std::string earth_model_path =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  // -------------------------------------------------------------------------
  // Step 2: SEM discretisation parameters
  // -------------------------------------------------------------------------

  /// Maximum angular degree for the spherical harmonic expansion.
  int lval = params.lmax();

  /// Number of Gauss–Lobatto–Legendre (GLL) quadrature points per element.
  int NQ = 5;

  /// Maximum element width (non-dimensional) used when building the radial
  /// mesh.
  double maxstep = 0.05;

  // -------------------------------------------------------------------------
  // Step 3: Frequency-domain solver parameters
  // -------------------------------------------------------------------------

  double dt = params.time_step_sec();    ///< Time step (s).
  double tout = params.t_out() / 60.0;   ///< Total output time (minutes).

  /// Frequency resolution of the raw spectrum (mHz).
  double df0 = 1.0;

  /// Width of the Hann taper applied in the initial freq→time step.
  double wtb = 0.05;

  double t1 = 0.0;    ///< Start time for the time window (s).
  double t2 = tout;   ///< End time for the time window (s).

  /// Exponent controlling the steepness of the cosine taper.
  int qex = 4;

  // -------------------------------------------------------------------------
  // Step 4: Initialise the Earth model and frequency class
  // -------------------------------------------------------------------------

  prem_norm<double> norm_class;

  /// 1D PREM model, non-dimensionalised using @p norm_class.
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  /**
   * @brief Frequency axis and Hann-taper helper.
   *
   * The constructor computes the non-dimensional angular frequency vector
   * and the corresponding spectral taper used throughout the pipeline.
   *
   * Arguments: f1, f2 (passband corners), f11, f12, f21, f22 (taper corners),
   * dt, tout, df0, wtb, t1, t2, qex, TimeNorm.
   */
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());

  auto vec_w = myff.w();   ///< Non-dimensional angular frequency vector.

  // -------------------------------------------------------------------------
  // Step 5: Define the earthquake source and compute the raw spectrum
  // -------------------------------------------------------------------------

  /// Centroid-moment-tensor (CMT) source read from the parameter file.
  auto cmt = SourceInfo::EarthquakeCMT(params);

  /**
   * @brief Sparse SEM frequency-domain solver.
   *
   * @c Spectra() iterates over all angular degrees @p l up to @p lval,
   * assembles the stiffness/inertia matrices, solves the linear system at
   * each frequency, and accumulates the receiver seismograms.
   *
   * @param myff        Frequency axis and taper.
   * @param prem        1D Earth model.
   * @param cmt         Earthquake CMT source.
   * @param params      Simulation parameters (receivers, tolerances, etc.).
   * @param NQ          GLL quadrature order.
   * @param sr_info     Source–receiver geometry.
   * @param rel_error   Relative tolerance for the iterative solver.
   *
   * @returns Matrix of complex spectra: rows = receivers × components,
   *          cols = frequencies.
   */
  SPARSESPEC::Sparse_F_Spec mytest;
  MATRIX vec_raw = mytest.Spectra(myff, prem, cmt, params, NQ, sr_info,
                                  params.relative_error());

  // -------------------------------------------------------------------------
  // Step 6: Physical normalisation
  // -------------------------------------------------------------------------

  /**
   * @brief Scale factor to convert from non-dimensional units to SI.
   *
   * | output_type | units          | norm_factor                  |
   * |-------------|----------------|------------------------------|
   * | 0           | displacement   | LengthNorm (m)               |
   * | 1           | velocity       | LengthNorm / TimeNorm (m/s)  |
   * | 2           | acceleration   | LengthNorm / TimeNorm² (m/s²)|
   */
  double accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  double norm_factor = 1.0;

  if (params.output_type() == 0) {
    norm_factor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    norm_factor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    norm_factor = accel_norm;
  }

  vec_raw *= norm_factor;

  // -------------------------------------------------------------------------
  // Step 7: Time-domain filtering pipeline
  // -------------------------------------------------------------------------

  /**
   * @details The filtering pipeline proceeds as follows:
   *
   *  (a) @c freq2time     — inverse FFT the raw spectrum to the time domain,
   * and multiply by exp(ep * t). (b) @c fulltime2freq — re-FFT with a Hann
   * taper (width @p wtb = 0.05) to remove aliasing. (c) @c
   * filtfreq2time — apply the bandpass filter defined by [f11, f12, f21, f22]
   * and inverse FFT. (d) @c fulltime2freq — final re-FFT with the
   * user-specified Hann taper (width @p hann_w = 0.5) for smooth spectra.
   */

  /// Hann-window half-width (wrt whole time series), ie 0.5 = taper over the
  /// first/last 50% of the time series.
  double hann_w = 0.5;

  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.05);
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);
  auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

  return 0;
}
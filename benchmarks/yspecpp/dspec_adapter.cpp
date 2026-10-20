#include <DSpecM1D/All>
#include <DSpecM1D/FrequencyTools>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  if (argc < 6 || argc > 8) {
    std::cerr << "usage: dspec_adapter INPUT_PARAMS OUTPUT_PREFIX NQ NSKIP MAXSTEP [paper2|paper_suite [dynamic|single]]\n";
    return 2;
  }
  try {
    const auto started = std::chrono::steady_clock::now();
    const bool paper2 = argc == 7 && std::string(argv[6]) == "paper2";
    const bool paperSuite = argc >= 7 && std::string(argv[6]) == "paper_suite";
    const std::string suiteMesh = paperSuite && argc == 8 ? argv[7] : "dynamic";
    if ((argc == 7 && !paper2 && !paperSuite) ||
        (argc == 8 && (!paperSuite ||
                       (suiteMesh != "dynamic" && suiteMesh != "single")))) {
      throw std::invalid_argument("unknown optional adapter mode");
    }
    InputParametersNew params(argv[1], std::stoi(argv[3]), std::stoi(argv[4]),
                              std::stod(argv[5]));
    SPARSESPEC::SparseFSpec solver;
    const bool forceCowling = params.inputParameters().gravitation() == 1;
    Eigen::MatrixXcd raw;
    if (paper2 || (paperSuite && suiteMesh == "dynamic")) {
      // The paper comparison uses the frequency-chunked dynamic multi-SEM API.
      raw = solver.spectra(params);
    } else if (paperSuite) {
      // Figure-specific fixed mesh cap, currently used by published Fig. 4.
      Full1D::SEM sem(params);
      SPARSESPEC::SpectraRunContext request(
          params.freqFull(), params.cmt(), params.inputParameters(),
          params.tref(), params.nskip());
      raw = solver.spectra(request, sem, forceCowling);
    } else {
      Full1D::SEM sem(params);
      SPARSESPEC::SpectraRunContext request(
          params.freqFull(), params.cmt(), params.inputParameters(),
          params.tref(), params.nskip());
      raw = solver.spectra(request, sem, forceCowling);
    }
    const auto &grid = params.freqFull();
    const auto timeNorm = params.earthModel().TimeNorm();
    const auto unit = params.normFactor();
    const Eigen::MatrixXd times = (paper2 || paperSuite)
        ? processfunctions::filtfreq2time(raw * unit, grid, true)
        : (processfunctions::freq2time(raw, grid) * unit).eval();

    std::ofstream spectra(std::string(argv[2]) + ".spectra.csv");
    std::ofstream series(std::string(argv[2]) + ".timeseries.csv");
    std::ofstream meta(std::string(argv[2]) + ".meta.json");
    if (!spectra || !series || !meta) {
      throw std::runtime_error("could not open adapter output");
    }
    const int receiverCount = params.inputParameters().num_receivers();
    std::ofstream allSpectra;
    std::ofstream allSeries;
    if (receiverCount > 1) {
      allSpectra.open(std::string(argv[2]) + ".all_receivers.spectra.csv");
      allSeries.open(std::string(argv[2]) + ".all_receivers.timeseries.csv");
      if (!allSpectra || !allSeries) {
        throw std::runtime_error("could not open all-receiver adapter output");
      }
      allSpectra << "frequency_bin,frequency_hz,omega_real_rad_s,omega_imag_rad_s,receiver_index,component,real,imag\n"
                 << std::setprecision(17);
      allSeries << "time_s" << std::setprecision(17);
      for (int receiver = 0; receiver < receiverCount; ++receiver) {
        for (const char *component : {"Z", "N", "E"})
          allSeries << ",r" << receiver << '_' << component << "_SI";
      }
      allSeries << '\n';
    }
    spectra << "frequency_bin,frequency_hz,omega_real_rad_s,omega_imag_rad_s,component,real,imag\n"
            << std::setprecision(17);
    const auto angular = grid.w();
    for (int index = grid.i1(); index < grid.i2(); ++index) {
      const auto fHz = grid.f(index) / timeNorm;
      const auto realOmega = angular[index] / timeNorm;
      const auto imagOmega = -grid.ep() / timeNorm;
      for (int component = 0; component < 3; ++component) {
        const auto value = raw(component, index) * unit * timeNorm;
        const auto bin = static_cast<long long>(index);
        spectra << bin << ',' << fHz << ',' << realOmega << ',' << imagOmega << ','
                << component << ',' << value.real() << ',' << value.imag()
                << '\n';
      }
      if (receiverCount > 1) {
        for (int receiver = 0; receiver < receiverCount; ++receiver) {
          for (int component = 0; component < 3; ++component) {
            const auto value = raw(3 * receiver + component, index) * unit * timeNorm;
            allSpectra << index << ',' << fHz << ',' << realOmega << ','
                       << imagOmega << ',' << receiver << ',' << component << ','
                       << value.real() << ',' << value.imag() << '\n';
          }
        }
      }
    }
    series << "time_s,Z_SI,N_SI,E_SI\n" << std::setprecision(17);
    const auto outputSamples = std::min<int>(
        times.cols(), static_cast<int>(std::floor(
                         params.inputParameters().t_out() * 60.0 /
                         params.inputParameters().time_step_sec())) + 1);
    for (int index = 0; index < outputSamples; ++index) {
      series << index * grid.dt() * timeNorm;
      for (int component = 0; component < 3; ++component) {
        series << ',' << times(component, index);
      }
      series << '\n';
      if (receiverCount > 1) {
        allSeries << index * grid.dt() * timeNorm;
        for (int receiver = 0; receiver < receiverCount; ++receiver) {
          for (int component = 0; component < 3; ++component)
            allSeries << ',' << times(3 * receiver + component, index);
        }
        allSeries << '\n';
      }
    }
    const auto elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
            .count();
    meta << std::setprecision(17)
         << "{\n  \"solver\": \"DSpecM1D\",\n"
         << "  \"observable\": \""
         << (params.inputParameters().output_type() == 0 ? "displacement" :
             params.inputParameters().output_type() == 1 ? "velocity" :
             params.inputParameters().output_type() == 2 ? "acceleration" :
             "unknown")
         << " (ZNE)\",\n"
         << "  \"gravity_selector\": "
         << (forceCowling ? "\"cowling\"" : "\"full\"") << ",\n"
         << "  \"mode_type\": " << params.inputParameters().type() << ",\n"
         << "  \"degree_min\": " << params.inputParameters().lmin() << ",\n"
         << "  \"degree_max\": " << params.inputParameters().lmax() << ",\n"
         << "  \"mesh_strategy\": \""
         << ((paper2 || (paperSuite && suiteMesh == "dynamic"))
                 ? "dynamic frequency-chunked multi-SEM"
                 : paperSuite ? "single SEM, paper suite" : "single SEM")
         << "\",\n"
         << "  \"nq\": " << params.nq() << ",\n"
         << "  \"maxstep\": ";
    if (paper2 || (paperSuite && suiteMesh == "dynamic"))
      meta << "null,\n  \"maxstep_source\": \"derived per frequency chunk\",\n"
           << "  \"unused_constructor_maxstep_argument\": "
           << params.maxstep() << ",\n";
    else
      meta << params.maxstep() << ",\n";
    meta << "  \"nskip\": ";
    if (paper2 || (paperSuite && suiteMesh == "dynamic"))
      meta << "null,\n  \"nskip_source\": \"derived from frequency grid\",\n"
           << "  \"unused_constructor_nskip_argument\": "
           << params.nskip() << ",\n";
    else
      meta << params.nskip() << ",\n";
    meta
         << "  \"relative_error\": "
         << params.inputParameters().relative_error() << ",\n"
         << "  \"potential_and_tilt_corrections\": false,\n"
         << "  \"source_time_convention\": \"step moment; no extra filter\",\n"
         << "  \"attenuation\": "
         << (params.inputParameters().attenuation() ? "true" : "false")
         << ",\n"
         << "  \"time_norm_s\": " << timeNorm << ",\n"
         << "  \"length_norm_m\": " << params.earthModel().LengthNorm()
         << ",\n"
         << "  \"reference_period_s\": "
         << params.earthModel().TREF() * timeNorm << ",\n"
         << "  \"damping_e_folds\": 5,\n"
         << "  \"damping_epsilon_rad_s\": " << grid.ep() / timeNorm
         << ",\n"
         << "  \"frequency_first_index\": " << grid.i1() << ",\n"
         << "  \"frequency_last_index_exclusive\": " << grid.i2() << ",\n"
         << "  \"samples\": " << grid.nt() << ",\n"
         << "  \"output_samples\": " << outputSamples << ",\n"
         << "  \"frequency_spacing_hz\": " << grid.df() / timeNorm << ",\n"
         << "  \"sample_interval_s\": " << grid.dt() * timeNorm << ",\n"
         << "  \"time_series_unit_factor_si\": " << unit << ",\n"
         << "  \"spectrum_unit_factor_si_per_hz\": "
         << unit * timeNorm << ",\n"
         << "  \"spectra_filtered\": false,\n"
         << "  \"native_time_series_filtered\": "
         << ((paper2 || paperSuite) ? "true" : "false") << ",\n"
         << "  \"filter_mHz\": [" << params.inputParameters().f11() << ", "
         << params.inputParameters().f12() << ", "
         << params.inputParameters().f21() << ", "
         << params.inputParameters().f22() << "],\n"
         << "  \"runtime_s\": " << elapsed << "\n}\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "DSpecM1D adapter: " << error.what() << '\n';
    return 1;
  }
}

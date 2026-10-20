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
  if (argc != 6 && argc != 7) {
    std::cerr << "usage: dspec_adapter INPUT_PARAMS OUTPUT_PREFIX NQ NSKIP MAXSTEP [paper2]\n";
    return 2;
  }
  try {
    const auto started = std::chrono::steady_clock::now();
    const bool paper2 = argc == 7 && std::string(argv[6]) == "paper2";
    if (argc == 7 && !paper2) {
      throw std::invalid_argument("unknown optional adapter mode");
    }
    InputParametersNew params(argv[1], std::stoi(argv[3]), std::stoi(argv[4]),
                              std::stod(argv[5]));
    SPARSESPEC::SparseFSpec solver;
    const bool forceCowling = params.inputParameters().gravitation() == 1;
    Eigen::MatrixXcd raw;
    if (paper2) {
      // The paper comparison uses the frequency-chunked dynamic multi-SEM API.
      raw = solver.spectra(params);
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
    const Eigen::MatrixXd times = paper2
        ? processfunctions::filtfreq2time(raw * unit, grid, true)
        : (processfunctions::freq2time(raw, grid) * unit).eval();

    std::ofstream spectra(std::string(argv[2]) + ".spectra.csv");
    std::ofstream series(std::string(argv[2]) + ".timeseries.csv");
    std::ofstream meta(std::string(argv[2]) + ".meta.json");
    if (!spectra || !series || !meta) {
      throw std::runtime_error("could not open adapter output");
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
         << (paper2 ? "dynamic frequency-chunked multi-SEM" :
                      "single SEM") << "\",\n"
         << "  \"nq\": " << params.nq() << ",\n"
         << "  \"maxstep\": ";
    if (paper2)
      meta << "null,\n  \"maxstep_source\": \"derived per frequency chunk\",\n"
           << "  \"unused_constructor_maxstep_argument\": "
           << params.maxstep() << ",\n";
    else
      meta << params.maxstep() << ",\n";
    meta << "  \"nskip\": ";
    if (paper2)
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
         << (paper2 ? "true" : "false") << ",\n"
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

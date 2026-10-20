#include <YSpec/Input.hpp>
#include <YSpec/TimeSeries.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "usage: yspecpp_adapter INPUT_JSON OUTPUT_PREFIX\n";
    return 2;
  }
  try {
    const auto started = std::chrono::steady_clock::now();
    const auto input = YSpec::ReadInput(argv[1]);
    const auto &settings = input.settings;
    const auto model = YSpec::Model::ReadDeck(input.model);
    const auto &scales = model.GetScales();
    const auto timeScale = scales.Time();
    const auto quantityName = settings.quantity == YSpec::Quantity::Displacement
        ? "displacement"
        : settings.quantity == YSpec::Quantity::Velocity ? "velocity"
                                                          : "acceleration";
    const auto unit = settings.quantity == YSpec::Quantity::Displacement
        ? scales.Length()
        : settings.quantity == YSpec::Quantity::Velocity
              ? scales.Velocity()
              : scales.Acceleration();
    const auto grid = YSpec::FrequencyGrid(
        settings.fMin / 1000.0 * timeScale,
        settings.fMax / 1000.0 * timeScale,
        settings.timeStep / timeScale,
        settings.length * 60.0 / timeScale);
    const auto raw = YSpec::ComputeSpectra(model, settings, grid);
    const auto filterFrequency = [timeScale](double mHz) {
      return mHz / 1000.0 * timeScale;
    };
    const YSpec::CosineFilter filter{
        filterFrequency(settings.filter[0]),
        filterFrequency(settings.filter[1]),
        filterFrequency(settings.filter[2]),
        filterFrequency(settings.filter[3])};
    const auto geometry = YSpec::Geometry(settings.receivers.at(0),
                                          settings.source);
    const auto c = std::cos(geometry.rotation);
    const auto s = std::sin(geometry.rotation);
    std::vector<std::complex<double>> spectraZne(raw.size());
    for (std::size_t index = 0; index < grid.Count(); ++index) {
      const auto radial = raw[index];
      const auto theta = raw[grid.Count() + index];
      const auto phi = raw[2 * grid.Count() + index];
      spectraZne[index] = radial;
      spectraZne[grid.Count() + index] = -theta * c + phi * s;
      spectraZne[2 * grid.Count() + index] = theta * s + phi * c;
    }

    std::ofstream spectra(std::string(argv[2]) + ".spectra.csv");
    std::ofstream series(std::string(argv[2]) + ".timeseries.csv");
    std::ofstream meta(std::string(argv[2]) + ".meta.json");
    if (!spectra || !series || !meta) {
      throw std::runtime_error("could not open adapter output");
    }
    spectra << "frequency_bin,frequency_hz,omega_real_rad_s,omega_imag_rad_s,component,real,imag\n"
            << std::setprecision(17);
    for (std::size_t index = 0; index < grid.Count(); ++index) {
      const auto fHz = grid.Frequency(index) / timeScale;
      const auto omega = grid.Omega(index) / timeScale;
      for (std::size_t component = 0; component < 3; ++component) {
        const auto value = spectraZne[component * grid.Count() + index] *
                           unit * timeScale;
        const auto bin = static_cast<long long>(grid.First() + index);
        spectra << bin << ',' << fHz << ',' << omega.real() << ','
                << omega.imag() << ','
                << component << ',' << value.real() << ',' << value.imag()
                << '\n';
      }
    }

    const auto seriesNd = YSpec::InverseTransform(grid, filter, raw, 3);
    const auto samples = std::min<std::size_t>(
        grid.Samples(), static_cast<std::size_t>(
                            settings.length * 60.0 / settings.timeStep + 1e-9) +
                            1);
    series << "time_s,Z_SI,N_SI,E_SI\n" << std::setprecision(17);
    for (std::size_t index = 0; index < samples; ++index) {
      const auto z = seriesNd[index] * unit;
      const auto theta = seriesNd[grid.Samples() + index] * unit;
      const auto phi = seriesNd[2 * grid.Samples() + index] * unit;
      series << index * grid.TimeStep() * timeScale << ',' << z << ',' << -theta * c + phi * s
             << ',' << theta * s + phi * c << '\n';
    }

    const auto elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
            .count();
    meta << std::setprecision(17)
         << "{\n  \"solver\": \"YSpecpp\",\n"
         << "  \"observable\": \"" << quantityName << " (ZNE)\",\n"
         << "  \"attenuation\": " << (settings.attenuation ? "true" : "false")
         << ",\n"
         << "  \"gravity_selector\": \""
         << (settings.gravity == YSpec::Gravity::Cowling ? "cowling" :
             settings.gravity == YSpec::Gravity::None ? "none" : "full")
         << "\",\n"
         << "  \"degree_min\": " << settings.lMin << ",\n"
         << "  \"degree_max\": " << settings.lMax << ",\n"
         << "  \"solver_tolerance\": "
         << settings.solver.integration.tolerance << ",\n"
         << "  \"solver_start_decay\": " << settings.solver.startDecay << ",\n"
         << "  \"potential_and_tilt_corrections\": false,\n"
         << "  \"source_time_convention\": \"step moment; no extra filter\",\n"
         << "  \"time_norm_s\": " << timeScale << ",\n"
         << "  \"length_norm_m\": " << scales.Length() << ",\n"
         << "  \"reference_period_s\": "
         << ((2.0 * std::numbers::pi * timeScale) / model.ReferenceFrequency()) << ",\n"
         << "  \"damping_e_folds\": 5,\n"
         << "  \"damping_epsilon_rad_s\": "
         << grid.Damping() / timeScale << ",\n"
         << "  \"frequency_first_index\": " << grid.First() << ",\n"
         << "  \"frequency_last_index_inclusive\": " << grid.Last() << ",\n"
         << "  \"samples\": " << grid.Samples() << ",\n"
         << "  \"output_samples\": " << samples << ",\n"
         << "  \"frequency_spacing_hz\": " << grid.Spacing() / timeScale << ",\n"
         << "  \"sample_interval_s\": " << grid.TimeStep() * timeScale
         << ",\n"
         << "  \"time_series_unit_factor\": " << unit << ",\n"
         << "  \"spectrum_unit_factor_si_per_hz\": "
         << unit * timeScale << ",\n"
         << "  \"spectra_filtered\": false,\n"
         << "  \"native_time_series_filtered\": true,\n"
         << "  \"filter_mHz\": [" << settings.filter[0] << ", "
         << settings.filter[1] << ", " << settings.filter[2] << ", "
         << settings.filter[3] << "],\n"
         << "  \"runtime_s\": " << elapsed << "\n}\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "YSpecpp adapter: " << error.what() << '\n';
    return 1;
  }
}

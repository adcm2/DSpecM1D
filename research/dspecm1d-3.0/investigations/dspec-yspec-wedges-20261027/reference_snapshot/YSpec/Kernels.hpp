#pragma once

/**
 * @file Kernels.hpp
 * @brief Sensitivity kernels of measurements made on seismograms, by the
 * adjoint method.
 *
 * @details A measurement @f$J@f$ on one component of the seismogram of a
 * source at a receiver is linearised as @f$\delta J = \int j(t)\,\delta
 * u(t)\,dt@f$, with @f$j@f$ its adjoint source and @f$u@f$ the seismogram as
 * ComputeSeismograms gives it. Its sensitivity to the model is formed from two
 * fields at each of the complex frequencies of the calculation: the forward
 * field of the source, and the adjoint field, which is the response at the
 * same complex frequency to a force at the receiver along the component
 * measured, with a spectrum formed from @f$j@f$. No time reversal is needed.
 * The derivation, from the weak form of the equations of motion with the
 * equilibrium gravitational potential as a variable constrained by Poisson's
 * equation, is in docs/theory/adjoint-kernels.tex.
 *
 * The radial kernels give the measurement's sensitivity to spherically
 * symmetric perturbations of the model: of its density, at fixed moduli and
 * with the equilibrium gravity that follows from it, of its five Love moduli,
 * at fixed density, and of their inverse quality factors. They are evaluated
 * on a mesh in radius, with Gauss-Legendre nodes in elements between the
 * discontinuities of the model, the source and the receiver, and come with
 * the weights of the quadrature on that mesh.
 */

#include <array>
#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include "Calculation.hpp"
#include "Location.hpp"
#include "Model.hpp"
#include "Receiver.hpp"
#include "Seismograms.hpp"
#include "Source.hpp"

namespace YSpec {

/**
 * @brief A window in time, with cosine tapers inside it at each end.
 */
struct Window {
  double start;      ///< Its start (s), after the origin time.
  double end;        ///< Its end (s).
  double taper = 0;  ///< The length of each taper (s).

  /** @brief The value of the window at a time (s). */
  double operator()(double t) const;
};

/**
 * @brief A direction at a receiver, by its components up, north and east.
 */
struct Direction {
  double up = 0;     ///< The upward component.
  double north = 0;  ///< The northward component.
  double east = 0;   ///< The eastward component.
};

/** @brief The direction of a component of a seismogram. */
Direction Along(Component component);

/**
 * @brief The horizontal direction at a receiver away from a source, along the
 * great circle through them.
 */
Direction RadialDirection(const Location& source, const Location& receiver);

/**
 * @brief The horizontal direction at a receiver 90 degrees anticlockwise from
 * RadialDirection, seen from above, as for the transverse component of the
 * snapshots.
 */
Direction TransverseDirection(const Location& source, const Location& receiver);

/**
 * @brief A measurement on a seismogram, along a direction, through its
 * adjoint source.
 *
 * @details The seismogram measured is @f$u = \hat{\mathbf{d}}\cdot
 * \mathbf{u}@f$, the combination of the vertical, north and east components
 * that the direction gives, as ComputeSeismograms gives them for the same
 * calculation, recorded at the receiver by its instrument, if any. The
 * measurement changes by @f$\delta J = \int j(t)\,\delta u(t)\,dt@f$,
 * approximated by the sum over the samples of the seismogram.
 */
struct Measurement {
  Receiver receiver;    ///< Where it is made.
  Direction direction;  ///< The direction of the seismogram measured.
  /// @f$j@f$ at the samples of the seismogram, in units of the measurement
  /// per unit of the seismogram per second.
  std::vector<double> adjointSource;
};

/**
 * @brief The traveltime of a phase in a window, measured by cross-correlation
 * with the synthetic itself.
 *
 * @details A delay @f$\delta\tau@f$ changes the seismogram by
 * @f$-\dot u\,\delta\tau@f$, so to first order @f$\delta\tau = -\int W\dot
 * u\,\delta u\,dt / \int W\dot u^2\,dt@f$ (Dahlen et al. 2000), and the adjoint
 * source is @f$j = -W\dot u / \int W\dot u^2\,dt@f$. The rate @f$\dot u@f$
 * is formed from the spectrum of the synthetic, and is exact.
 * @param model The earth model.
 * @param calculation The calculation of the seismogram.
 * @param source The source.
 * @param receiver The receiver.
 * @param direction The direction of the seismogram measured.
 * @param window The window.
 * @throws std::invalid_argument if the synthetic is zero in the window, or
 * as ComputeSeismograms does.
 */
Measurement Traveltime(const Model& model, const Calculation& calculation,
                       const Source& source, const Receiver& receiver,
                       const Direction& direction, const Window& window);

/**
 * @brief The traveltimes of several phases, or components, at one receiver,
 * as Traveltime gives each, from one synthetic.
 * @param model The earth model.
 * @param calculation The calculation of the seismogram.
 * @param source The source.
 * @param receiver The receiver.
 * @param directions The direction of each seismogram measured.
 * @param windows The window of each.
 * @throws std::invalid_argument if there is not a window for each direction,
 * or as Traveltime does.
 */
std::vector<Measurement> Traveltimes(const Model& model,
                                     const Calculation& calculation,
                                     const Source& source,
                                     const Receiver& receiver,
                                     std::span<const Direction> directions,
                                     std::span<const Window> windows);

/**
 * @brief The parameters of a spherically symmetric model, as the radial
 * kernels give them.
 *
 * @details The inverse quality factors are those of the model: @c QA
 * attenuates A and C, @c QF attenuates F, and @c QMu attenuates L and N.
 */
enum class RadialParameter : std::size_t {
  Rho,  ///< Density, at fixed moduli.
  A,    ///< Love modulus A, at fixed density.
  C,    ///< Love modulus C.
  F,    ///< Love modulus F.
  L,    ///< Love modulus L.
  N,    ///< Love modulus N.
  QA,   ///< Inverse quality factor of A and C.
  QF,   ///< Inverse quality factor of F.
  QMu,  ///< Inverse quality factor of L and N.
  Count
};

/** @brief How finely the radial kernels are evaluated. */
struct RadialMesh {
  int order = 5;             ///< Gauss-Legendre nodes in each element.
  double maxElement = 50e3;  ///< The largest element (m).
};

/**
 * @brief Kernels for spherically symmetric perturbations of the model.
 *
 * @details For each parameter @f$X@f$,
 *
 * @f[ \delta J = \int K_X(r)\,\delta X(r)\,r^2\,dr
 *     + P_X\,\delta X(r_s), @f]
 *
 * the integral over the model, which the weights approximate,
 * @f$\sum_i w_i K_X(r_i)\,\delta X(r_i)@f$, and the second term a point part
 * at the radius of the source. A moment tensor makes the displacement jump
 * across that radius, by @f$\mathbf{M}\cdot\hat{\mathbf{r}}@f$ divided by
 * C or L there, so its strain has a part concentrated on it, and the kernels
 * of C, F and L, and of their quality factors, a part @f$P_X@f$ there; the
 * others have none.
 *
 * In a fluid only perturbations that keep it isotropic have meaning:
 * @f$\delta A = \delta C = \delta F@f$ and @f$\delta q_A = \delta q_F@f$, to
 * which the sum of the kernels responds. Those of L, N and their quality
 * factor are zero there.
 *
 * In SI units: a kernel is in units of the measurement per unit of its
 * parameter per cubic metre, a weight in cubic metres, and a point part in
 * units of the measurement per unit of its parameter.
 */
class RadialKernels {
 public:
  /**
   * @brief Zero kernels at the given radii (m), with their weights (m^3),
   * for a source at radius @p sourceRadius (m).
   */
  RadialKernels(std::vector<double> radii, std::vector<double> weights,
                double sourceRadius);

  /** @brief The number of radii. */
  std::size_t Size() const { return radii_.size(); }
  /** @brief The radii (m), from the centre. */
  std::span<const double> Radii() const { return radii_; }
  /** @brief The weights of the quadrature over radius (m^3). */
  std::span<const double> Weights() const { return weights_; }

  /** @brief The kernel of a parameter, at the radii. */
  std::span<double> operator[](RadialParameter p) {
    return values_[static_cast<std::size_t>(p)];
  }
  /** @brief The kernel of a parameter, at the radii. */
  std::span<const double> operator[](RadialParameter p) const {
    return values_[static_cast<std::size_t>(p)];
  }

  /** @brief The radius of the source (m). */
  double SourceRadius() const { return sourceRadius_; }
  /** @brief The point part of the kernel of a parameter at the source. */
  double& AtSource(RadialParameter p) {
    return atSource_[static_cast<std::size_t>(p)];
  }
  /** @brief The point part of the kernel of a parameter at the source. */
  double AtSource(RadialParameter p) const {
    return atSource_[static_cast<std::size_t>(p)];
  }

 private:
  static constexpr auto Parameters =
      static_cast<std::size_t>(RadialParameter::Count);
  std::vector<double> radii_;
  std::vector<double> weights_;
  std::array<std::vector<double>, Parameters> values_;
  double sourceRadius_;
  std::array<double, Parameters> atSource_{};
};

/**
 * @brief Computes the radial kernels of a measurement.
 * @param model The earth model.
 * @param calculation The calculation of the seismogram measured.
 * @param source The source.
 * @param measurement The measurement.
 * @param mesh How finely to evaluate the kernels.
 * @param progress Called as frequencies are completed, if given.
 * @throws std::invalid_argument if the adjoint source does not have a value
 * at every sample of the seismogram, if the receiver is not in the solid
 * region of the source (the adjoint field's source is there), if the
 * seismogram has the potential and tilt corrections, or as
 * ComputeSeismograms does.
 */
RadialKernels ComputeRadialKernels(const Model& model,
                                   const Calculation& calculation,
                                   const Source& source,
                                   const Measurement& measurement,
                                   const RadialMesh& mesh = {},
                                   const Progress& progress = {});

/**
 * @brief Computes the radial kernels of several measurements at one
 * receiver, for little more than the cost of one: the forward field is the
 * same for all, and the adjoint field is linear in its force.
 * @throws std::invalid_argument if the measurements are not all at one
 * receiver, or as the other overload does.
 */
std::vector<RadialKernels> ComputeRadialKernels(
    const Model& model, const Calculation& calculation, const Source& source,
    std::span<const Measurement> measurements, const RadialMesh& mesh = {},
    const Progress& progress = {});

/** @brief The sets of parameters in which radial kernels can be given. */
enum class ParameterSet {
  /// rho, A, C, F, L, N, qA, qF, qMu: the model's own, density at fixed
  /// moduli and the moduli at fixed density.
  Love,
  /// rho, kappa, mu, qkappa, qmu: density at fixed moduli, and the Voigt
  /// averages of the moduli, perturbed isotropically: A and C change by
  /// @f$\delta\kappa + \tfrac43\delta\mu@f$, F by @f$\delta\kappa -
  /// \tfrac23\delta\mu@f$, L and N by @f$\delta\mu@f$.
  Moduli,
  /// rho, vp, vs, qkappa, qmu: density at fixed velocities, and the speeds
  /// of the Voigt averages, @f$v_p^2 = (\kappa + \tfrac43\mu)/\rho@f$ and
  /// @f$v_s^2 = \mu/\rho@f$, each at fixed density and the other speed.
  Velocities,
  /// rho, vpv, vph, vsv, vsh, eta, qkappa, qmu: density at fixed velocities
  /// and eta, and the five parameters of transverse isotropy, with
  /// @f$A = \rho v_{ph}^2@f$, @f$C = \rho v_{pv}^2@f$, @f$L = \rho
  /// v_{sv}^2@f$, @f$N = \rho v_{sh}^2@f$ and @f$F = \eta(A - 2L)@f$.
  Anisotropic,
};

/** @brief Kernels of perturbations @f$\delta X@f$, or of @f$\delta\ln X@f$. */
enum class Scaling { Absolute, Relative };

/**
 * @brief Radial kernels of a set of parameters, by name.
 *
 * @details As RadialKernels: @f$\delta J = \int\sum_X K_X\,\delta X\,
 * r^2\,dr + \sum_X P_X\,\delta X(r_s)@f$, with @f$\delta X@f$ replaced by
 * @f$\delta\ln X@f$ for relative kernels, which are @f$X K_X@f$ and
 * @f$X P_X@f$. The inverse quality factors qkappa and qmu reach the model's
 * through Minos's effective quality factors of A, C and F.
 */
struct NamedKernels {
  std::vector<std::string> names;           ///< The parameters.
  std::vector<std::vector<double>> values;  ///< Each kernel, at the radii.
  std::vector<double> atSource;             ///< Each point part.
};

/**
 * @brief Radial kernels in a set of parameters.
 * @param kernels The kernels in the model's own parameters.
 * @param model The model they are for.
 * @param set The parameters wanted.
 * @param scaling Absolute or relative.
 */
NamedKernels Reparameterise(const RadialKernels& kernels, const Model& model,
                            ParameterSet set, Scaling scaling);

/** @brief A radius (m) in a layer of the model, at which kernels are sampled.
 */
struct KernelNode {
  Int layer;      ///< The layer, which settles the side of a discontinuity.
  double radius;  ///< The radius (m).
};

/** @brief A place on the sphere, at which kernels are sampled. */
struct KernelPlace {
  double colatitude;  ///< The geocentric colatitude (radians).
  double longitude;   ///< The longitude (radians).
  /// For places laid out by the path from the source to the receiver, the
  /// distance along its great circle from the source (degrees).
  double along = 0;
  /// And the distance from it (degrees), positive to the left of travel.
  double across = 0;
};

/** @brief Where kernels are sampled: at each node, the same places. */
struct KernelSampling {
  std::vector<KernelNode> nodes;    ///< The nodes in radius.
  std::vector<KernelPlace> places;  ///< The places on the sphere.
};

/** @brief How three-dimensional kernels are sampled through the earth. */
struct VolumeGrid {
  int degree = 64;            ///< The angular grid resolves degrees to this.
  double maxElement = 100e3;  ///< The largest interval between radii (m).
};

/**
 * @brief The whole earth, laid out as planetmodel lays out its model files:
 * in each layer, radii equally spaced from its bottom to its top, both
 * included, so that each side of a discontinuity has its own; at each, a
 * Gauss-Legendre grid of geocentric colatitudes, increasing, and equally
 * spaced longitudes from zero, colatitude by colatitude, followed by the
 * north and south poles.
 */
KernelSampling VolumeSampling(const Model& model, const VolumeGrid& grid);

/**
 * @brief How sections along and across the path from the source to the
 * receiver are laid out. Angles are in degrees.
 */
struct SectionOptions {
  double margin = 15;           ///< How far the path runs beyond its ends.
  double halfWidth = 15;        ///< How far a section across reaches from it.
  double spacing = 0.2;         ///< The spacing of places on the sections.
  double radialSpacing = 25e3;  ///< The largest interval between radii (m).
};

/**
 * @brief The section through the earth along the great circle through the
 * source and the receiver, from @c margin behind the source to as far beyond
 * the receiver, with radii through every layer as in VolumeSampling.
 */
KernelSampling PathSection(const Model& model, const Location& source,
                           const Location& receiver,
                           const SectionOptions& options = {});

/**
 * @brief The section through the earth across the path at its middle, along
 * the great circle perpendicular to it there, out to @c halfWidth on each
 * side.
 */
KernelSampling CrossSection(const Model& model, const Location& source,
                            const Location& receiver,
                            const SectionOptions& options = {});

/** @brief How maps around the path are laid out. Angles are in degrees. */
struct MapOptions {
  std::vector<double> depths;  ///< The depths of the maps (m).
  double margin = 10;          ///< How far they run beyond its ends.
  double halfWidth = 15;       ///< How far they reach from it.
  double spacing = 0.5;        ///< The spacing of places on them.
};

/**
 * @brief Maps at some depths of the region around the path, on a grid of
 * distances along and across it. A depth at a discontinuity is taken on its
 * upper side.
 */
KernelSampling PathMaps(const Model& model, const Location& source,
                        const Location& receiver, const MapOptions& options);

/**
 * @brief Kernels for laterally varying perturbations of the elastic
 * parameters, at fixed density, sampled at points.
 *
 * @details @f$\delta J = \int_V \sum_X K_X\,\delta X\,dV@f$, with
 * @f$\delta X@f$ the perturbations of the elastic parameters of a set and of
 * their inverse quality factors (or of their logarithms, for relative
 * kernels), at fixed density: perturbations of any shape that leave the
 * model in hydrostatic equilibrium. Density, whose lateral perturbations do
 * not, has none (see docs/theory/adjoint-kernels.tex).
 *
 * Each field is evaluated in its own frame, the forward field's with the
 * source at the pole and the adjoint field's with the receiver there, so
 * that only the orders a point source excites are summed, and the two are
 * turned into one frame at each point. At the source the kernels also have a
 * part concentrated there, the three-dimensional form of the point parts of
 * RadialKernels, which no sampling holds; and near the source and the
 * receiver they vary on the scale of the shortest wavelength summed.
 */
class SampledKernels {
 public:
  /** @brief Zero kernels of the named parameters on a sampling. */
  SampledKernels(std::vector<std::string> names, KernelSampling sampling);

  /** @brief The parameters. */
  const std::vector<std::string>& Names() const { return names_; }
  /** @brief The nodes in radius. */
  std::span<const KernelNode> Nodes() const { return sampling_.nodes; }
  /** @brief The places on the sphere. */
  std::span<const KernelPlace> Places() const { return sampling_.places; }

  /** @brief The kernel of a parameter at the places of a node. */
  std::span<double> Values(std::size_t parameter, std::size_t node) {
    return {values_.data() + Offset(parameter, node), Places().size()};
  }
  /** @brief The kernel of a parameter at the places of a node. */
  std::span<const double> Values(std::size_t parameter,
                                 std::size_t node) const {
    return {values_.data() + Offset(parameter, node), Places().size()};
  }

 private:
  std::vector<std::string> names_;
  KernelSampling sampling_;
  std::vector<double> values_;

  std::size_t Offset(std::size_t parameter, std::size_t node) const {
    return (parameter * Nodes().size() + node) * Places().size();
  }
};

/**
 * @brief Computes the three-dimensional kernels of several measurements at
 * one receiver on several samplings, in one pass.
 *
 * @details The forward field is the same for every measurement, and the
 * adjoint field is linear in its force: it is evaluated for an orthonormal
 * basis of the measurements' directions, one, two or three forces, and each
 * measurement's is a combination, so that further measurements cost little.
 * @param model The earth model.
 * @param calculation The calculation of the seismogram measured.
 * @param source The source.
 * @param measurements The measurements, all at one receiver.
 * @param set The parameters, of which density is left out.
 * @param scaling Absolute or relative.
 * @param samplings Where to sample the kernels.
 * @param progress Called as frequencies are completed, if given.
 * @return The kernels of each measurement on each sampling.
 * @throws std::invalid_argument as ComputeRadialKernels does.
 */
std::vector<std::vector<SampledKernels>> ComputeSampledKernels(
    const Model& model, const Calculation& calculation, const Source& source,
    std::span<const Measurement> measurements, ParameterSet set,
    Scaling scaling, std::span<const KernelSampling> samplings,
    const Progress& progress = {});

/**
 * @brief Computes the three-dimensional kernels of a measurement through the
 * whole earth, on VolumeSampling.
 * @throws std::invalid_argument as ComputeRadialKernels does, or if the grid
 * is unusable.
 */
SampledKernels ComputeVolumeKernels(
    const Model& model, const Calculation& calculation, const Source& source,
    const Measurement& measurement, ParameterSet set, Scaling scaling,
    const VolumeGrid& grid = {}, const Progress& progress = {});

}  // namespace YSpec

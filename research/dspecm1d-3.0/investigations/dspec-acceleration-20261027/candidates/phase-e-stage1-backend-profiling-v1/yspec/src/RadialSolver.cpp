#include <YSpec/RadialSolver.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

#include "Equations.hpp"
#include "RungeKutta.hpp"

// The radial equations as a two-point boundary value problem.
//
// For every system the same steps are taken. A solution regular at the bottom
// is integrated upwards from a starting level, and one satisfying the
// conditions at the top downwards from the top; each is kept at every source
// and receiver it passes. At a source the two give the determinant of the
// problem. Above the source the response is formed from the solution from
// above, and below it from the solution from below, each combined with the
// jump at the source. For a 2x2 system the combination is a constant; for the
// spheroidal systems, integrated as minors, it is a b-vector carried from the
// source to each receiver (Al-Attar & Woodhouse 2008).
//
// The systems differ in how they start, what they integrate and how the
// pieces combine, and each is a small class below with the same members.

namespace YSpec {

namespace {

//----------------------------------------------------------------------------
// The spheroidal equations, in the three treatments of gravity.
//----------------------------------------------------------------------------

// Without self-gravitation: 4x4 in a solid, 2x2 in a fluid.
template <bool Cowling>
struct WithoutPotential {
  static constexpr std::size_t SolidSize = 5;  // second-order minors
  static constexpr std::size_t FluidSize = 2;  // the solution itself
  static constexpr std::size_t BSize = 4;
  // In a fluid the b-vector is a constant, which needs no integration.
  static constexpr std::size_t FluidBSize = 1;
  static constexpr bool IntegrateFluidB = false;
  // The sign of the response above the source.
  static constexpr double SignAbove = 1;

  static Block4 SolidBlock(const Context& c, const ComplexMaterial& m) {
    return SolidBlock4(c, m, Cowling);
  }
  static Block2 FluidBlock(const Context& c, const ComplexMaterial& m) {
    return FluidBlock2(c, m, Cowling);
  }
  static Vector<5> SolidStart(const Context& c, const ComplexMaterial& m) {
    return SolidStart4(c, m);
  }
  static Vector<2> FluidStart(const Context& c, const ComplexMaterial& m) {
    return FluidStart2(c, m);
  }
  static Vector<5> SolidDerivative(const Block4& block,
                                   const Vector<5>& minors) {
    return MinorDerivative(block, minors);
  }
  static Vector<2> FluidDerivative(const Block2& block, const Vector<2>& y) {
    return Derivative(block, y);
  }
  // A unit jump in (U, V, P, S).
  static Vector<4> Jump(int j) {
    auto jump = Vector<4>{};
    jump[j] = 1;
    return jump;
  }
  static Vector<4> BVector(const Vector<5>& minors, const Vector<4>& jump) {
    return BStart(minors, jump);
  }
  // U, V, the potential, the flux and the tractions P and S, from the
  // solution vector.
  static std::array<Complex, 6> Components(const Vector<4>& y) {
    return {y[0], y[1], Complex{}, Complex{}, y[2], y[3]};
  }
  // The same in a fluid, from the solution and the constant.
  static std::array<Complex, 6> FluidComponents(const Vector<2>& solution,
                                                const Vector<1>& b,
                                                const Context& c,
                                                const ComplexMaterial& m) {
    const auto rU = solution[0] * b[0];
    const auto P = solution[1] * b[0];
    return {rU,        FluidConsoidal(c, m, rU, 0.0, P, Cowling),
            Complex{}, Complex{},
            P,         Complex{}};
  }
};

// With self-gravitation: 6x6 in a solid, 4x4 in a fluid.
struct WithPotential {
  static constexpr std::size_t SolidSize = 14;  // third-order minors
  static constexpr std::size_t FluidSize = 5;   // second-order minors
  static constexpr std::size_t BSize = 14;
  static constexpr std::size_t FluidBSize = 4;
  static constexpr bool IntegrateFluidB = true;
  static constexpr double SignAbove = -1;

  static Block6 SolidBlock(const Context& c, const ComplexMaterial& m) {
    return SolidBlock6(c, m);
  }
  static Block4 FluidBlock(const Context& c, const ComplexMaterial& m) {
    return FluidBlock4(c, m);
  }
  static Vector<14> SolidStart(const Context& c, const ComplexMaterial& m) {
    return SolidStart6(c, m);
  }
  static Vector<5> FluidStart(const Context& c, const ComplexMaterial& m) {
    return FluidStart4(c, m);
  }
  static Vector<14> SolidDerivative(const Block6& block,
                                    const Vector<14>& minors) {
    return MinorDerivative(block, minors);
  }
  static Vector<5> FluidDerivative(const Block4& block,
                                   const Vector<5>& minors) {
    return MinorDerivative(block, minors);
  }
  // A unit jump in (U, V, P, S), placed in (U, V, phi, P, S, psi): the
  // potential and its flux are continuous at the source.
  static Vector<6> Jump(int j) {
    constexpr auto position = std::array{0, 1, 3, 4};
    auto jump = Vector<6>{};
    jump[position[j]] = 1;
    return jump;
  }
  static Vector<14> BVector(const Vector<14>& minors, const Vector<6>& jump) {
    return BStart(minors, jump);
  }
  static std::array<Complex, 6> Components(const Vector<6>& y) {
    return {y[0], y[1], y[2], y[5], y[3], y[4]};
  }
  // The same in a fluid, whose variables are (U, phi, P, psi).
  static std::array<Complex, 6> FluidComponents(const Vector<5>& minors,
                                                const Vector<4>& b,
                                                const Context& c,
                                                const ComplexMaterial& m) {
    const auto y = Solution(minors, b);
    return {y[0], FluidConsoidal(c, m, y[0], y[1], y[2], true),
            y[1], y[3],
            y[2], Complex{}};
  }
};

// Minors carried through regions of either phase.
template <class Form>
struct Minors {
  ScaledVector<Form::SolidSize> solid;
  ScaledVector<Form::FluidSize> fluid;
  bool inFluid = false;

  void ChangePhase(bool toFluid) {
    if (toFluid == inFluid) return;
    if (toFluid) {
      fluid = {SolidToFluid(solid.value), solid.logScale};
      fluid.Normalise();
    } else {
      solid = {FluidToSolid(fluid.value), fluid.logScale};
      solid.Normalise();
    }
    inFluid = toFluid;
  }
};

// The b-vectors of the four unit jumps, carried through regions of either
// phase.
template <class Form>
struct BVectors {
  std::array<ScaledVector<Form::BSize>, 4> solid;
  std::array<ScaledVector<Form::FluidBSize>, 4> fluid;
  bool inFluid = false;

  void ChangePhase(bool toFluid) {
    if (toFluid == inFluid) return;
    for (auto j : {0, 1, 2, 3}) {
      if (toFluid) {
        fluid[j] = {BIntoFluid(solid[j].value), solid[j].logScale};
        fluid[j].Normalise();
      } else {
        solid[j] = {BIntoSolid(fluid[j].value), fluid[j].logScale};
        solid[j].Normalise();
      }
    }
    inFluid = toFluid;
  }
};

// A point at which a state is kept, and the slot it is kept in.
struct Stop {
  std::size_t slot;
  Point point;
};

}  // namespace

//----------------------------------------------------------------------------
// The solver.
//----------------------------------------------------------------------------

class RadialSolver::Impl {
 public:
  Impl(const Model& model, std::vector<double> sourceRadii,
       std::vector<Point> receivers, Gravity gravity, bool attenuation,
       SolverOptions options)
      : model_{&model},
        gravity_{gravity},
        attenuation_{attenuation},
        options_{options},
        integrator_{IntegrationOptions{.tolerance = options.tolerance}} {
    if (!(options.tolerance > 0 && options.tolerance < 1)) {
      throw std::invalid_argument("the tolerance must lie between 0 and 1");
    }
    if (!(options.startDecay > 0)) {
      throw std::invalid_argument("the start decay must be positive");
    }
    if (!(options.stepFactor > 0)) {
      throw std::invalid_argument("the step factor must be positive");
    }
    if (sourceRadii.empty() || receivers.empty()) {
      throw std::invalid_argument("there must be a source and a receiver");
    }
    for (auto r : sourceRadii) {
      if (!(r > 0)) {
        throw std::invalid_argument(
            "a source at the centre of the model is not supported");
      }
      sources_.push_back(SolidSide(model, r));
    }
    for (const auto& receiver : receivers) {
      if (!(receiver.r > 0)) {
        throw std::invalid_argument(
            "a receiver at the centre of the model is not supported");
      }
      if (receiver.layer < 0 || receiver.layer >= model.NumberOfLayers()) {
        throw std::invalid_argument("a receiver is in no layer of the model");
      }
      const auto& layer = model.GetLayer(receiver.layer);
      const auto tolerance = 1e-12 * model.SurfaceRadius();
      if (receiver.r < layer.Lower() - tolerance ||
          receiver.r > layer.Upper() + tolerance) {
        throw std::invalid_argument("a receiver is outside its layer");
      }
    }
    receivers_ = std::move(receivers);

    region_ = model.RegionIndex(sources_.front().layer);
    for (const auto& source : sources_) {
      if (model.GetLayer(source.layer).IsFluid()) {
        throw std::invalid_argument(
            "a source is in a fluid region of the model (an ocean?); it must "
            "be in a solid");
      }
    }
    for (const auto& source : sources_) {
      if (model.RegionIndex(source.layer) != region_) {
        throw std::invalid_argument(
            "the sources must lie in one solid region of the model");
      }
    }
    for (const auto& layer : model.Layers()) {
      auto alpha = 0.0;
      for (auto r : layer.Knots()) alpha = std::max(alpha, layer(r).Alpha());
      maxAlpha_.push_back(alpha);
    }
  }

  // At a boundary between solid and fluid, the point on the solid side.
  static Point SolidSide(const Model& model, double r) {
    const auto below = model.LayerIndex(r, Side::Below);
    const auto above = model.LayerIndex(r, Side::Above);
    if (model.GetLayer(below).IsFluid() && !model.GetLayer(above).IsFluid()) {
      return {above, r};
    }
    return {below, r};
  }

  const Model& GetModel() const { return *model_; }
  const std::vector<Point>& Sources() const { return sources_; }
  const std::vector<Point>& Receivers() const { return receivers_; }
  Gravity GetGravity() const { return gravity_; }
  bool Attenuation() const { return attenuation_; }
  const RungeKutta<DormandPrince54>& Integrator() const { return integrator_; }
  Int SourceRegionIndex() const { return region_; }
  const Region& SourceRegion() const { return model_->Regions()[region_]; }

  // Every source and receiver, sources first.
  std::vector<Stop> All() const {
    auto all = std::vector<Stop>{};
    for (auto i : std::views::iota(std::size_t{0}, sources_.size())) {
      all.push_back({i, sources_[i]});
    }
    for (auto i : std::views::iota(std::size_t{0}, receivers_.size())) {
      all.push_back({sources_.size() + i, receivers_[i]});
    }
    return all;
  }

  Point StartPoint(int l, double omega, Point top, Int lowestRegion) const {
    const auto& lowest = model_->Regions()[lowestRegion];
    const auto bottom = Point{lowest.first, lowest.lower};
    if (l == 0 || omega <= 0) return bottom;
    // The square of the ray parameter, l(l+1) / omega^2.
    const auto rayParameter2 =
        static_cast<double>(l) * (l + 1) / (omega * omega);
    // The square of the vertical slowness of the slowest wave, negative
    // where the solution is evanescent.
    auto slowness2 = [&](Int layer, double r) {
      if (r <= 0) return -std::numeric_limits<double>::infinity();
      const auto& here = model_->GetLayer(layer);
      const auto material = here(r);
      const auto speed = here.IsFluid() ? material.Alpha() : material.Beta();
      return 1 / (speed * speed) - rayParameter2 / (r * r);
    };
    auto decay = 0.0;
    auto region = model_->RegionIndex(top.layer);
    for (auto layer = top.layer; layer >= lowest.first; --layer) {
      if (model_->RegionIndex(layer) != region) {
        region = model_->RegionIndex(layer);
        decay = 0;
      }
      const auto knots = model_->GetLayer(layer).Knots();
      for (auto k = std::ssize(knots) - 1; k > 0; --k) {
        const auto upper = knots[k];
        const auto lower = knots[k - 1];
        if (lower >= top.r) continue;
        const auto qUpper = slowness2(layer, upper);
        const auto qLower = slowness2(layer, lower);
        if (qUpper < 0 && qLower < 0) {
          // The trapezium rule for the WKB decay, omega times the
          // evanescent slowness, over the segment.
          const auto rateUpper =
              std::isinf(qUpper) ? std::sqrt(-qLower) : std::sqrt(-qUpper);
          const auto rateLower =
              std::isinf(qLower) ? rateUpper : std::sqrt(-qLower);
          decay += 0.5 * omega * (rateUpper + rateLower) *
                   (std::min(upper, top.r) - lower);
        } else {
          decay = 0;
        }
        if (decay > options_.startDecay) return Point{layer, lower};
      }
    }
    return bottom;
  }

  // The first step in a layer.
  double FirstStep(Int layer, Complex omega) const {
    return options_.stepFactor * maxAlpha_[layer] / std::abs(omega);
  }

  // A start at the centre is moved off it, since the starting solutions
  // divide by the radius, though not above the lowest point wanted: the
  // solutions are integrated up from the start.
  Point OffCentre(Point start, Point lowest, Complex omega,
                  double& step) const {
    step = FirstStep(start.layer, omega);
    if (start.r > 0) return start;
    const auto knots = model_->GetLayer(start.layer).Knots();
    return {start.layer,
            std::min({0.1 * step, 0.1 * knots[1], 0.5 * lowest.r})};
  }

  // Calls integrate(layer, begin, end, step) for each layer from one point
  // to another, in order, with the interval of the layer between them.
  template <class IntegrateLayer>
  void WalkLayers(Point from, Point to, double& step, Complex omega,
                  IntegrateLayer&& integrate) const {
    const auto up = to.layer >= from.layer;
    for (auto layer = from.layer;; layer += up ? 1 : -1) {
      const auto& here = model_->GetLayer(layer);
      const auto begin =
          layer == from.layer ? from.r : (up ? here.Lower() : here.Upper());
      const auto end =
          layer == to.layer ? to.r : (up ? here.Upper() : here.Lower());
      if (step <= 0) step = FirstStep(layer, omega);
      integrate(layer, begin, end, step);
      if (layer == to.layer) break;
    }
  }

  // Integrates a state from one point to another, calling keep(slot, state)
  // at each of the stops, which lie between them in the order of travel.
  template <class State, class IntegrateLayer, class Keep>
  void WalkKeeping(Point from, Point to, std::span<const Stop> stops,
                   State& state, double& step, Complex omega,
                   IntegrateLayer&& integrate, Keep&& keep) const {
    auto inLayer = [&](Int layer, double begin, double end, double& step) {
      integrate(state, layer, begin, end, step);
    };
    auto here = from;
    for (const auto& stop : stops) {
      WalkLayers(here, stop.point, step, omega, inLayer);
      keep(stop.slot, state);
      here = stop.point;
    }
    WalkLayers(here, to, step, omega, inLayer);
  }

  // The responses of a system for every source and receiver: for each, the
  // components the system gives for each unit jump, in its scaled variables.
  template <class System>
  auto Solve(const System& system, Complex omega) const;

 private:
  const Model* model_;
  Gravity gravity_;
  bool attenuation_;
  SolverOptions options_;
  RungeKutta<DormandPrince54> integrator_;
  std::vector<Point> sources_;
  std::vector<Point> receivers_;
  Int region_;
  std::vector<double> maxAlpha_;

  // A radius in its layer, taking the solid side of a solid-fluid boundary.
};

namespace {

//----------------------------------------------------------------------------
// The systems.
//
// Each gives: the State it integrates and how (Integrate), where and how the
// integrations from below and from above start (FromBelow, FromAbove), which
// points it reaches (Reaches), the determinant at a source, the b-vectors
// that carry the unit jumps at a source out to the receivers (StartBs, and
// IntegrateBs if they need integrating), and the components of the response
// formed from a state and the b-vectors at a receiver (Solution).
//----------------------------------------------------------------------------

// What the 2x2 systems share. The solution from above is itself the response
// above the source, up to a constant, and likewise below, so the "b-vector"
// is a constant: the Wronskian of the solution at the source with the jump.
struct TwoByTwo {
  using State = ScaledVector<2>;
  using Bs = std::array<ScaledVector<1>, 2>;
  static constexpr int Jumps = 2;
  static constexpr std::size_t Components = 2;
  static constexpr bool IntegrateB = false;
  static constexpr double SignAbove = -1;
  static constexpr double SignBelow = -1;

  static ScaledComplex Determinant(const State& below, const State& above) {
    return {above.value[1] * below.value[0] - above.value[0] * below.value[1],
            below.logScale + above.logScale};
  }
  // For the jumps (1, 0) and (0, 1) at the source, with the solution there
  // from the other side; the same above the source and below it.
  static Bs StartBs(const State& atSource, bool) {
    auto bs = Bs{ScaledVector<1>{{atSource.value[1]}, atSource.logScale},
                 ScaledVector<1>{{-atSource.value[0]}, atSource.logScale}};
    for (auto& b : bs) b.Normalise();
    return bs;
  }
  static void IntegrateBs(Bs&, Int, double, double, double&) {}
  static ScaledVector<2> Solution(const State& state, const Bs& bs, int jump,
                                  Point) {
    const auto& b = bs[jump];
    return {{state.value[0] * b.value[0], state.value[1] * b.value[0]},
            state.logScale + b.logScale};
  }
};

// The toroidal system, confined to the solid region of the sources: elsewhere
// the response is zero.
template <class Solver>
class ToroidalSystem : public TwoByTwo {
 public:
  ToroidalSystem(const Solver& solver, const Context& context)
      : solver_{solver}, context_{context} {}

  void Integrate(State& state, Int layer, double begin, double end,
                 double& step) const {
    const auto& context = context_;
    solver_.Integrator().Integrate(
        [&context, layer](double r, const Vector<2>& y) {
          return Derivative(ToroidalBlock(context, context(layer, r)), y);
        },
        begin, end, state, step);
  }

  bool Reaches(Point point) const {
    return solver_.GetModel().RegionIndex(point.layer) ==
           solver_.SourceRegionIndex();
  }

  // A starting level in the region, below the deepest source.
  Point StartLevel(Point top) const {
    return solver_.StartPoint(context_.Degree(), context_.Omega().real(), top,
                              solver_.SourceRegionIndex());
  }

  // From the starting level, or from the bottom of the region, where a fluid
  // below leaves no traction.
  std::pair<Point, State> FromBelow(Point start, Point lowest,
                                    double& step) const {
    const auto& region = solver_.SourceRegion();
    if (start.r == region.lower && region.lower > 0) {
      return {start, State{{1.0, 0.0}}};
    }
    const auto point = solver_.OffCentre(start, lowest, context_.Omega(), step);
    auto state = State{ToroidalStart(context_, context_(point.layer, point.r))};
    state.Normalise();
    return {point, state};
  }

  // From the top of the region, where the traction vanishes.
  std::pair<Point, State> FromAbove() const {
    const auto& region = solver_.SourceRegion();
    return {Point{region.last, region.upper}, State{{1.0, 0.0}}};
  }

 private:
  const Solver& solver_;
  const Context& context_;
};

// The radial (l = 0) system, through the whole model: in a fluid it is the
// same system with no rigidity, and nothing slips at a boundary.
template <class Solver>
class RadialSystem : public TwoByTwo {
 public:
  RadialSystem(const Solver& solver, const Context& context)
      : solver_{solver},
        context_{context},
        gravity_{solver.GetGravity() != Gravity::None} {}

  void Integrate(State& state, Int layer, double begin, double end,
                 double& step) const {
    const auto& context = context_;
    const auto gravity = gravity_;
    solver_.Integrator().Integrate(
        [&context, layer, gravity](double r, const Vector<2>& y) {
          return Derivative(RadialBlock(context, context(layer, r), gravity),
                            y);
        },
        begin, end, state, step);
  }

  static bool Reaches(Point) { return true; }

  // From the centre.
  Point StartLevel(Point) const {
    return {0, solver_.GetModel().InnerRadius()};
  }
  std::pair<Point, State> FromBelow(Point start, Point lowest,
                                    double& step) const {
    const auto point = solver_.OffCentre(start, lowest, context_.Omega(), step);
    auto state =
        State{RadialStart(context_, context_(point.layer, point.r), gravity_)};
    state.Normalise();
    return {point, state};
  }

  // From the surface, where the traction vanishes.
  std::pair<Point, State> FromAbove() const {
    const auto& model = solver_.GetModel();
    return {Point{model.NumberOfLayers() - 1, model.SurfaceRadius()},
            State{{1.0, 0.0}}};
  }

 private:
  const Solver& solver_;
  const Context& context_;
  bool gravity_;
};

// The spheroidal systems, integrated as minors through the whole model, with
// the b-vectors carried through it too.
template <class Form, class Solver>
class SpheroidalSystem {
 public:
  using State = Minors<Form>;
  using Bs = BVectors<Form>;
  static constexpr int Jumps = 4;
  static constexpr std::size_t Components = 6;
  static constexpr bool IntegrateB = true;
  static constexpr double SignAbove = Form::SignAbove;
  static constexpr double SignBelow = -1;

  SpheroidalSystem(const Solver& solver, const Context& context)
      : solver_{solver}, context_{context} {}

  void Integrate(State& minors, Int layer, double begin, double end,
                 double& step) const {
    const auto& context = context_;
    minors.ChangePhase(IsFluid(layer));
    if (minors.inFluid) {
      solver_.Integrator().Integrate(
          [&context, layer](double r, const Vector<Form::FluidSize>& y) {
            return Form::FluidDerivative(
                Form::FluidBlock(context, context(layer, r)), y);
          },
          begin, end, minors.fluid, step);
    } else {
      solver_.Integrator().Integrate(
          [&context, layer](double r, const Vector<Form::SolidSize>& y) {
            return Form::SolidDerivative(
                Form::SolidBlock(context, context(layer, r)), y);
          },
          begin, end, minors.solid, step);
    }
  }

  static bool Reaches(Point) { return true; }

  // A starting level anywhere below the deepest source.
  Point StartLevel(Point top) const {
    return solver_.StartPoint(context_.Degree(), context_.Omega().real(), top,
                              0);
  }

  // From the starting level.
  std::pair<Point, State> FromBelow(Point start, Point lowest,
                                    double& step) const {
    const auto point = solver_.OffCentre(start, lowest, context_.Omega(), step);
    const auto material = context_(point.layer, point.r);
    auto minors = State{};
    minors.inFluid = IsFluid(point.layer);
    if (minors.inFluid) {
      minors.fluid.value = Form::FluidStart(context_, material);
      minors.fluid.Normalise();
    } else {
      minors.solid.value = Form::SolidStart(context_, material);
      minors.solid.Normalise();
    }
    return {point, minors};
  }

  // From the surface, where the traction vanishes.
  std::pair<Point, State> FromAbove() const {
    const auto& model = solver_.GetModel();
    const auto top = model.NumberOfLayers() - 1;
    auto minors = State{};
    minors.inFluid = IsFluid(top);
    minors.fluid.value[0] = 1;
    minors.solid.value[0] = 1;
    return {Point{top, model.SurfaceRadius()}, minors};
  }

  // The sources are in a solid.
  static ScaledComplex Determinant(const State& below, const State& above) {
    return {SecularDeterminant(below.solid.value, above.solid.value),
            below.solid.logScale + above.solid.logScale};
  }

  // For the four unit jumps at the source, with the minors there from the
  // other side; the same formula above the source and below it.
  static Bs StartBs(const State& atSource, bool) {
    auto bs = Bs{};
    for (auto j : {0, 1, 2, 3}) {
      bs.solid[j] = {Form::BVector(atSource.solid.value, Form::Jump(j)),
                     atSource.solid.logScale};
      bs.solid[j].Normalise();
    }
    return bs;
  }

  // The four b-vectors together: their coefficients are the same. In a fluid
  // without self-gravitation they are constants.
  void IntegrateBs(Bs& bs, Int layer, double begin, double end,
                   double& step) const {
    const auto& context = context_;
    bs.ChangePhase(IsFluid(layer));
    if (!bs.inFluid) {
      solver_.Integrator().IntegrateTogether(
          [&context, layer](double r, const Vectors<Form::BSize, Jumps>& b) {
            const auto block = Form::SolidBlock(context, context(layer, r));
            auto derivatives = Vectors<Form::BSize, Jumps>{};
            for (auto j = std::size_t{0}; j < Jumps; ++j) {
              derivatives[j] = BDerivative(block, b[j]);
            }
            return derivatives;
          },
          begin, end, bs.solid, step);
    } else if constexpr (Form::IntegrateFluidB) {
      solver_.Integrator().IntegrateTogether(
          [&context, layer](double r,
                            const Vectors<Form::FluidBSize, Jumps>& b) {
            const auto block = Form::FluidBlock(context, context(layer, r));
            auto derivatives = Vectors<Form::FluidBSize, Jumps>{};
            for (auto j = std::size_t{0}; j < Jumps; ++j) {
              derivatives[j] = BDerivative(block, b[j]);
            }
            return derivatives;
          },
          begin, end, bs.fluid, step);
    }
  }

  // At a receiver the minors and the b-vectors are in the same phase.
  ScaledVector<6> Solution(const State& minors, const Bs& bs, int jump,
                           Point point) const {
    if (minors.inFluid) {
      const auto& b = bs.fluid[jump];
      return {Form::FluidComponents(minors.fluid.value, b.value, context_,
                                    context_(point.layer, point.r)),
              minors.fluid.logScale + b.logScale};
    }
    const auto& b = bs.solid[jump];
    return {Form::Components(YSpec::Solution(minors.solid.value, b.value)),
            minors.solid.logScale + b.logScale};
  }

 private:
  const Solver& solver_;
  const Context& context_;

  bool IsFluid(Int layer) const {
    return solver_.GetModel().GetLayer(layer).IsFluid();
  }
};

}  // namespace

template <class System>
auto RadialSolver::Impl::Solve(const System& system, Complex omega) const {
  using State = typename System::State;
  using Bs = typename System::Bs;
  constexpr auto jumps = System::Jumps;
  using Raw =
      std::array<std::array<Complex, System::Components>, System::Jumps>;
  const auto sourceCount = sources_.size();
  const auto all = All();
  // Points are ordered by radius and, at one radius on either side of a
  // discontinuity, by layer: the lower layer's is passed first going up, and
  // the upper's first going down, and one in a lower layer than a start at
  // its radius is below it.
  auto place = [](const Point& point) {
    return std::pair{point.r, point.layer};
  };
  auto byRadius = [&](const Stop& stop) { return place(stop.point); };
  const auto lowestSource = std::ranges::min(sources_, {}, place);
  const auto highestSource = std::ranges::max(sources_, {}, place);

  // The points the system reaches: the responses elsewhere are zero, as they
  // are, to the accuracy of the start, below the level at which the solution
  // from below starts.
  const auto start = system.StartLevel(lowestSource);
  auto reached = std::vector<Stop>{};
  for (const auto& stop : all) {
    if (system.Reaches(stop.point) && place(stop.point) >= place(start)) {
      reached.push_back(stop);
    }
  }
  const auto lowest = std::ranges::min(reached, {}, byRadius).point;

  auto stopsWhere = [&](auto&& keep, bool ascending) {
    auto stops = std::vector<Stop>{};
    for (const auto& stop : reached) {
      if (keep(place(stop.point))) stops.push_back(stop);
    }
    std::ranges::stable_sort(stops, {}, byRadius);
    if (!ascending) std::ranges::reverse(stops);
    return stops;
  };
  auto integrate = [&system](State& state, Int layer, double begin, double end,
                             double& step) {
    system.Integrate(state, layer, begin, end, step);
  };

  // From below, up to the highest source, kept at every source and every
  // receiver no higher.
  auto fromBelow = std::vector<State>(all.size());
  {
    const auto stops =
        stopsWhere([&](auto at) { return at <= place(highestSource); }, true);
    auto step = 0.0;
    auto [first, state] = system.FromBelow(start, lowest, step);
    WalkKeeping(
        first, highestSource, stops, state, step, omega, integrate,
        [&](std::size_t slot, const State& state) { fromBelow[slot] = state; });
  }

  // From above, down to the lowest source, kept at every source and every
  // receiver no lower.
  auto fromAbove = std::vector<State>(all.size());
  {
    const auto stops =
        stopsWhere([&](auto at) { return at >= place(lowestSource); }, false);
    auto step = 0.0;
    auto [top, state] = system.FromAbove();
    WalkKeeping(
        top, lowestSource, stops, state, step, omega, integrate,
        [&](std::size_t slot, const State& state) { fromAbove[slot] = state; });
  }

  auto responses = ResponseTable<Raw>(sourceCount, receivers_.size());
  for (auto source : std::views::iota(std::size_t{0}, sourceCount)) {
    const auto point = sources_[source];
    const auto determinant =
        System::Determinant(fromBelow[source], fromAbove[source]);

    // Carries the jumps from the source out to some receivers, all above it
    // or all below, and forms the responses there from the states kept.
    auto spread = [&](std::span<const Stop> receivers, bool above) {
      if (receivers.empty()) return;
      const auto& atSource = above ? fromBelow[source] : fromAbove[source];
      const auto& states = above ? fromAbove : fromBelow;
      const auto sign = above ? System::SignAbove : System::SignBelow;
      auto bs = System::StartBs(atSource, above);
      auto form = [&](std::size_t slot, const Bs& b) {
        auto& raw = responses(source, slot - sourceCount);
        for (auto j : std::views::iota(0, jumps)) {
          const auto y = system.Solution(states[slot], b, j, all[slot].point);
          const auto scale = ScaledComplex{sign, y.logScale} / determinant;
          for (auto i : std::views::iota(std::size_t{0}, System::Components)) {
            raw[j][i] = Rescale(y.value[i] * scale.value, scale.logScale);
          }
        }
      };
      if constexpr (System::IntegrateB) {
        auto step = 0.0;
        WalkKeeping(
            point, receivers.back().point, receivers, bs, step, omega,
            [&system](Bs& b, Int layer, double begin, double end,
                      double& step) {
              system.IntegrateBs(b, layer, begin, end, step);
            },
            form);
      } else {
        for (const auto& receiver : receivers) form(receiver.slot, bs);
      }
    };

    // Receivers at the radius of the source count as above it.
    auto receivers =
        std::vector<Stop>(reached.begin() + sourceCount, reached.end());
    std::ranges::stable_sort(receivers, {}, byRadius);
    const auto split = std::ranges::partition_point(
        receivers,
        [&](const Stop& stop) { return place(stop.point) < place(point); });
    const auto below =
        std::vector<Stop>(std::make_reverse_iterator(split), receivers.rend());
    spread(std::span{split, receivers.end()}, true);
    spread(below, false);
  }
  return responses;
}

//----------------------------------------------------------------------------
// The public interface: the responses in physical variables. The scaled ones
// are y1 = r U, y2 = zeta r V (y1 = zeta r W for the toroidal system),
// y3 = r phi and y6 = psi = r (phi' + (l+1) phi / r + 4 pi G rho U) / (4 pi G).
//----------------------------------------------------------------------------

RadialSolver::RadialSolver(const Model& model, std::vector<double> sourceRadii,
                           std::vector<Point> receivers, Gravity gravity,
                           bool attenuation, SolverOptions options)
    : impl_{std::make_shared<const Impl>(model, std::move(sourceRadii),
                                         std::move(receivers), gravity,
                                         attenuation, options)} {}

RadialSolver::RadialSolver(const Model& model, std::vector<double> sourceRadii,
                           std::vector<double> receiverRadii, Gravity gravity,
                           bool attenuation, SolverOptions options)
    : RadialSolver(
          model, std::move(sourceRadii),
          [&] {
            auto points = std::vector<Point>{};
            for (auto r : receiverRadii) {
              if (!(r > 0)) {
                throw std::invalid_argument(
                    "a receiver at the centre of the model is not "
                    "supported");
              }
              if (r > model.SurfaceRadius()) {
                throw std::invalid_argument("a receiver is outside the model");
              }
              points.push_back(Impl::SolidSide(model, r));
            }
            return points;
          }(),
          gravity, attenuation, options) {}

const Model& RadialSolver::GetModel() const { return impl_->GetModel(); }
std::size_t RadialSolver::Sources() const { return impl_->Sources().size(); }
std::size_t RadialSolver::Receivers() const {
  return impl_->Receivers().size();
}
Point RadialSolver::Source(std::size_t i) const {
  return impl_->Sources().at(i);
}
Point RadialSolver::Receiver(std::size_t i) const {
  return impl_->Receivers().at(i);
}
Gravity RadialSolver::GetGravity() const { return impl_->GetGravity(); }
bool RadialSolver::Attenuation() const { return impl_->Attenuation(); }

Point RadialSolver::StartPoint(int l, double omega, Point top,
                               Int lowestRegion) const {
  return impl_->StartPoint(l, omega, top, lowestRegion);
}

ResponseTable<ToroidalResponse> RadialSolver::Toroidal(int l,
                                                       Complex omega) const {
  const auto context = Context(GetModel(), l, omega, Attenuation());
  const auto raw = impl_->Solve(ToroidalSystem(*impl_, context), omega);
  const auto zeta = context.Zeta();
  auto responses = ResponseTable<ToroidalResponse>(Sources(), Receivers());
  for (auto receiver : std::views::iota(std::size_t{0}, Receivers())) {
    const auto point = Receiver(receiver);
    const auto r = point.r;
    const auto m = context(point.layer, r);
    const auto solid = !GetModel().GetLayer(point.layer).IsFluid();
    for (auto source : std::views::iota(std::size_t{0}, Sources())) {
      auto& response = responses(source, receiver);
      for (auto j : {0, 1}) {
        // The scaled variables y1 = zeta r W and y2 = T.
        const auto& y = raw(source, receiver)[j];
        response.W[j] = y[0] / (zeta * r);
        if (solid) {
          response.DW[j] = y[1] / (zeta * r * m.L) + response.W[j] / r;
        }
      }
    }
  }
  return responses;
}

ResponseTable<RadialResponse> RadialSolver::Radial(Complex omega) const {
  const auto context = Context(GetModel(), 0, omega, Attenuation());
  const auto raw = impl_->Solve(RadialSystem(*impl_, context), omega);
  auto responses = ResponseTable<RadialResponse>(Sources(), Receivers());
  for (auto receiver : std::views::iota(std::size_t{0}, Receivers())) {
    const auto point = Receiver(receiver);
    const auto r = point.r;
    const auto m = context(point.layer, r);
    for (auto source : std::views::iota(std::size_t{0}, Sources())) {
      auto& response = responses(source, receiver);
      for (auto j : {0, 1}) {
        // The scaled variables y1 = r U and y2 = P.
        const auto& y = raw(source, receiver)[j];
        response.U[j] = y[0] / r;
        response.DU[j] = (y[1] / r - 2.0 * m.F * response.U[j] / r) / m.C;
      }
    }
  }
  return responses;
}

ResponseTable<SpheroidalResponse> RadialSolver::Spheroidal(
    int l, Complex omega) const {
  const auto context = Context(GetModel(), l, omega, Attenuation());
  auto solve = [&]() {
    switch (GetGravity()) {
      case Gravity::None:
        return impl_->Solve(
            SpheroidalSystem<WithoutPotential<false>, Impl>(*impl_, context),
            omega);
      case Gravity::Cowling:
        return impl_->Solve(
            SpheroidalSystem<WithoutPotential<true>, Impl>(*impl_, context),
            omega);
      case Gravity::Full:
        break;
    }
    return impl_->Solve(SpheroidalSystem<WithPotential, Impl>(*impl_, context),
                        omega);
  };
  const auto raw = solve();
  const auto zeta = context.Zeta();
  const auto zeta2 = context.Zeta2();
  const auto gravitating = GetGravity() == Gravity::Full;
  auto responses = ResponseTable<SpheroidalResponse>(Sources(), Receivers());
  for (auto receiver : std::views::iota(std::size_t{0}, Receivers())) {
    const auto point = Receiver(receiver);
    const auto r = point.r;
    const auto m = context(point.layer, r);
    const auto solid = !GetModel().GetLayer(point.layer).IsFluid();
    for (auto source : std::views::iota(std::size_t{0}, Sources())) {
      auto& response = responses(source, receiver);
      for (auto j : {0, 1, 2, 3}) {
        // The scaled variables y1, y2, y3 and y6, and the tractions P and S.
        const auto& y = raw(source, receiver)[j];
        const auto U = y[0] / r;
        const auto V = y[1] / (zeta * r);
        response.U[j] = U;
        response.V[j] = V;
        response.DU[j] = (y[4] / r - m.F * (2.0 * U - zeta2 * V) / r) / m.C;
        if (solid) {
          response.DV[j] = y[5] / (zeta * r * m.L) + (V - U) / r;
        }
        if (gravitating) {
          response.Phi[j] = y[2] / r;
          response.DPhi[j] = -4 * PiG * m.rho * U -
                             (l + 1.0) * response.Phi[j] / r +
                             4 * PiG * y[3] / r;
        }
      }
    }
  }
  return responses;
}

}  // namespace YSpec

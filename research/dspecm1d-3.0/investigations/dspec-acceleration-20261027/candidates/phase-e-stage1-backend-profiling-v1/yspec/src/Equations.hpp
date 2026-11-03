#pragma once

/**
 * @file Equations.hpp
 * @brief The radial equations for one degree and one complex frequency.
 *
 * @details Each system has the form
 *
 * @f[ \frac{d}{dr}\begin{pmatrix} u \\ t \end{pmatrix} =
 *     \begin{pmatrix} T & K \\ S & -T^{\mathsf{T}} \end{pmatrix}
 *     \begin{pmatrix} u \\ t \end{pmatrix}, @f]
 *
 * with @f$K@f$ and @f$S@f$ symmetric, and is specified by the entries of the
 * blocks (Block2, Block4, Block6). The spheroidal systems are integrated as
 * minor vectors: second-order minors (five independent components) for the
 * 4x4 systems and third-order minors (fourteen, after the identities of a
 * solution regular at the centre) for the 6x6 system. The source is carried to
 * the receiver by the "b-vector" system. The orderings and the formulas are
 * those of the original YSpec (Al-Attar & Woodhouse 2008), transcribed
 * unchanged.
 *
 * The variables are scaled by the radius, as in the original code. With
 * @f$\zeta = \sqrt{l(l+1)}@f$ and @f$U, V, W@f$ the radial, consoidal and
 * toroidal displacement coefficients:
 *
 * - toroidal: @f$y_1 = \zeta r W@f$, and @f$y_2@f$ the matching traction;
 * - spheroidal, 6x6: @f$(y_1, \dots, y_6) = (rU, \zeta r V, r\phi, P, S,
 *   \psi)@f$, with @f$\phi@f$ the perturbation of the potential, @f$P@f$
 *   and @f$S@f$ the matching tractions, and the flux
 *   @f$\psi = r(\phi' + (l+1)\phi/r + 4\pi G\rho U)/(4\pi G)@f$, which
 *   is continuous and vanishes at a free surface;
 * - spheroidal, 4x4: @f$(rU, \zeta r V, P, S)@f$;
 * - radial: @f$y_1 = rU@f$, and @f$y_2@f$ the matching traction.
 *
 * Only the displacements, the potential and the flux at the receiver are
 * used.
 */

#include <YSpec/Model.hpp>
#include <YSpec/Units.hpp>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "Bessel.hpp"
#include "RungeKutta.hpp"

namespace YSpec {

/**
 * @brief The material at a radius at complex frequency: the moduli made
 * complex by the dispersion of a constant Q.
 */
struct ComplexMaterial {
  double r;       ///< The radius.
  double rho;     ///< Density.
  double g;       ///< Gravitational acceleration.
  Complex A;      ///< Love modulus A.
  Complex C;      ///< Love modulus C.
  Complex F;      ///< Love modulus F.
  Complex L;      ///< Love modulus L.
  Complex N;      ///< Love modulus N.
  Complex kappa;  ///< Voigt average bulk modulus.
  Complex mu;     ///< Voigt average shear modulus.
};

/**
 * @brief The degree and frequency of one solve, with the quantities derived
 * from them.
 *
 * @details The moduli at complex frequency @f$\omega@f$ are those of a constant
 * Q with logarithmic dispersion,
 *
 * @f[ M(\omega) = M_0\left[1 + q\,\frac{2}{\pi}
 *     \ln\frac{i\omega}{\omega_0}\right], @f]
 *
 * which is exact rather than first order in @f$q = 1/Q@f$.
 */
class Context {
 public:
  /**
   * @param model The earth model.
   * @param l The degree.
   * @param omega The complex angular frequency, non-dimensional, with a
   * negative imaginary part.
   * @param attenuation Whether to include attenuation.
   */
  Context(const Model& model, int l, Complex omega, bool attenuation = true)
      : model_{&model},
        l_{l},
        zeta2_{static_cast<double>(l) * (l + 1)},
        zeta_{std::sqrt(zeta2_)},
        omega_{omega},
        omega2_{omega * omega},
        omegaInverse2_{1.0 / omega2_},
        dispersion_{attenuation ? 2 / std::numbers::pi *
                                      std::log(Complex(0, 1) * omega /
                                               model.ReferenceFrequency())
                                : Complex{}} {}

  /** @brief The model. */
  const Model& GetModel() const { return *model_; }
  /** @brief The degree @f$l@f$. */
  int Degree() const { return l_; }
  /** @brief @f$l(l+1)@f$. */
  double Zeta2() const { return zeta2_; }
  /** @brief @f$\sqrt{l(l+1)}@f$. */
  double Zeta() const { return zeta_; }
  /** @brief The complex frequency. */
  Complex Omega() const { return omega_; }
  /** @brief The square of the frequency. */
  Complex Omega2() const { return omega2_; }
  /** @brief The inverse square of the frequency. */
  Complex OmegaInverse2() const { return omegaInverse2_; }
  /** @brief @f$(2/\pi)\ln(i\omega/\omega_0)@f$, or zero without attenuation. */
  Complex Dispersion() const { return dispersion_; }

  /** @brief The material at a radius in a layer, at this frequency. */
  ComplexMaterial operator()(Int layer, double r) const {
    return Complexify((*model_)(layer, r), r);
  }

  /** @brief A real material made complex at this frequency. */
  ComplexMaterial Complexify(const Material& m, double r) const {
    // Each modulus M becomes M (1 + q D), with q its inverse quality factor
    // and D the dispersion.
    auto at = [this](double modulus, double q) {
      return modulus * (1.0 + q * dispersion_);
    };
    return {.r = r,
            .rho = m.rho,
            .g = m.g,
            .A = at(m.A, m.qA),
            .C = at(m.C, m.qA),
            .F = at(m.F, m.qF),
            .L = at(m.L, m.qMu),
            .N = at(m.N, m.qMu),
            .kappa = at(m.kappa, m.qKappa),
            .mu = at(m.mu, m.qMu)};
  }

 private:
  const Model* model_;
  int l_;
  double zeta2_;
  double zeta_;
  Complex omega_;
  Complex omega2_;
  Complex omegaInverse2_;
  Complex dispersion_;
};

//----------------------------------------------------------------------------
// The blocks of the coefficient matrices.
//----------------------------------------------------------------------------

/** @brief A 2x2 system: @f$T = t@f$, @f$K = k@f$, @f$S = s@f$. */
struct Block2 {
  Complex t;  ///< T.
  Complex k;  ///< K.
  Complex s;  ///< S.
};

/**
 * @brief A 4x4 system, with the entries of T, K and S that are not zero.
 * Those that are real in both the solid and the fluid forms are held as
 * real, which halves the cost of the products they enter.
 */
struct Block4 {
  Complex t11;  ///< T(1,1).
  Complex t12;  ///< T(1,2).
  double t21;   ///< T(2,1).
  double t22;   ///< T(2,2).
  Complex k11;  ///< K(1,1).
  Complex k22;  ///< K(2,2).
  Complex s11;  ///< S(1,1).
  Complex s12;  ///< S(1,2) = S(2,1).
  Complex s22;  ///< S(2,2).
};

/** @brief The 6x6 self-gravitating system in a solid. */
struct Block6 {
  Complex t11;  ///< T(1,1).
  Complex t12;  ///< T(1,2).
  double t21;   ///< T(2,1).
  double t22;   ///< T(2,2).
  double t31;   ///< T(3,1).
  double t33;   ///< T(3,3).
  Complex k11;  ///< K(1,1).
  Complex k22;  ///< K(2,2).
  double k33;   ///< K(3,3).
  Complex s11;  ///< S(1,1).
  Complex s12;  ///< S(1,2) = S(2,1).
  double s13;   ///< S(1,3) = S(3,1).
  Complex s22;  ///< S(2,2).
  double s23;   ///< S(2,3) = S(3,2).
};

/** @brief The toroidal system. */
inline Block2 ToroidalBlock(const Context& c, const ComplexMaterial& m) {
  const auto r = m.r;
  return {.t = 2 / r,
          .k = 1.0 / m.L,
          .s = -c.Omega2() * m.rho + (c.Zeta2() - 2) * m.N / (r * r)};
}

/**
 * @brief The radial (l = 0) system. With gravity, the perturbation of the
 * potential has been eliminated: for l = 0 its gradient cancels the density
 * term of the gradient of gravity.
 */
inline Block2 RadialBlock(const Context& c, const ComplexMaterial& m,
                          bool gravity) {
  const auto r = m.r;
  const auto fOverC = m.F / m.C;
  // A - N - F^2/C, the modulus of the horizontal stress once the vertical
  // stress is eliminated.
  const auto lateral = m.A - m.N - m.F * fOverC;
  const auto g = gravity ? m.g : 0.0;
  return {.t = (1.0 - 2.0 * fOverC) / r,
          .k = 1.0 / m.C,
          .s = -c.Omega2() * m.rho + 4.0 * (lateral - m.rho * g * r) / (r * r)};
}

/** @brief The spheroidal system in a solid, with full self-gravitation. */
inline Block6 SolidBlock6(const Context& c, const ComplexMaterial& m) {
  const auto r = m.r;
  const auto l = static_cast<double>(c.Degree());
  const auto zeta = c.Zeta();
  const auto zeta2 = c.Zeta2();
  const auto omega2 = c.Omega2();
  const auto rho = m.rho;
  const auto g = m.g;
  const auto fOverC = m.F / m.C;
  const auto lateral = m.A - m.N - m.F * fOverC;
  return {
      .t11 = (1.0 - 2.0 * fOverC) / r,
      .t12 = zeta * fOverC / r,
      .t21 = -zeta / r,
      .t22 = 2 / r,
      .t31 = -4 * PiG * rho,
      .t33 = -l / r,
      .k11 = 1.0 / m.C,
      .k22 = 1.0 / m.L,
      .k33 = 4 * PiG,
      .s11 = -omega2 * rho + 4.0 * (lateral - rho * g * r) / (r * r),
      .s12 = zeta * (rho * g * r - 2.0 * lateral) / (r * r),
      .s13 = -rho * (l + 1) / r,
      .s22 = -omega2 * rho + (zeta2 * (lateral + m.N) - 2.0 * m.N) / (r * r),
      .s23 = rho * zeta / r};
}

/**
 * @brief The spheroidal system in a solid, without gravity or in the Cowling
 * approximation.
 */
inline Block4 SolidBlock4(const Context& c, const ComplexMaterial& m,
                          bool cowling) {
  const auto r = m.r;
  const auto zeta = c.Zeta();
  const auto zeta2 = c.Zeta2();
  const auto omega2 = c.Omega2();
  const auto rho = m.rho;
  const auto fOverC = m.F / m.C;
  const auto lateral = m.A - m.N - m.F * fOverC;
  // The weight term, rho g r, which only the Cowling approximation keeps.
  const auto weight = cowling ? rho * m.g * r : 0.0;
  return {
      .t11 = (1.0 - 2.0 * fOverC) / r,
      .t12 = zeta * fOverC / r,
      .t21 = -zeta / r,
      .t22 = 2 / r,
      .k11 = 1.0 / m.C,
      .k22 = 1.0 / m.L,
      .s11 = -omega2 * rho + 4.0 * (lateral - weight) / (r * r),
      .s12 = zeta * (weight - 2.0 * lateral) / (r * r),
      .s22 = -omega2 * rho + (zeta2 * (lateral + m.N) - 2.0 * m.N) / (r * r)};
}

/** @brief The spheroidal system in a fluid, with full self-gravitation. */
inline Block4 FluidBlock4(const Context& c, const ComplexMaterial& m) {
  const auto r = m.r;
  const auto l = static_cast<double>(c.Degree());
  const auto omega2 = c.Omega2();
  const auto rho = m.rho;
  const auto g = m.g;
  // l(l+1) / (omega^2 r): in a fluid the tangential displacement follows from
  // the pressure, and enters only through this.
  const auto tangential = c.Zeta2() * c.OmegaInverse2() / r;
  return {.t11 = (g * tangential - 1.0) / r,
          .t12 = tangential / r,
          .t21 = -4 * PiG * rho,
          .t22 = -l / r,
          .k11 = 1.0 / m.C - tangential / (r * rho),
          .k22 = 4 * PiG,
          .s11 = -omega2 * rho + rho * g * (g * tangential - 4.0) / r,
          .s12 = rho * (g * tangential - (l + 1)) / r,
          .s22 = rho * tangential / r};
}

/**
 * @brief The spheroidal system in a fluid, without gravity or in the Cowling
 * approximation.
 */
inline Block2 FluidBlock2(const Context& c, const ComplexMaterial& m,
                          bool cowling) {
  const auto r = m.r;
  const auto rho = m.rho;
  const auto g = cowling ? m.g : 0.0;
  const auto tangential = c.Zeta2() * c.OmegaInverse2() / r;
  return {.t = (g * tangential - 1.0) / r,
          .k = 1.0 / m.C - tangential / (r * rho),
          .s = -c.Omega2() * rho + rho * g * (g * tangential - 4.0) / r};
}

//----------------------------------------------------------------------------
// Derivatives.
//----------------------------------------------------------------------------

/** @brief The derivative of the solution of a 2x2 system. */
inline Vector<2> Derivative(const Block2& block, const Vector<2>& y) {
  return {block.t * y[0] + block.k * y[1], block.s * y[0] - block.t * y[1]};
}

/**
 * @brief The derivative of the second-order minors of a 4x4 system.
 * @param block The entries of the coefficient matrix.
 * @param m The minors.
 */
inline Vector<5> MinorDerivative(const Block4& block, const Vector<5>& m) {
  const auto& [t11, t12, t21, t22, k11, k22, s11, s12, s22] = block;
  // Sums and differences of the diagonal entries of T.
  const auto a11 = t11 + t22;
  const auto a22 = t11 - t22;
  auto dmdr = Vector<5>{};
  dmdr[0] = a11 * m[0] + k22 * m[2] - k11 * m[3];
  dmdr[1] = s12 * m[0] - t21 * m[2] + t12 * m[3];
  dmdr[2] = s22 * m[0] - 2.0 * t12 * m[1] + a22 * m[2] + k11 * m[4];
  dmdr[3] = -s11 * m[0] + 2.0 * t21 * m[1] - a22 * m[3] - k22 * m[4];
  dmdr[4] = -2.0 * s12 * m[1] + s11 * m[2] - s22 * m[3] - a11 * m[4];
  return dmdr;
}

/**
 * @brief The derivative of the reduced third-order minors of the 6x6 system.
 * @param block The entries of the coefficient matrix.
 * @param m The minors.
 */
inline Vector<14> MinorDerivative(const Block6& block, const Vector<14>& m) {
  const auto& [t11, t12, t21, t22, t31, t33, k11, k22, k33, s11, s12, s13, s22,
               s23] = block;
  // Signed sums of the diagonal entries of T.
  const auto a123 = t11 + t22 + t33;
  const auto a126 = t11 + t22 - t33;
  const auto a135 = t11 - t22 + t33;
  const auto a234 = -t11 + t22 + t33;
  auto dmdr = Vector<14>{};
  dmdr[0] = k33 * m[3] - k22 * m[5] + k11 * m[9] + m[0] * a123;
  dmdr[1] = -k22 * m[6] + m[0] * s13 - m[2] * t21 + m[1] * t22 - m[3] * t31;
  dmdr[2] = -k11 * m[10] + m[0] * s23 + m[2] * t11 - m[1] * t12;
  dmdr[3] = k22 * m[8] - k11 * m[11] + m[3] * a126;
  dmdr[4] = -k33 * m[7] - m[0] * s12 + m[9] * t12 - m[5] * t21 + m[2] * t31 +
            m[4] * t33;
  dmdr[5] =
      -k33 * m[8] - k11 * m[12] - m[0] * s22 - 2.0 * m[4] * t12 + m[5] * a135;
  dmdr[6] = m[2] * s12 + m[5] * s13 - m[1] * s22 - m[4] * s23 + m[10] * t12 -
            m[6] * t22 + m[8] * t31;
  dmdr[7] = m[3] * s12 - m[2] * s13 - m[1] * s23 + m[11] * t12 - m[8] * t21 -
            m[7] * t33;
  dmdr[8] = k11 * m[13] + m[3] * s22 - 2.0 * m[2] * s23 - 2.0 * m[7] * t12 -
            m[8] * a234;
  dmdr[9] = -k33 * m[11] + k22 * m[12] + m[0] * s11 + 2.0 * m[4] * t21 -
            2.0 * m[1] * t31 + m[9] * a234;
  dmdr[10] = -m[2] * s11 + m[1] * s12 - m[4] * s13 - m[9] * s23 - m[10] * t11 +
             m[6] * t21 - m[7] * t31;
  dmdr[11] = -k22 * m[13] - m[3] * s11 + 2.0 * m[1] * s13 + 2.0 * m[7] * t21 -
             m[11] * a135;
  dmdr[12] = k33 * m[13] - m[5] * s11 + 2.0 * m[4] * s12 + m[9] * s22 +
             2.0 * m[6] * t31 - m[12] * a126;
  dmdr[13] = m[8] * s11 - 2.0 * m[7] * s12 + 2.0 * m[6] * s13 - m[11] * s22 +
             2.0 * m[10] * s23 - m[13] * a123;
  return dmdr;
}

/** @brief The derivative of the b-vector of a 4x4 system. */
inline Vector<4> BDerivative(const Block4& block, const Vector<4>& b) {
  const auto& [t11, t12, t21, t22, k11, k22, s11, s12, s22] = block;
  auto dbdr = Vector<4>{};
  dbdr[0] = -s11 * b[2] - s12 * b[3] - t11 * b[0] - t21 * b[1];
  dbdr[1] = -s12 * b[2] - s22 * b[3] - t12 * b[0] - t22 * b[1];
  dbdr[2] = -k11 * b[0] + t11 * b[2] + t12 * b[3];
  dbdr[3] = -k22 * b[1] + t21 * b[2] + t22 * b[3];
  return dbdr;
}

/** @brief The derivative of the b-vector of the 6x6 system. */
inline Vector<14> BDerivative(const Block6& block, const Vector<14>& b) {
  const auto& [t11, t12, t21, t22, t31, t33, k11, k22, k33, s11, s12, s13, s22,
               s23] = block;
  auto dbdr = Vector<14>{};
  dbdr[0] = s11 * b[6] + s12 * (b[7] - b[2]) + s13 * b[8] - s22 * b[3] -
            s23 * b[4] - (t11 + t22) * b[0] + t31 * b[5];
  dbdr[1] = s11 * b[9] + s12 * b[10] - s13 * (2.0 * b[2] + b[7]) - s23 * b[3] -
            t21 * b[5] - (t11 + t33) * b[1];
  dbdr[2] = s12 * b[11] + s13 * b[12] + t12 * b[3] - t21 * b[6] - t31 * b[9];
  dbdr[3] = -k22 * b[0] - s11 * b[11] + s13 * b[13] + t21 * (b[2] - b[7]) +
            (t22 - t11) * b[3] - t31 * b[10];
  dbdr[4] = -k33 * b[1] - s11 * b[12] - s12 * b[13] - t21 * b[8] +
            t31 * (2.0 * b[2] + b[7]) + (t33 - t11) * b[4];
  dbdr[5] = s12 * b[9] - s13 * b[6] + s22 * b[10] - s23 * (b[2] + 2.0 * b[7]) -
            t12 * b[1] - (t22 + t33) * b[5];
  dbdr[6] = k11 * b[0] + s22 * b[11] + s23 * b[12] - t12 * (b[2] - b[7]) +
            (t11 - t22) * b[6];
  dbdr[7] = -s12 * b[11] + s23 * b[13] - t12 * b[3] + t21 * b[6];
  dbdr[8] = -k33 * b[5] - s12 * b[12] - s22 * b[13] - t12 * b[4] - t22 * b[8] +
            t31 * b[6] + t33 * b[8];
  dbdr[9] = k11 * b[1] + s23 * b[11] + t12 * b[10] + (t11 - t33) * b[9];
  dbdr[10] = k22 * b[5] - s13 * b[11] + t21 * b[9] + (t22 - t33) * b[10];
  dbdr[11] = -k11 * b[3] + k22 * b[6] + (t11 + t22) * b[11];
  dbdr[12] = -k11 * b[4] + k33 * b[9] + t12 * b[13] + (t11 + t33) * b[12];
  dbdr[13] = -k22 * b[8] + k33 * b[10] + t21 * b[12] - t31 * b[11] +
             (t22 + t33) * b[13];
  return dbdr;
}

//----------------------------------------------------------------------------
// Minors of starting solutions.
//----------------------------------------------------------------------------

namespace EquationDetails {

inline Complex Det2(Complex a11, Complex a12, Complex a21, Complex a22) {
  return a11 * a22 - a12 * a21;
}

// The minor of rows i and j of the 4x2 matrix whose columns are y1 and y2.
inline Complex Minor2(const Vector<4>& y1, const Vector<4>& y2, int i, int j) {
  return Det2(y1[i], y2[i], y1[j], y2[j]);
}

// The minor of rows i, j and k of the 6x3 matrix with columns y1, y2, y3.
inline Complex Minor3(const Vector<6>& y1, const Vector<6>& y2,
                      const Vector<6>& y3, int i, int j, int k) {
  return y1[i] * Det2(y2[j], y3[j], y2[k], y3[k]) -
         y2[i] * Det2(y1[j], y3[j], y1[k], y3[k]) +
         y3[i] * Det2(y1[j], y2[j], y1[k], y2[k]);
}

}  // namespace EquationDetails

/** @brief The five independent second-order minors of two 4-vectors. */
inline Vector<5> Minors(const Vector<4>& y1, const Vector<4>& y2) {
  using EquationDetails::Minor2;
  return {Minor2(y1, y2, 0, 1), Minor2(y1, y2, 0, 2), Minor2(y1, y2, 0, 3),
          Minor2(y1, y2, 1, 2), Minor2(y1, y2, 2, 3)};
}

/** @brief The fourteen independent third-order minors of three 6-vectors. */
inline Vector<14> Minors(const Vector<6>& y1, const Vector<6>& y2,
                         const Vector<6>& y3) {
  auto m = [&](int i, int j, int k) {
    return EquationDetails::Minor3(y1, y2, y3, i - 1, j - 1, k - 1);
  };
  return {m(1, 2, 3), m(1, 2, 4), m(1, 2, 5), m(1, 2, 6), m(1, 3, 4),
          m(1, 3, 5), m(1, 4, 5), m(1, 4, 6), m(1, 5, 6), m(2, 3, 4),
          m(2, 4, 5), m(2, 4, 6), m(3, 4, 5), m(4, 5, 6)};
}

//----------------------------------------------------------------------------
// Starting solutions: those regular at the centre of a homogeneous sphere.
//----------------------------------------------------------------------------

// In each, ratio is the ratio of spherical Bessel functions
// x j_{l+1}(x) / j_l(x) at x = k r, for the wavenumber k of the solution.

/** @brief The toroidal starting solution. */
inline Vector<2> ToroidalStart(const Context& c, const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto r = m.r;
  const auto kr2 = c.Omega2() * m.rho * r * r / m.mu;
  const auto ratio = BesselRatio(c.Degree(), kr2);
  return {c.Zeta() * r, c.Zeta() * m.mu * ((l - 1) - ratio)};
}

/** @brief The radial starting solution, with or without gravity. */
inline Vector<2> RadialStart(const Context& c, const ComplexMaterial& m,
                             bool gravity) {
  const auto r = m.r;
  const auto pModulus = m.kappa + 4.0 * m.mu / 3.0;
  // Gravity adds 16 pi G rho / 3 to omega^2.
  const auto shift = gravity ? 16 * PiG * m.rho / 3 : 0.0;
  const auto k2 = m.rho * (c.Omega2() + shift) / pModulus;
  const auto ratio = BesselRatio(0, k2 * r * r);
  return {ratio, r * pModulus * k2 - 4.0 * m.mu * ratio / r};
}

/**
 * @brief The spheroidal starting solutions in a solid, with
 * full self-gravitation.
 *
 * @details The three solutions of a homogeneous self-gravitating sphere
 * regular at its centre (Takeuchi & Saito 1972): two with wavenumbers that
 * solve a quadratic in @f$k^2@f$, gravity coupling the P and S waves, and one
 * without spatial oscillation. @c f and @c h are theirs.
 */
inline std::array<Vector<6>, 3> SolidStartSolutions6(const Context& c,
                                                     const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto zeta = c.Zeta();
  const auto zeta2 = c.Zeta2();
  const auto omega2 = c.Omega2();
  const auto r = m.r;
  const auto rho = m.rho;
  const auto alpha2 = (m.kappa + 4.0 * m.mu / 3.0) / rho;
  const auto beta2 = m.mu / rho;
  const auto kS2 = omega2 / beta2;
  const auto kP2 = (omega2 + 16 * PiG * rho / 3) / alpha2;
  const auto coupling =
      std::pow(8 * PiG * zeta * rho / 3, 2) / (alpha2 * beta2);

  auto y = std::array<Vector<6>, 3>{};
  for (auto i : {0, 1}) {
    const auto sign = i == 0 ? 1.0 : -1.0;
    const auto k2 =
        0.5 * (kS2 + kP2) +
        0.5 * sign * std::sqrt((kS2 - kP2) * (kS2 - kP2) + coupling);
    const auto f = 0.75 * beta2 * (k2 - kS2) / (PiG * rho);
    const auto h = f - (l + 1);
    const auto ratio = BesselRatio(c.Degree(), k2 * r * r);
    y[i][0] = l * h - f * ratio;
    y[i][1] = zeta * (h + ratio);
    y[i][2] = -4 * PiG * rho * r * f;
    y[i][3] = -rho * alpha2 * f * r * k2 + 2.0 * m.mu * l * (l - 1) * h / r +
              2.0 * m.mu * (2.0 * f + zeta2) * ratio / r;
    y[i][4] = m.mu * zeta *
              (r * k2 + 2.0 * (l - 1) * h / r - 2.0 * (f + 1.0) * ratio / r);
    y[i][5] = -rho * (zeta2 + (l + 1) * f);
  }
  y[2][0] = l;
  y[2][1] = zeta;
  y[2][2] = (omega2 - 4 * PiG * rho * l / 3) * r;
  y[2][3] = 2.0 * l * (l - 1) * m.mu / r;
  y[2][4] = 2.0 * zeta * (l - 1) * m.mu / r;
  y[2][5] =
      0.25 * ((2 * l + 1) * omega2 - 8 * PiG * l * (l - 1) * rho / 3) / PiG;
  return y;
}

/** @brief The minors of the solutions of SolidStartSolutions6. */
inline Vector<14> SolidStart6(const Context& c, const ComplexMaterial& m) {
  const auto y = SolidStartSolutions6(c, m);
  return Minors(y[0], y[1], y[2]);
}

/**
 * @brief The spheroidal starting solutions in a solid, without
 * gravity (used in the Cowling approximation too): a P and an S solution
 * (Takeuchi & Saito 1972).
 */
inline std::array<Vector<4>, 2> SolidStartSolutions4(const Context& c,
                                                     const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto zeta = c.Zeta();
  const auto omega2 = c.Omega2();
  const auto r = m.r;
  const auto rho = m.rho;
  const auto alpha2 = (m.kappa + 4.0 * m.mu / 3.0) / rho;
  const auto beta2 = m.mu / rho;
  auto p = Vector<4>{};
  {
    const auto kr2 = omega2 * r * r / alpha2;
    const auto ratio = BesselRatio(c.Degree(), kr2);
    p[0] = l - ratio;
    p[1] = zeta;
    p[2] = (-rho * alpha2 * kr2 + 2.0 * m.mu * (l * (l - 1) + 2.0 * ratio)) / r;
    p[3] = 2.0 * zeta * m.mu * ((l - 1) - ratio) / r;
  }
  auto s = Vector<4>{};
  {
    const auto kr2 = omega2 * r * r / beta2;
    const auto ratio = BesselRatio(c.Degree(), kr2);
    s[0] = -l * (l + 1);
    s[1] = zeta * (ratio - (l + 1));
    s[2] = -2.0 * m.mu * (l * (l * l - 1) - l * (l + 1) * ratio) / r;
    s[3] = zeta * m.mu * (kr2 - 2.0 * (l * l - 1) - 2.0 * ratio) / r;
  }
  return {p, s};
}

/** @brief The minors of the solutions of SolidStartSolutions4. */
inline Vector<5> SolidStart4(const Context& c, const ComplexMaterial& m) {
  const auto [p, s] = SolidStartSolutions4(c, m);
  return Minors(p, s);
}

/**
 * @brief The spheroidal starting solutions in a fluid, with
 * full self-gravitation (Takeuchi & Saito 1972, whose @c gamma, @c f and
 * @c h these are).
 */
inline std::array<Vector<4>, 2> FluidStartSolutions4(const Context& c,
                                                     const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto zeta2 = c.Zeta2();
  const auto omega2 = c.Omega2();
  const auto r = m.r;
  const auto rho = m.rho;
  const auto alpha2 = m.kappa / rho;
  const auto gamma = 4 * PiG * rho / 3;
  const auto k2 =
      (omega2 + 4 * gamma - zeta2 * gamma * gamma * c.OmegaInverse2()) / alpha2;
  const auto f = -omega2 / gamma;
  const auto h = f - (l + 1);
  const auto ratio = BesselRatio(c.Degree(), k2 * r * r);
  const auto oscillating =
      Vector<4>{l * h - f * ratio, -4 * PiG * rho * r * f,
                -rho * alpha2 * f * r * k2, -rho * (zeta2 + (l + 1) * f)};
  const auto steady = Vector<4>{
      l, (omega2 - 4 * PiG * rho * l / 3) * r, 0.0,
      0.25 * ((2 * l + 1) * omega2 - 8 * PiG * l * (l - 1) * rho / 3) / PiG};
  return {oscillating, steady};
}

/** @brief The minors of the solutions of FluidStartSolutions4. */
inline Vector<5> FluidStart4(const Context& c, const ComplexMaterial& m) {
  const auto [oscillating, steady] = FluidStartSolutions4(c, m);
  return Minors(oscillating, steady);
}

/**
 * @brief The spheroidal starting solution in a fluid, without gravity (used
 * in the Cowling approximation too).
 */
inline Vector<2> FluidStart2(const Context& c, const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto r = m.r;
  const auto alpha2 = m.kappa / m.rho;
  const auto kr2 = c.Omega2() * r * r / alpha2;
  const auto ratio = BesselRatio(c.Degree(), kr2);
  return {l - ratio, -m.rho * alpha2 * kr2 / r};
}

//----------------------------------------------------------------------------
// Solid-fluid boundaries, the secular determinant, the b-vector and the
// solution.
//----------------------------------------------------------------------------

/** @brief Third-order minors in a solid to second-order minors in a fluid. */
inline Vector<5> SolidToFluid(const Vector<14>& solid) {
  return {-solid[5], -solid[6], solid[8], -solid[12], solid[13]};
}

/** @brief Second-order minors in a solid to the fluid solution (no gravity or
 * Cowling). */
inline Vector<2> SolidToFluid(const Vector<5>& solid) {
  return {-solid[2], -solid[4]};
}

/** @brief Second-order minors in a fluid to third-order minors in a solid. */
inline Vector<14> FluidToSolid(const Vector<5>& fluid) {
  auto solid = Vector<14>{};
  solid[0] = fluid[0];
  solid[1] = fluid[1];
  solid[3] = fluid[2];
  solid[9] = -fluid[3];
  solid[11] = -fluid[4];
  return solid;
}

/** @brief The fluid solution to second-order minors in a solid (no gravity or
 * Cowling). */
inline Vector<5> FluidToSolid(const Vector<2>& fluid) {
  auto solid = Vector<5>{};
  solid[0] = fluid[0];
  solid[3] = -fluid[1];
  return solid;
}

/** @brief The secular determinant from minors integrated from below and from
 * above (4x4 systems). */
inline Complex SecularDeterminant(const Vector<5>& below,
                                  const Vector<5>& above) {
  return below[4] * above[0] + 2.0 * below[1] * above[1] + below[3] * above[2] +
         below[2] * above[3] + below[0] * above[4];
}

/** @brief The secular determinant from minors integrated from below and from
 * above (6x6 system). */
inline Complex SecularDeterminant(const Vector<14>& below,
                                  const Vector<14>& above) {
  return -below[13] * above[0] + 2.0 * below[6] * above[1] +
         2.0 * below[10] * above[2] + below[12] * above[3] +
         2.0 * below[7] * above[4] + below[11] * above[5] -
         2.0 * below[1] * above[6] - 2.0 * below[4] * above[7] -
         below[9] * above[8] + below[8] * above[9] -
         2.0 * below[2] * above[10] - below[5] * above[11] -
         below[3] * above[12] + below[0] * above[13];
}

/**
 * @brief The starting b-vector at the source from the minors and a jump (4x4
 * systems).
 * @param m The minors at the source, from below or from above.
 * @param s The jump in the solution at the source.
 */
inline Vector<4> BStart(const Vector<5>& m, const Vector<4>& s) {
  return {-m[4] * s[1] - m[1] * s[2] - m[3] * s[3],
          m[4] * s[0] - m[2] * s[2] + m[1] * s[3],
          m[1] * s[0] + m[2] * s[1] - m[0] * s[3],
          m[3] * s[0] - m[1] * s[1] + m[0] * s[2]};
}

/**
 * @brief The starting b-vector at the source from the minors and a jump (6x6
 * system).
 * @param m The minors at the source, from below or from above.
 * @param s The jump in the solution at the source.
 */
inline Vector<14> BStart(const Vector<14>& m, const Vector<6>& s) {
  auto b = Vector<14>{};
  b[0] = m[13] * s[2] - m[6] * s[3] - m[10] * s[4] - m[12] * s[5];
  b[1] = -m[13] * s[1] - m[7] * s[3] - m[11] * s[4] + m[10] * s[5];
  b[2] = m[6] * s[1] + m[7] * s[2] + m[1] * s[4] + m[4] * s[5];
  b[3] = m[10] * s[1] + m[11] * s[2] - m[1] * s[3] + m[9] * s[5];
  b[4] = m[12] * s[1] - m[10] * s[2] - m[4] * s[3] - m[9] * s[4];
  b[5] = m[13] * s[0] - m[8] * s[3] + m[7] * s[4] - m[6] * s[5];
  b[6] = -m[6] * s[0] + m[8] * s[2] + m[2] * s[4] + m[5] * s[5];
  b[7] = -m[10] * s[0] - m[7] * s[2] - m[2] * s[3] - m[4] * s[5];
  b[8] = -m[12] * s[0] + m[6] * s[2] - m[5] * s[3] + m[4] * s[4];
  b[9] = -m[7] * s[0] - m[8] * s[1] + m[3] * s[4] - m[2] * s[5];
  b[10] = -m[11] * s[0] + m[7] * s[1] - m[3] * s[3] + m[1] * s[5];
  b[11] = m[1] * s[0] + m[2] * s[1] + m[3] * s[2] - m[0] * s[5];
  b[12] = m[4] * s[0] + m[5] * s[1] - m[2] * s[2] + m[0] * s[4];
  b[13] = m[9] * s[0] - m[4] * s[1] + m[1] * s[2] - m[0] * s[3];
  return b;
}

/**
 * @brief The solution from the minors and the b-vector (4x4 systems).
 * @param m The minors, integrated from the side away from the source.
 * @param b The b-vector, integrated from the source.
 */
inline Vector<4> Solution(const Vector<5>& m, const Vector<4>& b) {
  return {m[0] * b[1] + m[1] * b[2] + m[2] * b[3],
          -m[0] * b[0] + m[3] * b[2] - m[1] * b[3],
          -m[1] * b[0] - m[3] * b[1] + m[4] * b[3],
          -m[2] * b[0] + m[1] * b[1] - m[4] * b[2]};
}

/**
 * @brief The solution from the minors and the b-vector (6x6 system).
 * @param m The minors, integrated from the side away from the source.
 * @param b The b-vector, integrated from the source.
 */
inline Vector<6> Solution(const Vector<14>& m, const Vector<14>& b) {
  auto y = Vector<6>{};
  y[0] = m[2] * b[2] + m[0] * b[5] + m[1] * b[6] + 2.0 * m[2] * b[7] +
         m[3] * b[8] + m[4] * b[9] + m[5] * b[10] + m[6] * b[11] +
         m[7] * b[12] + m[8] * b[13];
  y[1] = -m[0] * b[1] - m[1] * (2.0 * b[2] + b[7]) - m[2] * b[3] - m[3] * b[4] -
         m[4] * b[10] - m[7] * b[13] + m[9] * b[9] + m[10] * b[11] +
         m[11] * b[12];
  y[2] = m[0] * b[0] - m[1] * b[8] + m[2] * b[4] + m[4] * (b[7] - b[2]) -
         m[5] * b[3] + m[6] * b[13] - m[9] * b[6] - m[10] * b[12] +
         m[12] * b[11];
  y[3] = m[1] * b[0] + m[4] * b[1] - m[6] * b[3] - m[7] * b[4] + m[9] * b[5] -
         m[10] * (b[2] + 2.0 * b[7]) - m[11] * b[8] - m[12] * b[10] +
         m[13] * b[13];
  y[4] = m[2] * b[0] - m[4] * b[5] + m[5] * b[1] + m[6] * (2.0 * b[2] + b[7]) +
         m[7] * b[8] - m[8] * b[4] + m[10] * b[6] + m[12] * b[9] -
         m[13] * b[12];
  y[5] = m[1] * b[5] - m[2] * b[1] + m[3] * b[0] + m[6] * b[10] +
         m[7] * (b[2] - b[7]) + m[8] * b[3] - m[10] * b[9] + m[11] * b[6] +
         m[13] * b[11];
  return y;
}

//----------------------------------------------------------------------------
// The b-vector across a boundary between solid and fluid, and the consoidal
// displacement in a fluid.
//
// The b-vector is the solution of the adjoint system that turns the minors
// into the response. For the 4x4 systems it is a covector c, and with y1 and
// y2 the solutions whose minors are integrated the response is
// y = y1 (c.y2) - y2 (c.y1). For the 6x6 system it is a 2-form B, and with y1,
// y2, y3 the solutions the response is the sum of y1 B(y2, y3) and its cyclic
// permutations. Of its fifteen components the fourteen kept are, as the
// formula of Solution shows, B01 = b0, B02 = b1, B04 = b3, B05 = b4,
// B12 = b5, B13 = b6, B15 = b8, B23 = b9, B24 = b10, B34 = b11, B35 = b12 and
// B45 = b13, while B03, B14 and B25 enter only through B14 - B25 = b2 + 2 b7
// and B03 - B25 = 2 b2 + b7 (indices from zero, in the order of the 6x6
// variables).
//
// At a boundary, FluidToSolid gives the minors of the fluid solutions, with
// V = 0, together with the solution e_V that slips; SolidToFluid gives those
// of the combinations of the solid solutions with S = 0. The response keeps
// its continuous components and has S = 0 on the solid side, and so:
//
// - into a fluid, c_k = B(e_V, e_k) for each fluid variable k, and for the
//   4x4 system the fluid's constant is c(e_V);
// - into a solid, B(x, e_S) = -c.x for every x, and for the 4x4 system
//   c = -b e_S, with b the fluid's constant.
//----------------------------------------------------------------------------

/** @brief The b-vector of the 4x4 system into the constant of the fluid. */
inline Vector<1> BIntoFluid(const Vector<4>& b) { return {b[1]}; }

/** @brief The constant of the fluid into the b-vector of the 4x4 system. */
inline Vector<4> BIntoSolid(const Vector<1>& b) { return {0, 0, 0, -b[0]}; }

/** @brief The b-vector of the 6x6 system into that of the fluid's 4x4. */
inline Vector<4> BIntoFluid(const Vector<14>& b) {
  return {-b[0], b[5], b[6], b[8]};
}

/** @brief The b-vector of the fluid's 4x4 system into that of the 6x6. */
inline Vector<14> BIntoSolid(const Vector<4>& b) {
  auto solid = Vector<14>{};
  solid[3] = -b[0];
  solid[10] = -b[1];
  solid[11] = -b[2];
  solid[13] = b[3];
  return solid;
}

/**
 * @brief The consoidal displacement in a fluid, @f$\zeta r V@f$, which is not
 * a variable of the fluid's systems.
 *
 * @details In a fluid the tangential traction vanishes, and its equation,
 * @f$0 = \zeta\rho g\,y_1/r - \omega^2\rho\,y_2 + \zeta\rho\,y_3/r - \zeta
 * y_4/r@f$ in the variables of the 6x6 system, gives @f$y_2@f$.
 * @param c The degree and frequency.
 * @param m The material.
 * @param rU The radial displacement, @f$y_1@f$.
 * @param rPhi The potential, @f$y_3@f$, zero without self-gravitation.
 * @param P The radial traction, @f$y_4@f$.
 * @param gravity Whether the equilibrium gravity enters.
 */
inline Complex FluidConsoidal(const Context& c, const ComplexMaterial& m,
                              Complex rU, Complex rPhi, Complex P,
                              bool gravity) {
  const auto g = gravity ? m.g : 0.0;
  return c.Zeta() * (m.rho * (g * rU + rPhi) - P) / (c.Omega2() * m.rho * m.r);
}

}  // namespace YSpec

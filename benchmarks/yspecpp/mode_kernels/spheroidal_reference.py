"""Independent homogeneous-sphere full-gravity Ritz reference for one spheroidal mode."""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.linalg import eigh, solve

R, RHO, VS, VP, G = 6_371_000.0, 5515.0, 4500.0, 8000.0, 6.6723e-11
ELL = 2
LAM, K = ELL * (ELL + 1), np.sqrt(ELL * (ELL + 1))
GAMMA = 4 * np.pi * G * RHO * R**2 / VS**2
A0, C0, F0, L0, N0 = (VP / VS)**2, (VP / VS)**2, (VP / VS)**2 - 2, 1.0, 1.0


def _basis(x, order):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    z = 2 * x * x - 1
    p = np.column_stack([Legendre.basis(j)(z) for j in range(order)])
    dpz = np.column_stack([Legendre.basis(j).deriv()(z) for j in range(order)])
    dp = dpz * (4 * x[:, None])
    high = x[:, None] ** (ELL + 1) * p
    dhigh = (ELL + 1) * x[:, None] ** ELL * p + x[:, None] ** (ELL + 1) * dp
    pot = x[:, None] ** ELL * p
    dpot = ELL * x[:, None] ** (ELL - 1) * p + x[:, None] ** ELL * dp
    zero = np.zeros_like(high)
    U = np.column_stack((ELL * x ** (ELL - 1), high, zero))
    V = np.column_stack((K * x ** (ELL - 1), zero, high))
    dU = np.column_stack((ELL * (ELL - 1) * x ** (ELL - 2), dhigh, zero))
    dV = np.column_stack((K * (ELL - 1) * x ** (ELL - 2), zero, dhigh))
    return U, V, dU, dV, pot, dpot


def _gram(a, b, w):
    return a.T @ (w[:, None] * b)


def _density_tail(x, ucoef, vcoef, order):
    # The integrand is a regular polynomial; fixed Gauss integration is independent of either solver.
    nodes, weights = leggauss(80)
    x = np.atleast_1d(np.asarray(x, dtype=float))
    result = np.empty_like(x)
    for first in range(0, x.size, 256):
        last = min(first + 256, x.size); xx = x[first:last]
        q = xx[:, None] + (1 - xx[:, None]) * (nodes[None, :] + 1) / 2
        U, V, *_ = _basis(q.ravel(), order)
        uq = (U @ ucoef).reshape(q.shape); vq = (V @ vcoef).reshape(q.shape)
        result[first:last] = (1 - xx) / 2 * np.sum(weights[None, :] * uq * (2 * uq - K * vq) / q, axis=1)
    return result


def solve_reference(order=10, curve_x=None):
    nodes, weights = leggauss(180)
    x, w = (nodes + 1) / 2, weights / 2
    U, V, dU, dV, P, dP = _basis(x, order)
    t, shear = 2 * U - K * V, x[:, None] * dV - V + K * U
    mass = _gram(U, U, w * x * x) + _gram(V, V, w * x * x)
    elastic = (_gram(x[:, None] * dU, x[:, None] * dU, w * C0) +
               _gram(t, t, w * (A0 - N0)) + _gram(shear, shear, w * L0) +
               _gram(V, V, w * (LAM - 2) * N0))
    cross = _gram(x[:, None] * dU, t, w * F0)
    elastic += cross + cross.T
    gravity = _gram(U, U, w * GAMMA * x * x)
    cross_g = _gram(U, t, -w * (GAMMA * x / 3) * x)
    gravity += cross_g + cross_g.T
    coupling = _gram(U, dP, w * x * x) + _gram(V, P, w * K * x)
    p_surface = _basis(np.array([1.0]), order)[4][0]
    potential = (_gram(dP, dP, w * x * x) + _gram(P, P, w * LAM) +
                 (ELL + 1) * np.outer(p_surface, p_surface)) / GAMMA
    stiffness = elastic + gravity - coupling @ solve(potential, coupling.T, assume_a="pos")
    values, vectors = eigh(stiffness, mass)
    positive = np.flatnonzero(values > 0)
    if not positive.size:
        raise RuntimeError("full-gravity Ritz pencil has no positive spheroidal eigenvalue")
    mode = vectors[:, positive[0]]
    potential_coef = -solve(potential, coupling.T @ mode, assume_a="pos")
    uc, vc, pc = mode, mode, potential_coef
    omega_nd = float(np.sqrt(values[positive[0]]))
    I = float(mode @ mass @ mode)

    def profiles(at):
        B = _basis(at, order)
        u, v, du, dv, p, dp = (B[0] @ uc, B[1] @ vc, B[2] @ uc, B[3] @ vc,
                               B[4] @ pc, B[5] @ pc)
        tt, ss = 2 * u - K * v, at * dv - v + K * u
        kernels = {
            "A": tt**2, "C": (at * du)**2, "F": 2 * at * du * tt,
            "L": ss**2, "N": (LAM - 2) * v**2 - tt**2,
        }
        g = GAMMA * at / 3
        tail = _density_tail(at, uc, vc, order)
        kernels["rho"] = (2 * GAMMA * at**2 * u**2 - 2 * g * at * u * tt +
            2 * (at**2 * u * dp + K * at * v * p) - omega_nd**2 * at**2 * (u**2 + v**2) -
            2 * GAMMA * at**2 * tail)
        factors = {"A": A0, "C": C0, "F": F0, "L": L0, "N": N0, "rho": 1.0}
        # Fractional-parameter derivative integrands per x, in physical rad/s.
        return {name: VS / R * factors[name] * value / (2 * omega_nd * I)
                for name, value in kernels.items()}

    at = np.array([1.0])
    Bu, Bv, Bdu, Bdv, Bp, Bdp = _basis(at, order)
    us, vs, dus, dvs, ps, dps = (float((Bu @ uc)[0]), float((Bv @ vc)[0]), float((Bdu @ uc)[0]),
                                  float((Bdv @ vc)[0]), float((Bp @ pc)[0]), float((Bdp @ pc)[0]))
    traction_r = C0 * dus + F0 * (2 * us - K * vs)
    traction_s = dvs - vs + K * us
    potential_flux = dps + GAMMA * us + (ELL + 1) * ps
    scale = max(abs(C0 * dus), abs(F0 * (2 * us - K * vs)), abs(dvs - vs + K * us), 1e-30)
    boundary = max(abs(traction_r), abs(traction_s), abs(potential_flux)) / scale
    if curve_x is None:
        curve_x = x
    curves = profiles(np.asarray(curve_x, dtype=float))
    quadrature_profiles = profiles(x)
    integral = {name: float(np.dot(w, quadrature_profiles[name])) for name in curves}
    return {"order": order, "omega_nd": omega_nd, "omega": omega_nd * VS / R,
            "frequency_mHz": omega_nd * VS / R / (2 * np.pi) * 1e3,
            "mass_norm": I, "boundary_residual": boundary,
            "traction_r_residual": float(traction_r), "traction_s_residual": float(traction_s),
            "potential_flux_residual": float(potential_flux), "curves": curves,
            "derivatives": integral, "mode_at_quadrature": (U @ uc, V @ vc),
            "quadrature_x": x, "quadrature_w": w}


def compare_orders(low, high):
    x, w = high["quadrature_x"], high["quadrature_w"]
    U8, V8 = low["mode_at_quadrature"]
    U10, V10 = high["mode_at_quadrature"]
    # Reevaluate the lower-order mode on the higher-order quadrature is done by callers via its stored sampled fields.
    # Both references use the same quadrature, so this is the mass-weighted signed overlap.
    a = np.r_[U8, V8]; b = np.r_[U10, V10]
    mass_inner = np.sum(w * x * x * (U8 * U10 + V8 * V10))
    norm_a = np.sqrt(np.sum(w * x * x * (U8**2 + V8**2)))
    norm_b = np.sqrt(np.sum(w * x * x * (U10**2 + V10**2)))
    overlap = abs(mass_inner) / (norm_a * norm_b)
    return {"frequency_relative_change": abs(low["omega"] - high["omega"]) / high["omega"],
            "displacement_overlap": float(overlap)}

if __name__ == "__main__":
    for n in (8, 10):
        r = solve_reference(n)
        print(n, r["frequency_mHz"], r["boundary_residual"], r["derivatives"])

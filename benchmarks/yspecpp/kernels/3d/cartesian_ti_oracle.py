#!/usr/bin/env python3
"""Independent Cartesian rank-four TI contractions and sparse SEM assembly.

No DSpec spin, field, source-treatment, or kernel helpers are imported here.
The module consumes plain resolved mesh/mode arrays exported by the C++ bridge.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import argparse
import inspect
import json
from math import factorial, pi
import os
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import scipy
from scipy.io import mmread
from scipy.sparse import block_diag, coo_matrix, csr_matrix
from scipy.sparse.linalg import LinearOperator, gmres, splu

BASELINE_SPECTRUM_REL_TOL = 1.0e-9


def _gmres_tolerance_keyword() -> str:
    parameters = inspect.signature(gmres).parameters
    if "rtol" in parameters:
        return "rtol"
    if "tol" in parameters:
        return "tol"
    raise RuntimeError("installed SciPy GMRES exposes neither rtol nor tol")


@dataclass(frozen=True)
class Mode:
    degree: int
    order: int
    kind: str  # "U", "V", or "W"
    family: str = "S"  # "R", "S", or "T" solver block


@dataclass(frozen=True)
class Element:
    index: int
    nodes_model: np.ndarray
    dr_weights_model: np.ndarray
    derivative_model: np.ndarray  # [quadrature node, local basis node]
    dofs: tuple[np.ndarray, ...]  # one local-to-global map per Mode
    field_dofs: dict[tuple[str, int, str], np.ndarray]


@dataclass(frozen=True)
class ModalOracle:
    modes: tuple[Mode, ...]
    elements: tuple[Element, ...]
    offsets: dict[tuple[str, int, int], tuple[int, int]]
    matrices: dict[tuple[int, str, int], csr_matrix]
    source_local_matrices: dict[tuple[int, str, int], csr_matrix]
    source_element: int
    rhs: dict[int, np.ndarray]
    lift: dict[int, np.ndarray]
    receiver: dict[int, np.ndarray]
    metadata: dict[str, object]
    size: int


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _family_orders(family: str, degree: int) -> range | tuple[int, ...]:
    return (0,) if family == "R" else range(-degree, degree + 1)


def _barycentric_derivative(nodes: np.ndarray) -> np.ndarray:
    nodes = np.asarray(nodes, dtype=float)
    if nodes.ndim != 1 or len(nodes) < 2 or len(np.unique(nodes)) != len(nodes):
        raise ValueError("SEM element nodes must be a distinct one-dimensional array")
    weights = np.ones(len(nodes))
    for i in range(len(nodes)):
        weights[i] = 1.0/np.prod(nodes[i]-np.delete(nodes, i))
    derivative = np.zeros((len(nodes), len(nodes)))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i != j:
                derivative[i, j] = weights[j]/(weights[i]*(nodes[i]-nodes[j]))
        derivative[i, i] = -np.sum(derivative[i])
    return derivative


def _factor_h0(matrix: csr_matrix):
    """Factor the same symmetric row-equilibrated H0 used by DSpecM1D."""
    row_max = np.asarray(abs(matrix).max(axis=1).toarray()).ravel()
    if np.any(~np.isfinite(row_max)) or np.any(row_max <= 0.0):
        raise ValueError("H0 must have finite, nonzero rows for equilibration")
    scale = 1.0/np.sqrt(row_max)
    equilibrated = matrix.multiply(scale[:, None]).multiply(scale[None, :]).tocsc()
    factor = splu(equilibrated, permc_spec="COLAMD", options={"Equil": False})
    return factor, scale


def _solve_h0_factor(entry, rhs: np.ndarray,
                     transpose: str | None = None) -> np.ndarray:
    factor, scale = entry
    if rhs.ndim == 1:
        return scale*factor.solve(scale*rhs, trans=transpose or "N")
    return scale[:, None]*factor.solve(scale[:, None]*rhs,
                                       trans=transpose or "N")


def load_modal_oracle(oracle_dir: str | Path) -> ModalOracle:
    """Load sparse H0, source, receiver, and SEM maps exported by DSpecM1D."""
    root = Path(oracle_dir)
    metadata = json.loads((root/"meta.json").read_text(encoding="utf-8"))
    scale = float(metadata["length_scale_m"])
    matrix_paths = sorted(root.glob("H0_k*_*.mtx"))
    matrices: dict[tuple[int, str, int], csr_matrix] = {}
    for path in matrix_paths:
        match = re.fullmatch(r"H0_k(\d+)_([RST])_l(\d+)\.mtx", path.name)
        if not match:
            raise ValueError(f"unrecognized modal matrix filename {path.name}")
        key = (int(match.group(1)), match.group(2), int(match.group(3)))
        matrices[key] = csr_matrix(mmread(path), dtype=np.complex128)
    if not matrices:
        raise ValueError("no exported DSpec H0 matrices were found")
    local_matrices: dict[tuple[int, str, int], csr_matrix] = {}
    for path in sorted(root.glob("source_local_H0_k*_*.mtx")):
        match = re.fullmatch(r"source_local_H0_k(\d+)_([RST])_l(\d+)\.mtx", path.name)
        if not match:
            raise ValueError(f"unrecognized source-local matrix filename {path.name}")
        key = (int(match.group(1)), match.group(2), int(match.group(3)))
        local_matrices[key] = csr_matrix(mmread(path), dtype=np.complex128)
    frequency_bins = sorted({key[0] for key in matrices})
    if frequency_bins != [1, 2, 3]:
        raise ValueError(f"expected H0 blocks for bins 1,2,3, found {frequency_bins}")

    template = {key[1:]: matrix for key, matrix in matrices.items() if key[0] == frequency_bins[0]}
    for bin_index in frequency_bins[1:]:
        for family_degree, reference in template.items():
            current = matrices.get((bin_index, *family_degree))
            if current is None or current.shape != reference.shape:
                raise ValueError("modal H0 block dimensions change across active bins")

    block_keys = []
    for family, degree in sorted(template, key=lambda item: (item[1], item[0])):
        for order in _family_orders(family, degree):
            block_keys.append((family, degree, order))
    offsets: dict[tuple[str, int, int], tuple[int, int]] = {}
    cursor = 0
    for family, degree, order in block_keys:
        size = template[(family, degree)].shape[0]
        offsets[(family, degree, order)] = (cursor, size)
        cursor += size

    modes = []
    for family, degree, order in block_keys:
        kinds = ("U",) if family == "R" else (("U", "V") if family == "S" else ("W",))
        modes.extend(Mode(degree, order, kind, family) for kind in kinds)
    dof_rows = _csv_rows(root/"mesh_dofs.csv")
    dof_lookup = {}
    for row in dof_rows:
        key = (row["family"], int(row["degree"]), int(row["element"]),
               int(row["node"]), row["field"])
        dof_lookup[key] = int(row["index"])
    element_rows = _csv_rows(root.parent/"mesh_nodes.csv")
    grouped_elements: dict[int, list[dict[str, str]]] = {}
    for row in element_rows:
        grouped_elements.setdefault(int(row["element"]), []).append(row)
    elements = []
    for element_index, rows in sorted(grouped_elements.items()):
        rows.sort(key=lambda row: int(row["node"]))
        nodes = np.asarray([float(row["radius_m"])/scale for row in rows])
        dr = np.asarray([float(row["dr_m"])/scale for row in rows])
        local_maps = []
        for mode in modes:
            offset, size = offsets[(mode.family, mode.degree, mode.order)]
            field = mode.kind
            local = [dof_lookup[(mode.family, mode.degree, element_index,
                                 int(row["node"]), field)] for row in rows]
            if any(index >= size for index in local):
                raise ValueError("exported radial DOF index exceeds its H0 block")
            local_maps.append(offset+np.asarray(local, dtype=np.int64))
        local_field_maps = {}
        for family, degree in template:
            offset, size = offsets[(family, degree, 0)]
            fields = {"R": ("U", "phi"),
                      "S": ("U", "V", "phi"),
                      "T": ("W",)}[family]
            for field in fields:
                local = [dof_lookup[(family, degree, element_index,
                                     int(row["node"]), field)] for row in rows]
                if any(index >= size for index in local):
                    raise ValueError("exported field map exceeds its H0 block")
                local_field_maps[(family, degree, field)] = (
                    offset+np.asarray(local, dtype=np.int64))
        elements.append(Element(element_index, nodes, dr,
                                _barycentric_derivative(nodes), tuple(local_maps),
                                local_field_maps))

    rhs = {bin_index: np.zeros(cursor, dtype=np.complex128) for bin_index in frequency_bins}
    lift = {bin_index: np.zeros(cursor, dtype=np.complex128) for bin_index in frequency_bins}
    for row in _csv_rows(root/"source.csv"):
        bin_index = int(row["bin"])
        block = offsets[(row["family"], int(row["degree"]), int(row["order"]))]
        dof = block[0]+int(row["dof"])
        rhs[bin_index][dof] += complex(float(row["rhs_real"]), float(row["rhs_imag"]))
        lift[bin_index][dof] += complex(float(row["lift_real"]), float(row["lift_imag"]))
    receiver = {bin_index: np.zeros((3, cursor), dtype=np.complex128)
                for bin_index in frequency_bins}
    for row in _csv_rows(root/"receiver.csv"):
        bin_index = int(row["bin"])
        block = offsets[(row["family"], int(row["degree"]), int(row["order"]))]
        dof = block[0]+int(row["dof"])
        component = int(row["component"])
        receiver[bin_index][component, dof] += complex(
            float(row["coefficient_real"]), float(row["coefficient_imag"]))
    if local_matrices:
        expected_local = {(k, family, degree)
                          for k in frequency_bins
                          for family, degree in template}
        if set(local_matrices) != expected_local:
            raise ValueError("source-local H0 exports do not cover every modal block and bin")
        nq = int(metadata["source_local_nq"])
        for (_, family, _), matrix in local_matrices.items():
            fields = 2 if family == "R" else (3 if family == "S" else 1)
            if matrix.shape != (fields*nq, fields*nq):
                raise ValueError("source-local H0 matrix has inconsistent element dimensions")
    source_element = int(metadata.get("source_local_element", metadata["source_element"]))
    return ModalOracle(tuple(modes), tuple(elements), offsets, matrices,
                       local_matrices, source_element, rhs, lift, receiver,
                       metadata, cursor)


def _windowed_spectrum(raw: np.ndarray, length_scale_m: float,
                       time_scale_s: float, epsilon_s_inv: float,
                       nt: int = 128, dt_s: float = 20.0,
                       dt_model_s: float | None = None,
                       df_model: float | None = None,
                       epsilon_model: float | None = None) -> np.ndarray:
    """Match WindowedDisplacementMap.forward on the full real FFT grid."""
    if raw.shape != (3, nt//2+1):
        raise ValueError("raw receiver spectrum must use the full nonnegative grid")
    time_model = np.fft.irfft(raw, n=nt, axis=1)*nt
    dt_model = dt_s/time_scale_s if dt_model_s is None else dt_model_s
    damping_rate = epsilon_s_inv*time_scale_s if epsilon_model is None else epsilon_model
    df_model = 1.0/(nt*dt_model) if df_model is None else df_model
    dt_si = dt_model*time_scale_s
    damping = np.exp(damping_rate*dt_model*np.arange(nt))
    window = np.zeros(nt)
    nwindow = 91
    window[:nwindow] = 0.5-0.5*np.cos(2*pi*np.arange(nwindow)/(nwindow-1))
    time = time_model*df_model*damping[None, :]
    time *= length_scale_m*window[None, :]
    return np.fft.rfft(time, n=nt, axis=1)*dt_si


def _window_adjoint(weights: np.ndarray, length_scale_m: float,
                    time_scale_s: float, epsilon_s_inv: float,
                    nt: int = 128, dt_s: float = 20.0,
                    dt_model_s: float | None = None,
                    df_model: float | None = None,
                    epsilon_model: float | None = None) -> np.ndarray:
    """Match WindowedDisplacementMap.adjoint under Re(x^dagger y)."""
    if weights.shape != (3, nt//2+1):
        raise ValueError("measurement weights must use the full nonnegative grid")
    dt_model = dt_s/time_scale_s if dt_model_s is None else dt_model_s
    dt_si = dt_model*time_scale_s
    df_model = 1.0/(nt*dt_model) if df_model is None else df_model
    damping_rate = epsilon_s_inv*time_scale_s if epsilon_model is None else epsilon_model
    eta = np.full(nt//2+1, 2.0)
    eta[0] = 1.0
    eta[-1] = 1.0
    q = weights/eta[None, :]
    time = np.fft.irfft(q, n=nt, axis=1)*nt
    damping = np.exp(damping_rate*dt_model*np.arange(nt))
    window = np.zeros(nt)
    nwindow = 91
    window[:nwindow] = 0.5-0.5*np.cos(2*pi*np.arange(nwindow)/(nwindow-1))
    time *= dt_si*damping[None, :]*window[None, :]
    spectrum = np.fft.rfft(time, n=nt, axis=1)
    return length_scale_m*df_model*spectrum*eta[None, :]


def solve_modal_baseline(system: ModalOracle, nt: int = 128,
                         dt_s: float = 20.0, epsilon_s_inv: float = 1/360,
                         q_by_objective: dict[str, np.ndarray] | None = None
                         ) -> dict[str, object]:
    """Solve exported independent H0 blocks and evaluate the two objectives."""
    length_scale = float(system.metadata["length_scale_m"])
    time_scale = float(system.metadata["time_scale_s"])
    dt_model = float(system.metadata["dt_model_s"])
    df_model = float(system.metadata["df_model"])
    epsilon_model = float(system.metadata["epsilon_model"])
    raw_by_objective = {"J_Z": np.zeros((3, nt//2+1), dtype=np.complex128),
                        "J_N": np.zeros((3, nt//2+1), dtype=np.complex128)}
    adjoint_by_objective = {}
    displacement_by_bin = {}
    continuous_by_bin = {}
    factors = {}
    for bin_index in (1, 2, 3):
        continuous = np.zeros(system.size, dtype=np.complex128)
        for family, degree, order in system.offsets:
            offset, size = system.offsets[(family, degree, order)]
            matrix = system.matrices[(bin_index, family, degree)]
            factor_key = (bin_index, family, degree)
            if factor_key not in factors:
                factors[factor_key] = _factor_h0(matrix)
            factor = factors[factor_key]
            continuous[offset:offset+size] = _solve_h0_factor(
                factor, system.rhs[bin_index][offset:offset+size])
        displacement = continuous + system.lift[bin_index]
        z = complex(2*pi*bin_index*df_model, -epsilon_model)
        fz = -1j/z
        raw = system.receiver[bin_index]@displacement*fz
        raw_by_objective["J_Z"][:, bin_index] = raw
        raw_by_objective["J_N"][:, bin_index] = raw
        displacement_by_bin[bin_index] = displacement
        continuous_by_bin[bin_index] = continuous
    if q_by_objective is None:
        q = np.zeros((3, nt//2+1), dtype=np.complex128)
        q[0, 1], q[0, 2] = 0.7+0.2j, -0.25+0.45j
        qz = q.copy()
        qn = np.zeros_like(q)
        qn[1] = q[0]
    else:
        qz = np.asarray(q_by_objective["J_Z"], dtype=np.complex128)
        qn = np.asarray(q_by_objective["J_N"], dtype=np.complex128)
        if qz.shape != (3, nt//2+1) or qn.shape != qz.shape:
            raise ValueError("fixture q maps must have shape 3 by nt/2+1")
    result = {}
    for name, weights in (("J_Z", qz), ("J_N", qn)):
        processed = _windowed_spectrum(raw_by_objective[name], length_scale,
                                        time_scale, epsilon_s_inv, nt, dt_s,
                                        dt_model, df_model, epsilon_model)
        signal = float(np.linalg.norm(weights)*np.linalg.norm(processed))
        if not np.isfinite(signal) or signal <= 0:
            raise ValueError(f"baseline objective {name} has zero signal scale")
        result[name] = {"raw": raw_by_objective[name], "processed": processed,
                        "J0": float(np.vdot(weights, processed).real), "S": signal}
        adjoint_by_objective[name] = _window_adjoint(
            weights, length_scale, time_scale, epsilon_s_inv, nt, dt_s,
            dt_model, df_model, epsilon_model)
    return {"objectives": result, "q": {"J_Z": qz, "J_N": qn},
            "qhat": adjoint_by_objective, "h0_factors": factors,
            "displacement": displacement_by_bin,
            "continuous": continuous_by_bin}


def _solve_h0_vector(system: ModalOracle, bin_index: int,
                     vector: np.ndarray, factors: dict,
                     transpose: str | None = None) -> np.ndarray:
    result = np.zeros(system.size, dtype=np.complex128)
    for family, degree, order in system.offsets:
        offset, size = system.offsets[(family, degree, order)]
        factor = factors[(bin_index, family, degree)]
        result[offset:offset+size] = _solve_h0_factor(
            factor, vector[offset:offset+size], transpose)
    return result


def assembled_coupling_structure(system: ModalOracle,
                                 matrix: csr_matrix) -> dict[str, object]:
    """Measure actual sparse S/T and off-diagonal (l,m) operator blocks."""
    coo = matrix.tocoo()
    blocks = list(system.offsets.items())
    ends = np.asarray([offset+size for _, (offset, size) in blocks])
    row_block = np.searchsorted(ends, coo.row, side="right")
    col_block = np.searchsorted(ends, coo.col, side="right")
    row_family = np.asarray([blocks[i][0][0] for i in row_block])
    col_family = np.asarray([blocks[i][0][0] for i in col_block])
    row_degree = np.asarray([blocks[i][0][1] for i in row_block])
    col_degree = np.asarray([blocks[i][0][1] for i in col_block])
    row_order = np.asarray([blocks[i][0][2] for i in row_block])
    col_order = np.asarray([blocks[i][0][2] for i in col_block])
    cross = ((row_family == "S") & (col_family == "T") |
             (row_family == "T") & (col_family == "S"))
    offdiag = (row_degree != col_degree) | (row_order != col_order)
    cross_norm = float(np.linalg.norm(coo.data[cross]))
    offdiag_norm = float(np.linalg.norm(coo.data[offdiag]))
    if cross_norm <= 1.0e-12 or offdiag_norm <= 1.0e-12:
        raise ValueError("assembled Cartesian ΔH lacks required S/T or off-diagonal coupling")
    return {"spheroidal_toroidal_frobenius_norm": cross_norm,
            "offdiagonal_lm_frobenius_norm": offdiag_norm,
            "assembled_delta_h_nnz": int(matrix.nnz)}


def _global_h0(system: ModalOracle, bin_index: int) -> csr_matrix:
    blocks = [system.matrices[(bin_index, family, degree)]
              for family, degree, _ in system.offsets]
    return block_diag(blocks, format="csr", dtype=np.complex128)


def _solve_coupled(system: ModalOracle, bin_index: int, delta_h: csr_matrix,
                   h: float, rhs: np.ndarray, factors: dict,
                   rtol: float = 1.0e-13,
                   true_residual_gate: float = 3.0e-13) -> tuple[np.ndarray, dict]:
    """Solve a sparse coupled modal system with cached DAD H0 factors."""
    h0 = _global_h0(system, bin_index)
    operator = (h0+h*delta_h).tocsr()
    preconditioner = LinearOperator(
        operator.shape,
        matvec=lambda value: _solve_h0_vector(
            system, bin_index, np.asarray(value), factors),
        dtype=np.complex128)
    iterations = [0]
    solution, info = gmres(
        operator, rhs, M=preconditioner, **{_gmres_tolerance_keyword(): rtol},
        atol=0.0,
        restart=40, maxiter=500,
        callback=lambda _: iterations.__setitem__(0, iterations[0]+1),
        callback_type="pr_norm")
    residual = operator@solution-rhs
    relative = float(np.linalg.norm(residual)/max(np.linalg.norm(rhs), 1.0e-300))
    if info != 0 or not np.isfinite(relative) or relative > true_residual_gate:
        raise ValueError(
            f"coupled GMRES failed at h={h:g}, bin={bin_index}: "
            f"info={info}, true residual={relative:.3e}")
    return solution, {"iterations": iterations[0],
                      "true_relative_residual": relative,
                      "preconditioned_info": int(info)}, residual


def _objective_from_displacements(system: ModalOracle,
                                  displacement_by_bin: dict[int, np.ndarray],
                                  q_by_objective: dict[str, np.ndarray],
                                  nt: int, dt_s: float,
                                  epsilon_s_inv: float) -> dict[str, object]:
    length_scale = float(system.metadata["length_scale_m"])
    time_scale = float(system.metadata["time_scale_s"])
    df_model = float(system.metadata["df_model"])
    dt_model = float(system.metadata["dt_model_s"])
    epsilon_model = float(system.metadata["epsilon_model"])
    raw = np.zeros((3, nt//2+1), dtype=np.complex128)
    for bin_index, displacement in displacement_by_bin.items():
        z = complex(2*pi*bin_index*df_model, -epsilon_model)
        raw[:, bin_index] = ((-1j/z)*system.receiver[bin_index]@displacement)
    results = {}
    for objective, q in q_by_objective.items():
        processed = _windowed_spectrum(
            raw, length_scale, time_scale, epsilon_s_inv, nt, dt_s,
            dt_model, df_model, epsilon_model)
        results[objective] = {
            "J": float(np.vdot(q, processed).real),
            "raw": raw, "processed": processed}
    return results


def finite_difference_objectives(
        system: ModalOracle, baseline: dict[str, object],
        projection: dict[str, object], angular_points: np.ndarray,
        moduli_pa: dict[str, float], pressure_scale_pa: float,
        length_scale_m: float, q_by_objective: dict[str, np.ndarray],
        nt: int, dt_s: float, epsilon_s_inv: float,
        support: str, directions: Iterable[str] | None = None,
        memory_cap_bytes: int = 1 << 30,
        progress_path: Path | None = None) -> dict[str, object]:
    """Dormant six-step solve implementation; CLI remains release-gated."""
    if support not in {"regular", "source_local"}:
        raise ValueError("support must be regular or source_local")
    modes = system.modes
    gram = angular_gram(modes, angular_points)
    source_center = float(system.metadata["source_radius_m"])
    if support == "regular":
        support_m, profile = (4.2e6, 6.2e6), regular_weight
        active_elements = system.elements
        selected_directions = tuple(moduli_pa.keys() if directions is None else directions)
    else:
        support_m = (source_center-10_000.0, source_center+10_000.0)
        profile = lambda r: source_local_weight(r, source_center)
        active_elements = system.elements
        selected_directions = tuple(("C", "F", "L") if directions is None else directions)
    nq = len(system.elements[0].nodes_model)
    pair_count = sum(_angularly_coupled(a, b) for a in modes for b in modes)
    active_count = sum(
        float(np.min(e.nodes_model)*length_scale_m) >= support_m[0]-1.0e-7 and
        float(np.max(e.nodes_model)*length_scale_m) <= support_m[1]+1.0e-7
        for e in active_elements)
    delta_nnz_upper = active_count*nq*nq*pair_count
    delta_bytes_upper = 4*delta_nnz_upper*(4+4+16)
    factor_bytes = sum(int(factor.L.nnz+factor.U.nnz)*20
                       for factor, _ in baseline["h0_factors"].values())
    gmres_bytes = (40+6)*system.size*16
    vector_bytes = 12*system.size*16
    unique_blocks = sorted({key[:2] for key in system.offsets})
    h0_nnz_upper = sum(
        system.matrices[(1, family, degree)].nnz*
        sum(1 for key in system.offsets if key[:2] == (family, degree))
        for family, degree in unique_blocks)
    h0_sparse_bytes_upper = sum(
        (matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)*
        sum(1 for key in system.offsets if key[:2] == (family, degree))
        for (bin_index, family, degree), matrix in system.matrices.items()
        if bin_index == 1)
    operator_bytes = 2*h0_sparse_bytes_upper
    estimated_peak = (delta_bytes_upper+factor_bytes+gmres_bytes+
                      vector_bytes+operator_bytes)
    if estimated_peak > memory_cap_bytes:
        raise MemoryError(f"coupled oracle estimate {estimated_peak} exceeds 1 GiB cap")
    base_lift = system.lift
    projection_ratio = {"A": 0.0, "C": moduli_pa["C"]/moduli_pa["C"],
                        "F": 0.0, "L": moduli_pa["L"]/moduli_pa["L"],
                        "N": 0.0}
    raw_results = {}
    residuals = {}
    delta_nnz = {}
    source_local_element = system.elements[system.source_element]
    source_only = tuple(e for e in system.elements
                        if e.index == source_local_element.index)
    for direction in selected_directions:
        modulus = moduli_pa[direction]
        delta_h = assemble_delta_h(
            modes, gram, system.elements, direction, modulus,
            pressure_scale_pa, length_scale_m, system.size,
            support_m=support_m, radial_profile=profile)
        delta_nnz[direction] = int(delta_h.nnz)
        delta_upper = None
        if support == "source_local":
            delta_upper = assemble_delta_h(
                modes, gram, source_only, direction, modulus,
                pressure_scale_pa, length_scale_m, system.size,
                support_m=support_m, radial_profile=profile)
        for h in (0.08, 0.04, 0.02, 0.01, 0.005, 0.0025,
                  -0.08, -0.04, -0.02, -0.01, -0.005, -0.0025):
            displacement_by_bin, solve_metrics = {}, {}
            for bin_index in (1, 2, 3):
                rhs = system.rhs[bin_index].copy()
                if support == "source_local":
                    g, *_ = source_local_lift(
                        system, bin_index, projection, h, direction,
                        projection_ratio[direction])
                    g0 = base_lift[bin_index]
                    explicit0 = rhs+source_element_action(system, bin_index, g0)
                    rhs = explicit0+source_force_delta(
                        system, bin_index, projection, direction, h, modulus,
                        moduli_pa)
                    rhs -= source_element_action(system, bin_index, g)
                    rhs -= h*(delta_upper@g)
                else:
                    g = base_lift[bin_index]
                u, metrics, residual = _solve_coupled(
                    system, bin_index, delta_h, h, rhs,
                    baseline["h0_factors"])
                for objective in ("J_Z", "J_N"):
                    signal = float(baseline["objectives"][objective]["S"])
                    lam0 = baseline["adjoints"][objective][bin_index]
                    metrics[objective] = {
                        "lambda0_adjoint_residual_over_S":
                            float(abs(np.vdot(lam0, residual))/signal)}
                displacement_by_bin[bin_index] = u+g
                solve_metrics[str(bin_index)] = metrics
            objectives = _objective_from_displacements(
                system, displacement_by_bin, q_by_objective,
                nt, dt_s, epsilon_s_inv)
            raw_results[(direction, h)] = objectives
            residuals[f"{direction}_{h:+.6g}"] = solve_metrics
            if progress_path is not None:
                progress_record = {
                    "schema": "dspec-3d-fd-progress-v1",
                    "support": support, "modulus": direction, "h": h,
                    "objectives_m_s": {
                        name: float(value["J"])
                        for name, value in objectives.items()},
                    "solver_residuals": solve_metrics}
                with progress_path.open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps(_plain(progress_record), sort_keys=True)+"\n")
                    stream.flush()
                    os.fsync(stream.fileno())
    return {"objectives": raw_results, "residuals": residuals,
            "resource_estimate": {
                "total_dofs": system.size,
                "structural_nonzeros_H0": int(h0_nnz_upper),
                "active_elements": active_count,
                "allowed_mode_pairs": pair_count,
                "delta_upper_nnz": delta_nnz_upper,
                "structural_nonzeros_DeltaH_by_modulus": delta_nnz,
                "delta_sparse_conversion_bytes": delta_bytes_upper,
                "h0_sparse_operator_bytes_estimate": h0_sparse_bytes_upper,
                "coupled_operator_and_h0_bytes_estimate": operator_bytes,
                "cached_h0_factor_bytes_estimate": factor_bytes,
                "gmres_workspace_bytes_estimate": gmres_bytes,
                "solution_vector_bytes_estimate": vector_bytes,
                "estimated_peak_bytes": estimated_peak,
                "memory_cap_bytes": memory_cap_bytes}}


def regular_born_predictions(system: ModalOracle,
                             baseline: dict[str, object],
                             angular_points: np.ndarray,
                             moduli_pa: dict[str, float],
                             pressure_scale_pa: float,
                             length_scale_m: float,
                             memory_cap_bytes: int = 1 << 30) -> dict[str, object]:
    """Independent full-modal derivative for each regular TI modulus direction."""
    modes = system.modes
    gram = angular_gram(modes, angular_points)
    factors = baseline["h0_factors"]
    active = []
    for element in system.elements:
        radii = np.asarray(element.nodes_model)*length_scale_m
        tol = 8*np.finfo(float).eps*6.2e6
        if np.min(radii) >= 4.2e6-tol and np.max(radii) <= 6.2e6+tol:
            active.append(element)
    if not active:
        raise ValueError("no complete SEM element lies in the regular support")
    nq = len(active[0].nodes_model)
    pair_count = sum(_angularly_coupled(a, b) for a in modes for b in modes)
    delta_upper_nnz = len(active)*nq*nq*pair_count
    delta_coo_conversion_bytes = 4*delta_upper_nnz*(4+4+16)
    factor_bytes = sum(int(factor.L.nnz+factor.U.nnz)*20
                       for factor, _ in factors.values())
    gmres_bytes = (40+6)*system.size*16
    vector_bytes = 12*system.size*16
    estimated_peak = (delta_coo_conversion_bytes+factor_bytes+
                      gmres_bytes+vector_bytes+gram.nbytes)
    if estimated_peak > memory_cap_bytes:
        raise MemoryError(f"oracle peak estimate {estimated_peak} exceeds {memory_cap_bytes} bytes")
    predictions = {objective: {} for objective in ("J_Z", "J_N")}
    delta_nnz = {}
    assembled_structure = {}
    adjoints = {objective: {} for objective in ("J_Z", "J_N")}
    time_scale = float(system.metadata["time_scale_s"])
    df_model = float(system.metadata["df_model"])
    epsilon_model = float(system.metadata["epsilon_model"])
    nt, dt_s, epsilon = 128, 20.0, 1/360
    for bin_index in (1, 2, 3):
        z = complex(2*pi*bin_index*df_model, -epsilon_model)
        fz = -1j/z
        for objective in ("J_Z", "J_N"):
            qhat = baseline["qhat"][objective][:, bin_index]
            rhs_adj = system.receiver[bin_index].conj().T @ (
                np.conjugate(fz)*qhat)
            adjoints[objective][bin_index] = _solve_h0_vector(
                system, bin_index, rhs_adj, factors, "H")
    for direction, modulus in moduli_pa.items():
        delta_h = assemble_delta_h(modes, gram, system.elements, direction,
                                   modulus, pressure_scale_pa, length_scale_m,
                                   system.size)
        delta_nnz[direction] = int(delta_h.nnz)
        if direction in {"L", "N"}:
            assembled_structure[direction] = assembled_coupling_structure(
                system, delta_h)
        for objective in ("J_Z", "J_N"):
            value = 0.0
            for bin_index in (1, 2, 3):
                u0 = baseline["continuous"][bin_index]
                lam = adjoints[objective][bin_index]
                value -= np.vdot(lam, delta_h@u0).real
            predictions[objective][direction] = float(value)
    return {"predictions": predictions, "angular_gram": gram,
            "delta_nnz": delta_nnz,
            "assembled_coupling_structure": assembled_structure,
            "adjoints": adjoints,
            "resource_estimate": {"total_dofs": system.size,
                                  "allowed_mode_pairs": pair_count,
                                  "delta_upper_nnz": delta_upper_nnz,
                                  "delta_coo_conversion_bytes": delta_coo_conversion_bytes,
                                  "h0_factor_bytes_estimate": factor_bytes,
                                  "gmres_basis_bytes_estimate": gmres_bytes,
                                  "estimated_peak_bytes": estimated_peak,
                                  "memory_cap_bytes": memory_cap_bytes}}


def _plain(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()
                if key not in {"h0_factors", "displacement", "raw", "processed",
                               "q", "qhat", "angular_gram", "adjoints"}}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def run_baseline(dspec_dir: str | Path, output_dir: str | Path,
                 fixture_path: str | Path) -> dict[str, object]:
    """Create the review packet for the independent baseline and Born oracle."""
    dspec_dir, output_dir = Path(dspec_dir), Path(output_dir)
    fixture = json.loads(Path(fixture_path).read_text(encoding="utf-8"))
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"oracle output must be new or empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    system = load_modal_oracle(dspec_dir/"oracle")
    solver = fixture["solver"]
    nt, dt_s = int(solver["nt"]), float(solver["dt_s"])
    epsilon = float(solver["epsilon_s_inv"])
    q_by_objective = {}
    for objective in ("J_Z", "J_N"):
        q = np.zeros((3, nt//2+1), dtype=np.complex128)
        for key, pair in fixture["solver"]["objectives"][objective]["nonzero_q"].items():
            component_name, bin_text = key.split(",")
            component = {"Z": 0, "N": 1, "E": 2}[component_name]
            q[component, int(bin_text)] = complex(float(pair[0]), float(pair[1]))
        q_by_objective[objective] = q
    baseline = solve_modal_baseline(system, nt, dt_s, epsilon, q_by_objective)
    theta_x, theta_w = np.polynomial.legendre.leggauss(8)
    nphi = 24
    angular_points = np.asarray([
        (np.arccos(x), 2*pi*j/nphi, weight*2*pi/nphi)
        for x, weight in zip(theta_x, theta_w) for j in range(nphi)
    ], dtype=float)
    rho = float(fixture["model"]["rho_kg_m3"])
    vp = float(fixture["model"]["vp_m_s"])
    vs = float(fixture["model"]["vs_m_s"])
    a0, l0 = rho*vp*vp, rho*vs*vs
    moduli_pa = {"A": a0, "C": a0, "F": a0-2*l0, "L": l0, "N": l0}
    born = regular_born_predictions(
        system, baseline, angular_points, moduli_pa,
        float(system.metadata["pressure_scale_pa"]),
        float(system.metadata["length_scale_m"]))
    baseline["adjoints"] = born["adjoints"]
    source_projection = source_projection_grams(system, angular_points)
    source_born = source_local_born_predictions(
        system, baseline, source_projection, angular_points, moduli_pa,
        float(system.metadata["pressure_scale_pa"]),
        float(system.metadata["length_scale_m"]))
    source_lift_zero_error = 0.0
    for bin_index in (1, 2, 3):
        for direction in ("C", "L"):
            projected, *_ = source_local_lift(
                system, bin_index, source_projection, 0.0, direction)
            source_lift_zero_error = max(
                source_lift_zero_error,
                float(np.max(np.abs(projected-system.lift[bin_index]))))
    if source_lift_zero_error > 2.0e-12:
        raise ValueError("zero-step source projection does not reproduce baseline lift")
    coupling = coupling_structure(system.modes, born["angular_gram"])
    production_rows = _csv_rows(dspec_dir/"integrals.csv")
    production_integrals = {
        (row["objective"], row["modulus"]): float(row["I_regular_m_s"])
        for row in production_rows
    }
    born_bridge = {}
    for objective in ("J_Z", "J_N"):
        signal = float(baseline["objectives"][objective]["S"])
        for modulus, prediction in born["predictions"][objective].items():
            reference = production_integrals[(objective, modulus)]
            scale = max(abs(prediction), abs(reference), 1.0e-14*signal)
            relative_error = abs(prediction-reference)/scale
            informative = abs(reference) > 1.0e-6*signal
            born_bridge[f"{objective}_{modulus}"] = {
                "cartesian_born_m_s_per_h": prediction,
                "dspec_integral_m_s_per_h": reference,
                "informative": informative,
                "relative_error": relative_error if informative else None,
                "scaled_absolute_error": abs(prediction-reference)/signal,
                "passed": (relative_error <= 1.0e-6) if informative else None
            }
            if informative and relative_error > 1.0e-6:
                raise ValueError(
                    f"regular Cartesian Born {objective}/{modulus} differs from "
                    f"DSpec production integral by {relative_error:.3e} > 1e-6")
    source_born_bridge = {}
    source_integral_path = dspec_dir/"source_local_integrals.csv"
    if source_integral_path.exists():
        source_rows = _csv_rows(source_integral_path)
        source_integrals = {(row["objective"], row["modulus"]): row
                            for row in source_rows}
        for objective in ("J_Z", "J_N"):
            signal = float(baseline["objectives"][objective]["S"])
            for modulus in fixture["perturbations"]["source_local"]["moduli"]:
                row = source_integrals[(objective, modulus)]
                components = {}
                for component, column in (
                        ("volume", "I_volume_m_s_per_h"),
                        ("contact", "I_contact_m_s_per_h"),
                        ("total", "Kpred_total_m_s_per_h")):
                    independent = source_born[component][objective][modulus]
                    production = float(row[column])
                    scale = max(abs(independent), abs(production),
                                1.0e-14*signal)
                    error = abs(independent-production)/scale
                    informative = abs(production) > 1.0e-6*signal
                    components[component] = {
                        "cartesian_m_s_per_h": independent,
                        "dspec_m_s_per_h": production,
                        "informative": informative,
                        "relative_error": error if informative else None,
                        "scaled_absolute_error": abs(independent-production)/signal,
                        "passed": error <= 1.0e-6 if informative else None}
                    if informative and error > 1.0e-6:
                        raise ValueError(
                            f"source-local {component} Born {objective}/{modulus} "
                            f"differs from DSpec by {error:.3e} > 1e-6")
                source_born_bridge[f"{objective}_{modulus}"] = components
    factors = baseline["h0_factors"]
    factor_nnz = sum(int(factor.L.nnz+factor.U.nnz)
                     for factor, _ in factors.values())
    factor_bytes = factor_nnz*20
    unique_blocks = sorted({key[:2] for key in system.offsets})
    h0_structural_nnz = sum(
        system.matrices[(1, family, degree)].nnz*
        sum(1 for key in system.offsets if key[:2] == (family, degree))
        for family, degree in unique_blocks)
    h0_sparse_bytes = sum(
        (matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)*
        sum(1 for key in system.offsets if key[:2] == (family, degree))
        for (bin_index, family, degree), matrix in system.matrices.items()
        if bin_index == 1)
    pair_count = int(born["resource_estimate"]["allowed_mode_pairs"])
    nq = len(system.elements[0].nodes_model)
    linear_resources = {}
    resource_bounds = {
        "regular": tuple(fixture["perturbations"]["regular_volume"]["support_m"]),
        "source_local": tuple(fixture["perturbations"]["source_local"]["support_m"])}
    for label, bounds in resource_bounds.items():
        active_count = sum(
            float(np.min(element.nodes_model)*float(system.metadata["length_scale_m"]))
            >= bounds[0]-1.0e-7 and
            float(np.max(element.nodes_model)*float(system.metadata["length_scale_m"]))
            <= bounds[1]+1.0e-7
            for element in system.elements)
        delta_upper_nnz = active_count*nq*nq*pair_count
        delta_conversion_bytes = 4*delta_upper_nnz*(4+4+16)
        gmres_bytes = (40+6)*system.size*16
        vector_bytes = 12*system.size*16
        operator_bytes = 2*h0_sparse_bytes
        peak = (delta_conversion_bytes+factor_bytes+gmres_bytes+
                vector_bytes+operator_bytes)
        if peak > 1 << 30:
            raise MemoryError(
                f"{label} coupled estimate {peak} exceeds frozen 1 GiB cap")
        linear_resources[label] = {
            "total_dofs": system.size,
            "structural_nonzeros_H0": int(h0_structural_nnz),
            "structural_nonzeros_DeltaH_upper_bound": int(delta_upper_nnz),
            "active_elements": int(active_count),
            "estimated_sparse_matrix_bytes": int(delta_conversion_bytes+operator_bytes),
            "estimated_factorization_bytes": int(factor_bytes),
            "gmres_workspace_bytes": int(gmres_bytes),
            "solution_vector_bytes": int(vector_bytes),
            "estimated_peak_bytes": int(peak),
            "memory_cap_bytes": 1 << 30,
            "global_dense_matrix_or_factorization": False}
    csr_bytes = sum(int(matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)
                    for matrix in system.matrices.values())
    def read_spectrum(path: Path, columns: tuple[str, str]) -> np.ndarray:
        rows = _csv_rows(path)
        components = {"Z": 0, "N": 1, "E": 2}
        spectrum = np.zeros((3, int(solver["nt"])//2+1), dtype=np.complex128)
        for row in rows:
            if row["objective"] != "J_Z":
                continue
            c, k = components[row["component"]], int(row["bin"])
            spectrum[c, k] = complex(float(row[columns[0]]), float(row[columns[1]]))
        return spectrum

    raw_reference = read_spectrum(dspec_dir/"spectra.csv",
                                  ("raw_real_m_s", "raw_imag_m_s"))
    processed_reference = read_spectrum(dspec_dir/"processed.csv",
                                        ("real_m_s", "imag_m_s"))
    raw_oracle = (baseline["objectives"]["J_Z"]["raw"]*
                  float(system.metadata["length_scale_m"])*
                  float(system.metadata["time_scale_s"]))
    processed_oracle = baseline["objectives"]["J_Z"]["processed"]
    spectrum_errors = {}
    for label, observed, reference in (("raw", raw_oracle, raw_reference),
                                      ("processed", processed_oracle,
                                       processed_reference)):
        for component, name in enumerate(("Z", "N", "E")):
            denom = max(float(np.linalg.norm(observed[component])),
                        float(np.linalg.norm(reference[component])), 1e-300)
            error = float(np.linalg.norm(observed[component]-reference[component])/denom)
            spectrum_errors[f"{label}_{name}"] = error
            if error > BASELINE_SPECTRUM_REL_TOL:
                raise ValueError(
                    f"independent oracle {label} {name} spectrum error {error:.3e} "
                    f"exceeds predeclared {BASELINE_SPECTRUM_REL_TOL:.1e} guard")

    packet = {
        "schema": "dspec-3d-independent-oracle-baseline-v1",
        "status": "BASELINE ONLY; independent pre-run review required before any +/-h",
        "fixture": str(Path(fixture_path).resolve()),
        "active_bins": [1, 2, 3],
        "modal_block_count": len(system.offsets),
        "mode_count": len(system.modes),
        "total_dofs": system.size,
        "objectives": {
            name: {"J0_m_s": data["J0"], "S_m_s": data["S"]}
            for name, data in baseline["objectives"].items()
        },
        "baseline_spectrum_agreement": {
            "relative_l2_by_quantity_component": spectrum_errors,
            "predeclared_tolerance": BASELINE_SPECTRUM_REL_TOL,
            "passed": True
        },
        "cartesian_strain_fd_check": cartesian_strain_fd_check(),
        "source_projection_checks": {
            **source_projection["check"],
            "zero_step_lift_max_error": source_lift_zero_error,
            "source_element": system.source_element,
            "source_element_bounds_m": system.metadata["source_element_bounds_m"],
            "receiver_outside_source_upper_element": True
        },
        "coupling_structure": {
            "angular_cartesian_gram": coupling,
            "assembled_regular_delta_h": born["assembled_coupling_structure"],
            "assembled_source_local_delta_h":
                source_born["assembled_coupling_structure"]},
        "regular_volume_born_m_s_per_h": born["predictions"],
        "source_local_born_m_s_per_h": source_born,
        "source_local_born_vs_dspec_integral": {
            "relative_tolerance": 1.0e-6,
            "directions": source_born_bridge,
            "production_file_present": source_integral_path.exists()
        },
        "regular_born_vs_dspec_integral": {
            "relative_tolerance": 1.0e-6,
            "directions": born_bridge
        },
        "delta_h_nnz": born["delta_nnz"],
        "resource_estimate": born["resource_estimate"],
        "linear_system_resource_estimate": linear_resources,
        "solver_controls": {
            "baseline": "SparseLU per unique degree/family H0 block with DAD scaling; reused across m",
            "baseline_equilibration": "D_i=1/sqrt(max_j |H0_ij|); factor D H0 D; SuperLU Equil=False",
            "perturbed": {"method": "GMRES", "relative_tolerance": 1e-13,
                          "relative_tolerance_keyword": _gmres_tolerance_keyword(),
                          "scipy_version": scipy.__version__, "atol": 0.0,
                          "restart": 40, "maxiter_cycles": 500,
                          "preconditioner": "cached sparse LU per H0 block",
                          "true_residual_gate": 3e-13,
                          "note": "frozen for later independent review; no perturbed solve was run"}
        },
        "solver_workspace": {
            "sparse_h0_template_bytes": csr_bytes,
            "cached_h0_lu_factor_nonzeros": factor_nnz,
            "dense_global_h0_or_factorization": False
        },
        "metadata": system.metadata
    }
    target = output_dir/"baseline.json"
    target.write_text(json.dumps(_plain(packet), indent=2, sort_keys=True)+"\n",
                      encoding="utf-8")
    return packet


def _richardson_plateau(estimates: list[float], signal: float) -> dict[str, object]:
    if len(estimates) != 5:
        raise ValueError("six fixed h values must produce five Richardson estimates")
    floor = 1.0e-6*signal
    candidates = []
    for start in range(3):
        triple = estimates[start:start+3]
        scale = max(abs(value) for value in triple)
        informative = all(abs(value) > floor for value in triple)
        spread = (max(triple)-min(triple))/scale if scale else float("inf")
        if informative and spread <= 1.0e-6:
            candidates.append((start, triple, spread))
    if candidates:
        start, triple, spread = candidates[-1]
        return {"selected_indices": [start, start+1, start+2],
                "values_m_s_per_h": triple,
                "median_m_s_per_h": float(np.median(triple)),
                "spread_over_max_abs": spread,
                "informative": True, "passed": True}
    weak = any(abs(value) <= floor for value in estimates)
    all_spread = max(estimates)-min(estimates)
    return {"selected_indices": None, "values_m_s_per_h": None,
            "median_m_s_per_h": None,
            "spread_over_S": float(all_spread/signal),
            "informative": not weak, "passed": False}


def _production_comparison(prediction: float, plateau: float | None,
                           signal: float) -> tuple[str, dict[str, object]]:
    """Compare the frozen SEM prediction directly with the selected FD plateau."""
    threshold = 1.0e-6*signal
    informative_class = "informative" if abs(prediction) > threshold else "weak"
    if plateau is None:
        return informative_class, {
            "informative": None,
            "relative_error": None, "scaled_absolute_error": None,
            "tolerance": 1.0e-6, "passed": None}
    scaled_absolute = abs(prediction-plateau)/signal
    informative = informative_class == "informative"
    relative = (abs(prediction-plateau)/max(abs(prediction), abs(plateau),
                                            1.0e-14*signal)
                if informative else None)
    passed = relative <= 1.0e-6 if informative else None
    return informative_class, {
        "informative": informative,
        "relative_error": relative,
        "scaled_absolute_error": scaled_absolute,
        "tolerance": 1.0e-6, "passed": passed}


def run_finite_differences(dspec_dir: str | Path, output_dir: str | Path,
                           fixture_path: str | Path) -> dict[str, object]:
    """Run frozen +/-h systems and save the independent plateau packet."""
    dspec_dir, output_dir = Path(dspec_dir), Path(output_dir)
    target = output_dir/"finite_differences.json"
    progress_path = output_dir/"finite_difference_progress.jsonl"
    if not (output_dir/"baseline.json").is_file():
        raise FileNotFoundError("a successful baseline.json is required before FD")
    if target.exists():
        raise FileExistsError(f"finite-difference output already exists: {target}")
    if progress_path.exists():
        raise FileExistsError(f"finite-difference progress already exists: {progress_path}")
    fixture = json.loads(Path(fixture_path).read_text(encoding="utf-8"))
    system = load_modal_oracle(dspec_dir/"oracle")
    solver = fixture["solver"]
    nt, dt_s = int(solver["nt"]), float(solver["dt_s"])
    epsilon = float(solver["epsilon_s_inv"])
    q_by_objective = {}
    for objective in ("J_Z", "J_N"):
        q = np.zeros((3, nt//2+1), dtype=np.complex128)
        for key, pair in solver["objectives"][objective]["nonzero_q"].items():
            component_name, bin_text = key.split(",")
            component = {"Z": 0, "N": 1, "E": 2}[component_name]
            q[component, int(bin_text)] = complex(float(pair[0]), float(pair[1]))
        q_by_objective[objective] = q
    baseline = solve_modal_baseline(system, nt, dt_s, epsilon, q_by_objective)
    theta_x, theta_w = np.polynomial.legendre.leggauss(8)
    angular_points = np.asarray([
        (np.arccos(x), 2*pi*j/24, weight*2*pi/24)
        for x, weight in zip(theta_x, theta_w) for j in range(24)
    ], dtype=float)
    rho = float(fixture["model"]["rho_kg_m3"])
    vp = float(fixture["model"]["vp_m_s"])
    vs = float(fixture["model"]["vs_m_s"])
    a0, l0 = rho*vp*vp, rho*vs*vs
    moduli_pa = {"A": a0, "C": a0, "F": a0-2*l0, "L": l0, "N": l0}
    pressure = float(system.metadata["pressure_scale_pa"])
    length = float(system.metadata["length_scale_m"])
    regular_born = regular_born_predictions(
        system, baseline, angular_points, moduli_pa, pressure, length)
    baseline["adjoints"] = regular_born["adjoints"]
    projection = source_projection_grams(system, angular_points)
    source_born = source_local_born_predictions(
        system, baseline, projection, angular_points, moduli_pa, pressure, length)
    source_directions = tuple(fixture["perturbations"]["source_local"]["moduli"])
    regular_directions = tuple(fixture["perturbations"]["regular_volume"].get(
        "moduli", ("A", "C", "F", "L", "N")))
    supports = {}
    residuals, resources = {}, {}
    source_integrals = {(row["objective"], row["modulus"]): row
                        for row in _csv_rows(dspec_dir/"source_local_integrals.csv")}
    regular_integrals = {(row["objective"], row["modulus"]): row
                         for row in _csv_rows(dspec_dir/"integrals.csv")}
    for support, directions in (("regular", regular_directions),
                                ("source_local", source_directions)):
        born_values = (regular_born["predictions"] if support == "regular"
                       else source_born["total"])
        production_rows = (regular_integrals if support == "regular"
                           else source_integrals)
        production_column = ("I_regular_m_s" if support == "regular"
                             else "Kpred_total_m_s_per_h")
        active = [direction for direction in directions
                  if any(abs(float(production_rows[(objective, direction)][production_column])) >
                         1.0e-6*float(baseline["objectives"][objective]["S"])
                         for objective in ("J_Z", "J_N"))]
        skipped = [direction for direction in directions if direction not in active]
        solved = finite_difference_objectives(
            system, baseline, projection, angular_points, moduli_pa,
            pressure, length, q_by_objective, nt, dt_s, epsilon,
            support, active, progress_path=progress_path)
        residuals[support] = solved["residuals"]
        resources[support] = solved["resource_estimate"]
        directions_report = {}
        for direction in directions:
            per_objective = {}
            for objective in ("J_Z", "J_N"):
                signal = float(baseline["objectives"][objective]["S"])
                pred = float(born_values[objective][direction])
                production_row = (regular_integrals if support == "regular"
                                  else source_integrals)[(objective, direction)]
                production_column = ("I_regular_m_s" if support == "regular"
                                     else "Kpred_total_m_s_per_h")
                production = float(production_row[production_column])
                informative_class = ("informative" if abs(production) > 1.0e-6*signal
                                     else "weak")
                record = {"signal_scale_S_m_s": signal,
                          "born_m_s_per_h": pred,
                          "informative_born": abs(pred) > 1.0e-6*signal,
                          "production_Kpred_m_s_per_h": production,
                          "informative_class": informative_class,
                          "informative_kernel": informative_class == "informative"}
                if direction not in active:
                    record["status"] = "weak; finite differences skipped by frozen rule"
                    record["D_h_m_s_per_h"] = None
                    record["richardson_m_s_per_h"] = None
                    record["plateau"] = None
                    record["born_comparison"] = {
                        "relative_error": None, "scaled_absolute_error": None,
                        "passed": None}
                    _, record["production_comparison"] = _production_comparison(
                        production, None, signal)
                else:
                    dvalues = []
                    h_objectives = []
                    for h in (0.08, 0.04, 0.02, 0.01, 0.005, 0.0025):
                        plus = float(solved["objectives"][(direction, h)][objective]["J"])
                        minus = float(solved["objectives"][(direction, -h)][objective]["J"])
                        dvalues.append((plus-minus)/(2.0*h))
                        h_objectives.append({"h": h, "J_plus_m_s": plus,
                                             "J_minus_m_s": minus})
                    richardson = [(4.0*dvalues[i+1]-dvalues[i])/3.0
                                  for i in range(5)]
                    plateau = _richardson_plateau(richardson, signal)
                    frozen = plateau["median_m_s_per_h"]
                    _, production_comparison = _production_comparison(
                        production, frozen, signal)
                    informative = (abs(pred) > 1.0e-6*signal and frozen is not None and
                                   abs(frozen) > 1.0e-6*signal)
                    relative = (abs(pred-frozen)/max(abs(pred), abs(frozen),
                                1.0e-14*signal) if informative else None)
                    record.update({
                        "status": "computed",
                        "h_objectives": h_objectives,
                        "D_h_m_s_per_h": dvalues,
                        "richardson_m_s_per_h": richardson,
                        "plateau": plateau,
                        "born_comparison": {
                            "informative": informative,
                            "relative_error": relative,
                            "scaled_absolute_error": (abs(pred-frozen)/signal
                                                       if frozen is not None else None),
                            "passed": (relative <= 1.0e-6
                                       if informative and relative is not None else None)},
                        "production_comparison": production_comparison})
                    if support == "source_local":
                        production = source_integrals[(objective, direction)]
                        record["source_prediction_split_m_s_per_h"] = {
                            "oracle_volume": source_born["volume"][objective][direction],
                            "oracle_contact": source_born["contact"][objective][direction],
                            "production_volume": float(production["I_volume_m_s_per_h"]),
                            "production_contact": float(production["I_contact_m_s_per_h"]),
                            "production_total": float(production["Kpred_total_m_s_per_h"])}
                per_objective[objective] = record
            directions_report[direction] = per_objective
        supports[support] = {"directions": directions_report,
                             "weak_directions_skipped": skipped}
    packet = {
        "schema": "dspec-3d-independent-oracle-finite-differences-v1",
        "status": "FINITE DIFFERENCES COMPUTED; independent review required",
        "fixture": str(Path(fixture_path).resolve()),
        "steps_h": [0.08, 0.04, 0.02, 0.01, 0.005, 0.0025],
        "richardson_rule": "R_i=(4*D(h[i+1])-D(h[i]))/3; choose finest passing informative consecutive triple using only R and S",
        "supports": supports,
        "coupling_structure": {
            "angular_cartesian_gram": coupling_structure(
                system.modes, regular_born["angular_gram"]),
            "assembled_regular_delta_h":
                regular_born["assembled_coupling_structure"],
            "assembled_source_local_delta_h":
                source_born["assembled_coupling_structure"]},
        "solver_controls": {
            "method": "preconditioned GMRES",
            "rtol": 1.0e-13, "atol": 0.0, "restart": 40,
            "relative_tolerance_keyword": _gmres_tolerance_keyword(),
            "scipy_version": scipy.__version__,
            "maxiter_cycles": 500, "true_relative_residual_gate": 3.0e-13,
            "preconditioner": "cached symmetric DAD SparseLU factors of H0; SuperLU Equil=False",
            "residual_diagnostic": "save |lambda0^H r|/S; diagnostic only, not a bound for H(h)"},
        "solver_residuals": residuals,
        "resources": resources,
        "memory_cap_bytes": 1 << 30}
    target.write_text(json.dumps(_plain(packet), indent=2, sort_keys=True)+"\n",
                      encoding="utf-8")
    return packet


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path,
                        help="DSpec adapter output directory containing oracle/")
    parser.add_argument("--output", type=Path,
                        help="new or empty directory for baseline.json")
    parser.add_argument("--fixture", type=Path, required=True,
                        help="resolved fixture.json used by the run")
    parser.add_argument("--finite-differences", action="store_true",
                        help="run the frozen six +/-h cases after approved pre-run review")
    args = parser.parse_args()
    if args.finite_differences:
        if args.baseline is None or args.output is None:
            parser.error("--baseline and --output are required for finite differences")
        packet = run_finite_differences(args.baseline, args.output, args.fixture)
        print(json.dumps({"status": packet["status"],
                          "output": str((args.output/"finite_differences.json").resolve()),
                          "progress": str((args.output/"finite_difference_progress.jsonl").resolve()),
                          "supports": list(packet["supports"])}, indent=2))
        return 0
    if args.baseline is None or args.output is None:
        parser.error("--baseline and --output are required for baseline mode")
    packet = run_baseline(args.baseline, args.output, args.fixture)
    print(json.dumps({"status": packet["status"],
                      "output": str((args.output/"baseline.json").resolve()),
                      "total_dofs": packet["total_dofs"],
                      "objectives": packet["objectives"]}, indent=2))
    return 0


def _associated_legendre_2(l: int, m: int, x: float) -> tuple[float, float, float]:
    """Return P_l^m(x), dP/dx, d2P/dx2 with Condon-Shortley phase."""
    if l < 0 or m < 0 or m > l or not (-1.0 < x < 1.0):
        raise ValueError("associated Legendre evaluation requires |x|<1 and 0<=m<=l")
    p, px, pxx = 1.0, 0.0, 0.0
    s = np.sqrt(1.0 - x*x)
    for k in range(1, m + 1):
        a = -(2.0*k - 1.0)
        p, px, pxx = (a*s*p,
                       a*(s*px - x/s*p),
                       a*(s*pxx - 2.0*x/s*px - p/s**3))
    if l == m:
        return p, px, pxx
    p1 = (2*m + 1)*x*p
    p1x = (2*m + 1)*(p + x*px)
    p1xx = (2*m + 1)*(2*px + x*pxx)
    if l == m + 1:
        return p1, p1x, p1xx
    for n in range(m + 2, l + 1):
        a = (2*n - 1)/(n - m)
        b = (n + m - 1)/(n - m)
        nxt = (a*x*p1 - b*p,
               a*(p1 + x*p1x) - b*px,
               a*(2*p1x + x*p1xx) - b*pxx)
        p, px, pxx = p1, p1x, p1xx
        p1, p1x, p1xx = nxt
    return p1, p1x, p1xx


def scalar_harmonic(l: int, m: int, theta: float, phi: float) -> np.ndarray:
    """[Y, Y_theta, Y_phi, Y_thetatheta, Y_thetaphi, Y_phiphi]."""
    if abs(m) > l or not (0.0 < theta < pi):
        raise ValueError("invalid scalar-harmonic arguments")
    am = abs(m)
    x, s = np.cos(theta), np.sin(theta)
    p, px, pxx = _associated_legendre_2(l, am, float(x))
    norm = np.sqrt((2*l + 1)/(4*pi) * factorial(l-am)/factorial(l+am))
    phase = np.exp(1j*am*phi)
    y = norm*p*phase
    yt = -norm*s*px*phase
    ytt = norm*(-x*px + s*s*pxx)*phase
    yp = 1j*am*y
    ytp = 1j*am*yt
    ypp = -(am*am)*y
    if m < 0:
        sign = -1.0 if am % 2 else 1.0
        return sign*np.conjugate(np.array([y, yt, yp, ytt, ytp, ypp]))
    return np.array([y, yt, yp, ytt, ytp, ypp])


def _spherical_strain(mode: Mode, amplitude: int, theta: float,
                      phi: float) -> np.ndarray:
    """Strain from a unit coefficient of (U/r,U',V/r,V',W/r,W')."""
    y, yt, yp, ytt, ytp, ypp = scalar_harmonic(
        mode.degree, mode.order, theta, phi)
    U = Up = V = Vp = W = Wp = 0.0j
    if mode.kind == "U":
        if amplitude == 0: U = 1.0
        else: Up = 1.0
    elif mode.kind == "V":
        if amplitude == 0: V = 1.0
        else: Vp = 1.0
    elif mode.kind == "W":
        if amplitude == 0: W = 1.0
        else: Wp = 1.0
    else:
        raise ValueError(f"unsupported displacement family {mode.kind!r}")
    st, ct = np.sin(theta), np.cos(theta)
    cot, csc = ct/st, 1.0/st
    if mode.kind == "W":
        ur = urt = urp = 0.0j
        ut, utp = -W*csc*yp, -Wp*csc*yp
        up, upp = W*yt, Wp*yt
        ut_t = W*(cot*csc*yp - csc*ytp)
        up_t = W*ytt
        ut_p = -W*csc*ypp
        up_p = W*ytp
    else:
        ur, urt, urp = U*y, U*yt, U*yp
        ut, utp = V*yt, Vp*yt
        up, upp = V*csc*yp, Vp*csc*yp
        ut_t = V*ytt
        # Differentiate csc(theta)*Y_phi completely.  The cot*csc term
        # contributes to epsilon_theta_phi and is needed for N sensitivity.
        up_t = V*(csc*ytp - cot*csc*yp)
        ut_p, up_p = V*ytp, V*csc*ypp
    strain = np.array([Up*y,
                       ur + ut_t,
                       ur + cot*ut + csc*up_p,
                       0.5*(utp - ut + urt),
                       0.5*(upp - up + csc*urp),
                       0.5*(up_t - cot*up + csc*ut_p)], dtype=np.complex128)
    if mode.kind == "W":
        # The formulas above use +rhat x grad(Y); the frozen convention is
        # B_T = -rhat x grad(Y), matching the physical north/east rows.
        strain = -strain
    return strain


def _cartesian_strain(spherical: np.ndarray, theta: float,
                      phi: float) -> np.ndarray:
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    Q = np.array([[st*cp, ct*cp, -sp],
                  [st*sp, ct*sp, cp],
                  [ct, -st, 0.0]])
    E = np.array([[spherical[0], spherical[3], spherical[4]],
                  [spherical[3], spherical[1], spherical[5]],
                  [spherical[4], spherical[5], spherical[2]]])
    return Q @ E @ Q.T


def cartesian_strain_fd_check() -> dict[str, object]:
    """Check spherical strain components against Cartesian displacement FD."""
    theta, phi, radius = 0.83, 0.41, 1.0
    direction = np.array([np.sin(theta)*np.cos(phi),
                          np.sin(theta)*np.sin(phi), np.cos(theta)])
    point = radius*direction
    step = 2.0e-5
    maxima = []
    for degree, order in ((2, 1), (3, -2)):
        for kind in ("U", "V", "W"):
            mode = Mode(degree, order, kind)
            for amplitude in (0, 1):
                def displacement(position: np.ndarray) -> np.ndarray:
                    r = float(np.linalg.norm(position))
                    th = float(np.arccos(position[2]/r))
                    ph = float(np.arctan2(position[1], position[0]))
                    y, yt, yp, *_ = scalar_harmonic(degree, order, th, ph)
                    st = np.sin(th)
                    dr = r-radius
                    if kind == "U":
                        ur = (1.0 if amplitude == 0 else dr)*y
                        uth = uph = 0.0j
                    elif kind == "V":
                        v = 1.0 if amplitude == 0 else dr
                        ur, uth, uph = 0.0j, v*yt, v*yp/st
                    else:
                        w = 1.0 if amplitude == 0 else dr
                        ur, uth, uph = 0.0j, w*yp/st, -w*yt
                    er = np.array([np.sin(th)*np.cos(ph),
                                   np.sin(th)*np.sin(ph), np.cos(th)])
                    eth = np.array([np.cos(th)*np.cos(ph),
                                    np.cos(th)*np.sin(ph), -np.sin(th)])
                    eph = np.array([-np.sin(ph), np.cos(ph), 0.0])
                    return ur*er + uth*eth + uph*eph

                gradient = np.empty((3, 3), dtype=np.complex128)
                for axis in range(3):
                    delta = np.zeros(3)
                    delta[axis] = step
                    gradient[:, axis] = (displacement(point+delta)-
                                         displacement(point-delta))/(2*step)
                finite_difference = 0.5*(gradient+gradient.T)
                predicted = _cartesian_strain(
                    _spherical_strain(mode, amplitude, theta, phi), theta, phi)
                maxima.append(float(np.max(np.abs(finite_difference-predicted))))
    error = max(maxima)
    if error > 2.0e-8:
        raise ValueError(f"Cartesian strain finite-difference check failed: {error:.3e}")
    return {"sample_count": len(maxima), "step": step,
            "max_absolute_component_error": error,
            "tolerance": 2.0e-8, "passed": True}


def _angularly_coupled(test: Mode, trial: Mode) -> bool:
    """Exact Gaunt support for chi's (L,M)=(2,1),(3,2) components."""
    dm = abs(test.order-trial.order)
    dl = abs(test.degree-trial.degree)
    total_degree = test.degree+trial.degree
    return dm in (1, 2) and any(dl <= degree <= total_degree and dm <= degree
                                  for degree in (2, 3))


def ti_tensor(direction: str, theta: float, phi: float) -> np.ndarray:
    """Cartesian C_ijkl for a unit-Pa increment of one TI modulus."""
    if direction not in {"A", "C", "F", "L", "N"}:
        raise ValueError(f"unknown TI modulus {direction!r}")
    n = np.array([np.sin(theta)*np.cos(phi),
                  np.sin(theta)*np.sin(phi), np.cos(theta)])
    P = np.eye(3) - np.outer(n, n)
    a, c, f, ell, nn = (direction == "A", direction == "C",
                        direction == "F", direction == "L",
                        direction == "N")
    out = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    out[i,j,k,l] = (
                        c*n[i]*n[j]*n[k]*n[l]
                        + f*(n[i]*n[j]*P[k,l] + P[i,j]*n[k]*n[l])
                        + (a - 2*nn)*P[i,j]*P[k,l]
                        + nn*(P[i,k]*P[j,l] + P[i,l]*P[j,k])
                        + ell*(n[i]*n[k]*P[j,l] + n[i]*n[l]*P[j,k]
                               + n[j]*n[k]*P[i,l] + n[j]*n[l]*P[i,k])
                    )
    return out


def cartesian_bilinear(test: np.ndarray, trial: np.ndarray,
                       direction: str, theta: float, phi: float) -> complex:
    return np.einsum('ij,ijkl,kl->', np.conjugate(test),
                     ti_tensor(direction, theta, phi), trial,
                     optimize=True)


def angular_gram(modes: Iterable[Mode], points: np.ndarray) -> np.ndarray:
    """Return [mode,mode,A/C/F/L/N,2,2] χ-weighted Cartesian contractions."""
    modes = tuple(modes)
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("angular points must have theta, phi, dOmega columns")
    gram = np.zeros((len(modes), len(modes), 5, 2, 2), dtype=np.complex128)
    directions = ("A", "C", "F", "L", "N")
    for theta, phi, domega in points:
        if domega <= 0 or not np.isfinite(domega):
            raise ValueError("angular weights must be finite and positive")
        chi = scalar_harmonic(2, 1, theta, phi)[0].real + 0.5*scalar_harmonic(3, 2, theta, phi)[0].real
        strain = np.asarray([
            _cartesian_strain(_spherical_strain(mode, amp, theta, phi),
                              theta, phi)
            for mode in modes for amp in range(2)
        ], dtype=np.complex128).reshape(2*len(modes), 9)
        conjugate_strain = np.conjugate(strain)
        weight = domega*chi
        for im, name in enumerate(directions):
            tensor = ti_tensor(name, theta, phi).reshape(9, 9)
            # Two small BLAS products replace mode-pair Python loops and 92
            # repeated rank-four contractions at each angular point.
            block = (conjugate_strain @ tensor) @ strain.T
            gram[:, :, im] += weight*block.reshape(len(modes), 2,
                                                    len(modes), 2).transpose(0, 2, 1, 3)
    return gram


def regular_weight(radius_m: float, r0_m: float = 4.2e6,
                   r1_m: float = 6.2e6) -> float:
    tol = 8*np.finfo(float).eps*max(abs(r0_m), abs(r1_m))
    if radius_m <= r0_m + tol or radius_m >= r1_m - tol:
        return 0.0
    if radius_m < r0_m or radius_m > r1_m:
        return 0.0
    return float(np.sin(pi*(radius_m-r0_m)/(r1_m-r0_m))**2)


def source_local_weight(radius_m: float, center_m: float,
                        halfwidth_m: float = 10_000.0) -> float:
    x = (radius_m-center_m)/halfwidth_m
    if abs(x) >= 1.0:
        return 0.0
    return float(np.cos(0.5*pi*x)**2)


def source_projection_grams(system: ModalOracle,
                            angular_points: np.ndarray) -> dict[str, object]:
    """Scalar and S/T vector harmonic projections for the finite source jump."""
    scalar = [(i, mode) for i, mode in enumerate(system.modes)
              if (mode.family, mode.kind) in (("R", "U"), ("S", "U"))]
    vector = [(i, mode) for i, mode in enumerate(system.modes)
              if (mode.family, mode.kind) in (("S", "V"), ("T", "W"))]
    xs = np.zeros((len(scalar), len(scalar)), dtype=np.complex128)
    xv = np.zeros((len(vector), len(vector)), dtype=np.complex128)
    gram = np.diag([mode.degree*(mode.degree+1.0) for _, mode in vector])
    for theta, phi, domega in np.asarray(angular_points, dtype=float):
        chi = (scalar_harmonic(2, 1, theta, phi)[0].real +
               0.5*scalar_harmonic(3, 2, theta, phi)[0].real)
        yvalues = [scalar_harmonic(mode.degree, mode.order, theta, phi)[0]
                   for _, mode in scalar]
        for i, yi in enumerate(yvalues):
            for j, yj in enumerate(yvalues):
                xs[i, j] += domega*chi*np.conjugate(yi)*yj
        bvalues = []
        csc = 1.0/np.sin(theta)
        for _, mode in vector:
            _, yt, yp, *_ = scalar_harmonic(
                mode.degree, mode.order, theta, phi)
            if mode.kind == "V":
                bvalues.append(np.array([yt, csc*yp]))
            else:
                bvalues.append(np.array([csc*yp, -yt]))
        for i, bi in enumerate(bvalues):
            for j, bj in enumerate(bvalues):
                xv[i, j] += domega*chi*np.vdot(bi, bj)
    gram_check = np.zeros_like(gram)
    for theta, phi, domega in np.asarray(angular_points, dtype=float):
        bvalues = []
        csc = 1.0/np.sin(theta)
        for _, mode in vector:
            _, yt, yp, *_ = scalar_harmonic(
                mode.degree, mode.order, theta, phi)
            bvalues.append(np.array([yt, csc*yp]) if mode.kind == "V" else
                           np.array([csc*yp, -yt]))
        for i, bi in enumerate(bvalues):
            for j, bj in enumerate(bvalues):
                gram_check[i, j] += domega*np.vdot(bi, bj).real
    scalar_hermitian = float(np.max(np.abs(xs-xs.conj().T)))
    vector_hermitian = float(np.max(np.abs(xv-xv.conj().T)))
    vector_gram_error = float(np.max(np.abs(gram_check-gram)))
    if max(scalar_hermitian, vector_hermitian, vector_gram_error) > 2.0e-11:
        raise ValueError("source harmonic projection Gram check failed")
    sidx = [i for i, (_, mode) in enumerate(vector) if mode.family == "S"]
    tidx = [i for i, (_, mode) in enumerate(vector) if mode.family == "T"]
    vector_cross = float(np.linalg.norm(xv[np.ix_(sidx, tidx)]))
    if vector_cross <= 1.0e-12:
        raise ValueError("source-local vector projection has no S/T coupling")
    return {"scalar": scalar, "vector": vector, "X_scalar": xs,
            "X_vector": xv, "G_vector": gram,
            "check": {"scalar_hermitian_error": scalar_hermitian,
                      "vector_hermitian_error": vector_hermitian,
                      "vector_base_gram_max_error": vector_gram_error,
                      "source_vector_ST_coupling_norm": vector_cross}}


def coupling_structure(modes: tuple[Mode, ...], gram: np.ndarray
                       ) -> dict[str, object]:
    result = {}
    sidx = [i for i, mode in enumerate(modes) if mode.family == "S"]
    tidx = [i for i, mode in enumerate(modes) if mode.family == "T"]
    offdiagonal = np.asarray([
        [(a.degree, a.order) != (b.degree, b.order) for b in modes]
        for a in modes], dtype=bool)
    for index, name in ((3, "L"), (4, "N")):
        cross = float(np.linalg.norm(gram[np.ix_(sidx, tidx, [index],
                                                  range(2), range(2))]))
        offdiag = float(np.linalg.norm(gram[:, :, index, :, :][offdiagonal]))
        if cross <= 1.0e-12 or offdiag <= 1.0e-12:
            raise ValueError(f"{name} Cartesian ΔH is missing required modal coupling")
        result[name] = {"spheroidal_toroidal_block_norm": cross,
                        "offdiagonal_lm_block_norm": offdiag}
    return result


def _element_field_dofs(system: ModalOracle, element: Element,
                        mode: Mode, field: str) -> np.ndarray:
    base = element.field_dofs[(mode.family, mode.degree, field)]
    mode_offset = system.offsets[(mode.family, mode.degree, mode.order)][0]
    base_offset = system.offsets[(mode.family, mode.degree, 0)][0]
    return base + mode_offset-base_offset


def source_local_lift(system: ModalOracle, bin_index: int,
                      projection: dict[str, object], h: float,
                      direction: str, modulus_ratio: float = 1.0
                      ) -> tuple[np.ndarray, np.ndarray,
                                                np.ndarray, np.ndarray,
                                                np.ndarray]:
    """Project source jumps, then embed only on the owning upper element."""
    if not system.source_local_matrices:
        raise ValueError("source-local baseline H_e exports are required")
    element = system.elements[system.source_element]
    base_lift = system.lift[bin_index]
    scalar = projection["scalar"]
    vector = projection["vector"]
    x_scalar = projection["X_scalar"]
    x_vector = projection["X_vector"]
    g_vector = projection["G_vector"]
    scalar_dofs = [_element_field_dofs(system, element, mode, "U")[0]
                   for _, mode in scalar]
    vector_field = ["V" if mode.family == "S" else "W" for _, mode in vector]
    vector_dofs = [_element_field_dofs(system, element, mode, field)[0]
                   for (_, mode), field in zip(vector, vector_field)]
    jr0 = base_lift[scalar_dofs].copy()
    jh0 = base_lift[vector_dofs].copy()
    jr = jr0
    jh = jh0
    if direction == "C":
        jr = np.linalg.solve(np.eye(len(jr0))+h*modulus_ratio*x_scalar, jr0)
    if direction == "L":
        jh = np.linalg.solve(g_vector+h*modulus_ratio*x_vector,
                             g_vector@jh0)
    lift = np.zeros(system.size, dtype=np.complex128)
    lift[scalar_dofs] = jr
    lift[vector_dofs] = jh
    return lift, jr0, jh0, jr, jh


def source_element_action(system: ModalOracle, bin_index: int,
                          vector: np.ndarray) -> np.ndarray:
    """Apply H_e only through the source-owning upper element rows."""
    element = system.elements[system.source_element]
    fields_by_family = {"R": ("U", "phi"),
                        "S": ("U", "V", "phi"),
                        "T": ("W",)}
    result = np.zeros(system.size, dtype=np.complex128)
    nq = len(element.nodes_model)
    for family, degree, order in system.offsets:
        local = np.zeros(len(fields_by_family[family])*nq,
                         dtype=np.complex128)
        dof_maps = []
        for field_index, field in enumerate(fields_by_family[family]):
            dofs = _element_field_dofs(
                system, element, Mode(degree, order, family_kind(family), family), field)
            dof_maps.append(dofs)
            local[field_index*nq:(field_index+1)*nq] = vector[dofs]
        matrix = system.source_local_matrices[(bin_index, family, degree)]
        action = matrix@local
        for field_index, dofs in enumerate(dof_maps):
            result[dofs] += action[field_index*nq:(field_index+1)*nq]
    return result


def source_local_born_predictions(
        system: ModalOracle, baseline: dict[str, object],
        projection: dict[str, object], angular_points: np.ndarray,
        moduli_pa: dict[str, float], pressure_scale_pa: float,
        length_scale_m: float) -> dict[str, object]:
    """Independent volume/contact derivative for the upper source patch."""
    modes = system.modes
    gram = angular_gram(modes, angular_points)
    center = float(system.metadata["source_radius_m"])
    bounds = tuple(float(value) for value in system.metadata[
        "source_element_bounds_m"])
    support = (center-10_000.0, center+10_000.0)
    if not np.allclose(bounds, (center, center+10_000.0), rtol=0.0, atol=1e-7):
        raise ValueError("source-owning upper element is not the upper support half")
    pressure_scale = float(system.metadata["pressure_scale_pa"])
    length_scale = float(system.metadata["length_scale_m"])
    volume = {objective: {} for objective in ("J_Z", "J_N")}
    contact = {objective: {} for objective in ("J_Z", "J_N")}
    total = {objective: {} for objective in ("J_Z", "J_N")}
    assembled_structure = {}
    baseline_moduli = {"A": moduli_pa["A"], "C": moduli_pa["C"],
                       "F": moduli_pa["F"], "L": moduli_pa["L"],
                       "N": moduli_pa["N"]}
    source_element = system.elements[system.source_element]
    source_only = tuple(element for element in system.elements
                        if element.index == source_element.index)
    for direction, modulus in moduli_pa.items():
        delta_global = assemble_delta_h(
            modes, gram, system.elements, direction, modulus, pressure_scale_pa,
            length_scale_m, system.size, support_m=support,
            radial_profile=lambda r: source_local_weight(r, center))
        delta_upper = assemble_delta_h(
            modes, gram, source_only, direction, modulus, pressure_scale_pa,
            length_scale_m, system.size, support_m=support,
            radial_profile=lambda r: source_local_weight(r, center))
        if direction in {"L", "N"}:
            assembled_structure[direction] = assembled_coupling_structure(
                system, delta_global)
        for objective in ("J_Z", "J_N"):
            v, c = 0.0, 0.0
            for bin_index in (1, 2, 3):
                u0 = baseline["continuous"][bin_index]
                g0 = system.lift[bin_index]
                lam = baseline["adjoints"][objective][bin_index]
                v -= np.vdot(lam, delta_global@u0).real
                v -= np.vdot(lam, delta_upper@g0).real
                if direction in {"C", "F"}:
                    fprime = _source_force_derivative(
                        system, bin_index, projection, direction, modulus,
                        baseline_moduli)
                else:
                    fprime = np.zeros(system.size, dtype=np.complex128)
                gprime = source_lift_derivative(
                    system, bin_index, projection, direction,
                    modulus/baseline_moduli.get(direction, 1.0))
                c += np.vdot(lam, fprime-source_element_action(
                    system, bin_index, gprime)).real
                qhat = baseline["qhat"][objective][:, bin_index]
                z = complex(2*pi*bin_index*float(system.metadata["df_model"]),
                            -float(system.metadata["epsilon_model"]))
                direct_raw = (-1j/z)*(system.receiver[bin_index]@gprime)
                c += np.vdot(qhat, direct_raw).real
            volume[objective][direction] = float(v)
            contact[objective][direction] = float(c)
            total[objective][direction] = float(v+c)
    return {"volume": volume, "contact": contact, "total": total,
            "source_support_m": list(support),
            "source_upper_element": source_element.index,
            "source_projection": {
                "scalar": "(I+h*(deltaC/C0)*X_scalar) j_r(h)=j_r0",
                "vector": "(G+h*(deltaL/L0)*X_vector) j_h(h)=G j_h0",
                "G_vector": "diag(l*(l+1)) over S and T modes",
                "F_source_model": float(moduli_pa["F"]/pressure_scale),
                "source_radius_model": center/length_scale,
                "explicit_F_load": "(-2,+lambda) r_model*(F_source_model*j_r) with chi multiplication X_scalar; perturbation Pa converted by pressure_scale",
                "contact_scatter": "upper source-owning element rows only"},
            "assembled_coupling_structure": assembled_structure,
            "direct_receiver_is_zero_for_fixture": True}


def source_lift_derivative(system: ModalOracle, bin_index: int,
                           projection: dict[str, object], direction: str,
                           modulus_ratio: float) -> np.ndarray:
    """Derivative of the projected source lift at h=0."""
    result = np.zeros(system.size, dtype=np.complex128)
    source = system.elements[system.source_element]
    if direction == "C":
        scalar = projection["scalar"]
        x = projection["X_scalar"]
        base = system.lift[bin_index]
        j0 = np.asarray([base[_element_field_dofs(system, source, mode, "U")[0]]
                         for _, mode in scalar])
        derivative = -modulus_ratio*(x@j0)
        for value, (_, mode) in zip(derivative, scalar):
            dof = _element_field_dofs(system, source, mode, "U")[0]
            result[dof] = value
    elif direction == "L":
        vector = projection["vector"]
        x, g = projection["X_vector"], projection["G_vector"]
        base = system.lift[bin_index]
        j0 = np.asarray([base[_element_field_dofs(
            system, source, mode, "V" if mode.family == "S" else "W")[0]]
                         for _, mode in vector])
        derivative = -modulus_ratio*np.linalg.solve(g, x@j0)
        for value, (_, mode) in zip(derivative, vector):
            field = "V" if mode.family == "S" else "W"
            dof = _element_field_dofs(system, source, mode, field)[0]
            result[dof] = value
    return result


def _source_force_derivative(system: ModalOracle, bin_index: int,
                             projection: dict[str, object], direction: str,
                             modulus_pa: float,
                             baseline_moduli_pa: dict[str, float]) -> np.ndarray:
    """Derivative of the explicit F*j_r source rows for C or F perturbations."""
    result = np.zeros(system.size, dtype=np.complex128)
    scalar = projection["scalar"]
    source = system.elements[system.source_element]
    base = system.lift[bin_index]
    jr0 = np.asarray([base[_element_field_dofs(system, source, mode, "U")[0]]
                      for _, mode in scalar])
    ratio = modulus_pa/baseline_moduli_pa[direction]
    pressure_scale = float(system.metadata["pressure_scale_pa"])
    if direction == "C":
        delta_fj_model = -ratio*(baseline_moduli_pa["F"]/pressure_scale)*(
            projection["X_scalar"]@jr0)
    elif direction == "F":
        delta_fj_model = ratio*(baseline_moduli_pa["F"]/pressure_scale)*(
            projection["X_scalar"]@jr0)
    else:
        return result
    rsource_model = (float(system.metadata["source_radius_m"])/
                     float(system.metadata["length_scale_m"]))
    for value, (_, mode) in zip(delta_fj_model, scalar):
        udof = _element_field_dofs(system, source, mode, "U")[0]
        result[udof] -= 2.0*rsource_model*value
        if mode.family == "S":
            vdof = _element_field_dofs(system, source, mode, "V")[0]
            result[vdof] += mode.degree*(mode.degree+1.0)*rsource_model*value
    return result


def source_force_delta(system: ModalOracle, bin_index: int,
                       projection: dict[str, object], direction: str,
                       h: float, modulus_pa: float,
                       baseline_moduli_pa: dict[str, float]) -> np.ndarray:
    """Reproject the explicit C/F source load without changing the event."""
    result = np.zeros(system.size, dtype=np.complex128)
    if direction not in {"C", "F"} or h == 0.0:
        return result
    scalar = projection["scalar"]
    base = system.lift[bin_index]
    source = system.elements[system.source_element]
    jr0 = np.asarray([base[_element_field_dofs(system, source, mode, "U")[0]]
                      for _, mode in scalar], dtype=np.complex128)
    pressure_scale = float(system.metadata["pressure_scale_pa"])
    if float(system.metadata["length_scale_m"]) <= 0 or pressure_scale <= 0:
        raise ValueError("source load model scales must be positive")
    if direction == "C":
        ratio = modulus_pa/baseline_moduli_pa["C"]
        _, _, _, jr, _ = source_local_lift(
            system, bin_index, projection, h, "C", ratio)
        delta_fj_model = (baseline_moduli_pa["F"]/pressure_scale)*(jr-jr0)
    else:
        delta_fj_model = h*(modulus_pa/pressure_scale)*(
            projection["X_scalar"]@jr0)
    rsource_model = (float(system.metadata["source_radius_m"])/
                     float(system.metadata["length_scale_m"]))
    for value, (_, mode) in zip(delta_fj_model, scalar):
        if mode.family == "R":
            dof = _element_field_dofs(system, source, mode, "U")[0]
            result[dof] -= 2.0*rsource_model*value
        else:
            udof = _element_field_dofs(system, source, mode, "U")[0]
            vdof = _element_field_dofs(system, source, mode, "V")[0]
            lam = mode.degree*(mode.degree+1.0)
            result[udof] -= 2.0*rsource_model*value
            result[vdof] += lam*rsource_model*value
    return result


def family_kind(family: str) -> str:
    return {"R": "U", "S": "U", "T": "W"}[family]


def assemble_delta_h(modes: tuple[Mode, ...], gram: np.ndarray,
                     elements: Iterable[Element], direction: str,
                     modulus_pa: float, pressure_scale_pa: float,
                     length_scale_m: float, matrix_size: int,
                     memory_cap_bytes: int = 1 << 30, *,
                     support_m: tuple[float, float] = (4.2e6, 6.2e6),
                     radial_profile=regular_weight) -> csr_matrix:
    """Assemble one unit-h Cartesian ΔH using sparse element connectivity.

    Matrices use DSpec model units. The physical modulus increment is scaled
    by pressure_scale_pa and the physical volume by length_scale_m**3.
    """
    if pressure_scale_pa <= 0 or length_scale_m <= 0 or matrix_size <= 0:
        raise ValueError("SEM SI scales must be positive")
    modes = tuple(modes)
    elements = tuple(elements)
    direction_index = ("A", "C", "F", "L", "N").index(direction)
    active = []
    for element in elements:
        radii = np.asarray(element.nodes_model, dtype=float)*length_scale_m
        lower = float(np.min(radii))
        upper = float(np.max(radii))
        tol = 8*np.finfo(float).eps*max(abs(support_m[0]), abs(support_m[1]))
        if lower < support_m[0]-tol or upper > support_m[1]+tol:
            continue
        active.append(element)
    n_q = len(active[0].nodes_model) if active else 0
    n_modes = len(modes)
    if not active:
        raise ValueError("no complete SEM element lies in the regular support")
    # Gaunt selection limits coupling to the actual L=2,M=1 and L=3,M=2
    # perturbation components. Each allowed pair has at most nq^2 entries.
    pair_count = sum(_angularly_coupled(a, b) for a in modes for b in modes)
    upper_nnz = len(active)*n_q*n_q*pair_count
    estimated_bytes = 4*upper_nnz*(4+4+16)
    if estimated_bytes > memory_cap_bytes:
        raise MemoryError(f"sparse ΔH estimate {estimated_bytes} exceeds 1 GiB cap")
    row = np.empty(upper_nnz, dtype=np.int32)
    col = np.empty(upper_nnz, dtype=np.int32)
    data = np.empty(upper_nnz, dtype=np.complex128)
    used = 0
    direction_scale = modulus_pa/pressure_scale_pa
    for element in active:
        rnodes = np.asarray(element.nodes_model, dtype=float)
        nq = len(rnodes)
        if len(element.dofs) != n_modes:
            raise ValueError("element local maps do not match harmonic modes")
        weights = (rnodes*rnodes*np.asarray(element.dr_weights_model, dtype=float)
                   * np.asarray([radial_profile(r*length_scale_m) for r in rnodes]))
        derivative = np.asarray(element.derivative_model, dtype=float)
        if derivative.shape != (nq, nq) or weights.shape != (nq,):
            raise ValueError("element radial quadrature arrays have inconsistent shapes")
        radial = (np.diag(1.0/rnodes), derivative)
        radial_blocks = tuple(
            radial[a].T @ (weights[:, None]*radial[b])
            for a in range(2) for b in range(2)
        )
        for ia in range(n_modes):
            map_a = np.asarray(element.dofs[ia], dtype=np.int64)
            if map_a.shape != (nq,):
                raise ValueError("test local map must have one global DOF per GLL node")
            for ib in range(n_modes):
                if not _angularly_coupled(modes[ia], modes[ib]):
                    continue
                map_b = np.asarray(element.dofs[ib], dtype=np.int64)
                if map_b.shape != (nq,):
                    raise ValueError("trial local map must have one global DOF per GLL node")
                localgram = gram[ia, ib, direction_index]
                block = sum(localgram[a, b]*radial_blocks[2*a+b]
                            for a in range(2) for b in range(2))
                block *= direction_scale
                nz_i, nz_j = np.nonzero(block)
                count = len(nz_i)
                if count:
                    row[used:used+count] = map_a[nz_i]
                    col[used:used+count] = map_b[nz_j]
                    data[used:used+count] = block[nz_i, nz_j]
                    used += count
    if used and max(int(row[:used].max()), int(col[:used].max())) >= matrix_size:
        raise ValueError("local-to-global map exceeds declared matrix size")
    matrix = coo_matrix((data[:used], (row[:used], col[:used])),
                        shape=(matrix_size, matrix_size)).tocsr()
    matrix.sum_duplicates()
    return matrix


if __name__ == "__main__":
    raise SystemExit(main())

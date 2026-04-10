from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

from .plot_common import loadtxt_or_exit


@dataclass
class DiffSummary:
    average_percent: float
    max_percent: float


@dataclass
class Ex1Dataset:
    frequency: np.ndarray
    yspec_z: np.ndarray
    specnm_z: np.ndarray
    dspecm_z: np.ndarray
    norm_z: float
    specnm_column: int


@dataclass
class Ex1DiffReport:
    fmin: float
    fmax: float
    specnm_column: int
    yspec_vs_dspecm: DiffSummary
    specnm_vs_dspecm: DiffSummary


class SpectraComparison:
    """Utility for loading spectra tables and computing relative-difference metrics."""

    def __init__(self, path, delimiter=";"):
        self.path = path
        self.data = loadtxt_or_exit(path, delimiter=delimiter)

    def window(self, lowidx, upidx):
        return self.data[lowidx:upidx, :]

    def frequency_window_indices(self, fmin, fmax):
        idx_min = np.searchsorted(self.data[:, 0], fmin)
        idx_max = np.searchsorted(self.data[:, 0], fmax, side="right")
        return idx_min, idx_max

    def ex1_dataset(self, fmin, fmax, yspec_col=12, dspecm_col=3, specnm_preferred_col=21, norm_scale=1.01):
        idx_min, idx_max = self.frequency_window_indices(fmin, fmax)
        window = self.data[idx_min:idx_max, :]

        specnm_col = specnm_preferred_col if self.data.shape[1] > specnm_preferred_col else yspec_col
        yspec_z = window[:, yspec_col]
        dspecm_z = window[:, dspecm_col]
        specnm_z = window[:, specnm_col]
        norm_z = float(np.max(np.abs(yspec_z)) * norm_scale)

        return Ex1Dataset(
            frequency=window[:, 0],
            yspec_z=yspec_z,
            specnm_z=specnm_z,
            dspecm_z=dspecm_z,
            norm_z=norm_z,
            specnm_column=specnm_col,
        )

    def ex1_diff_report(self, fmin, fmax, yspec_col=12, dspecm_col=3, specnm_preferred_col=21, norm_scale=1.01):
        dataset = self.ex1_dataset(
            fmin=fmin,
            fmax=fmax,
            yspec_col=yspec_col,
            dspecm_col=dspecm_col,
            specnm_preferred_col=specnm_preferred_col,
            norm_scale=norm_scale,
        )

        yspec_vs_dspecm = self.summarize_difference(dataset.yspec_z, dataset.dspecm_z, dataset.norm_z)
        specnm_vs_dspecm = self.summarize_difference(dataset.dspecm_z, dataset.specnm_z, dataset.norm_z)

        return Ex1DiffReport(
            fmin=fmin,
            fmax=fmax,
            specnm_column=dataset.specnm_column,
            yspec_vs_dspecm=yspec_vs_dspecm,
            specnm_vs_dspecm=specnm_vs_dspecm,
        )

    @staticmethod
    def normalization(*arrays, eps=1e-12):
        return max(np.max(np.abs(arr)) for arr in arrays) + eps

    @staticmethod
    def relative_difference_percent(a, b, norm):
        return np.abs(a - b) / norm * 100

    @classmethod
    def summarize_difference(cls, a, b, norm):
        diff = cls.relative_difference_percent(a, b, norm)
        return DiffSummary(
            average_percent=float(np.mean(diff)),
            max_percent=float(np.max(diff)),
        )

    @classmethod
    def summarize_pairs(cls, pairs: Dict[str, Tuple[np.ndarray, np.ndarray]], norm):
        return {name: cls.summarize_difference(a, b, norm) for name, (a, b) in pairs.items()}

    @staticmethod
    def l2_relative(reference, target):
        ref_norm = np.sqrt(np.sum(reference ** 2))
        return np.sqrt(np.sum((target - reference) ** 2)) / ref_norm

    @classmethod
    def l2_relative_many(cls, reference, series: Dict[str, np.ndarray]):
        return {name: cls.l2_relative(reference, values) for name, values in series.items()}

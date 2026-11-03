"""What is computed: the degrees, the band, the time series and the physics,
as the members of an input file give them."""

from dataclasses import dataclass

from . import _core

_GRAVITY = {"none": _core.Gravity.none, "cowling": _core.Gravity.cowling,
            "full": _core.Gravity.full}
_QUANTITY = {"displacement": _core.Quantity.displacement,
             "velocity": _core.Quantity.velocity,
             "acceleration": _core.Quantity.acceleration}


def _choice(name, value, choices):
    try:
        return choices[value]
    except KeyError:
        raise ValueError(f"{name} should be one of {', '.join(choices)}, "
                         f"not {value!r}") from None


@dataclass
class Calculation:
    """What is computed.

    Parameters
    ----------
    degrees
        The highest spherical harmonic degree summed, or the lowest and
        highest. Enough are needed for the fundamental surface waves at the
        highest frequency, about 10 for each mHz, even where only body waves
        are wanted.
    frequencies_mhz
        The lowest and highest frequencies computed (mHz).
    length_min
        The length of the seismograms (minutes).
    time_step_s
        The interval between samples (s); by default a quarter of the shortest
        period.
    filter_mhz
        The four corners (mHz) of a raised-cosine taper of the spectra: rising
        from zero to one between the first two, falling back between the last
        two. It has zero phase. None for no taper.
    gravity
        "none", "cowling" (the equilibrium gravity but no perturbation of the
        potential) or "full".
    attenuation
        Whether the model attenuates, with the dispersion that goes with it.
    quantity
        "displacement", "velocity" or "acceleration".
    corrections
        For acceleration, the change in gravity seen by an instrument that
        moves with the ground, and the tilt of its axis.
    tolerance
        The relative tolerance of the radial integrations.
    """

    degrees: int | tuple[int, int]
    frequencies_mhz: tuple[float, float]
    length_min: float
    time_step_s: float | None = None
    filter_mhz: tuple[float, float, float, float] | None = None
    gravity: str = "full"
    attenuation: bool = True
    quantity: str = "displacement"
    corrections: bool = False
    tolerance: float = 1e-5

    def _core(self):
        """The calculation for the library, in SI units, checked as an input
        file's is."""
        c = _core.Calculation()
        low, high = ((0, self.degrees) if isinstance(self.degrees, int)
                     else self.degrees)
        if not 0 <= low <= high:
            raise ValueError("degrees should satisfy 0 <= min <= max")
        c.degrees = _core.Degrees(int(low), int(high))
        f_min, f_max = self.frequencies_mhz
        if not 0 < f_min < f_max:
            raise ValueError("frequencies_mhz should satisfy 0 < min < max")
        c.band = _core.FrequencyBand(f_min * 1e-3, f_max * 1e-3)
        time_step = (0.25 / (f_max * 1e-3) if self.time_step_s is None
                     else self.time_step_s)
        if not (self.length_min > 0 and time_step > 0):
            raise ValueError("length_min and time_step_s should be positive")
        c.sampling = _core.Sampling(self.length_min * 60, time_step)
        if self.filter_mhz is not None:
            corners = [f * 1e-3 for f in self.filter_mhz]
            if len(corners) != 4 or not (0 <= corners[0] <= corners[1]
                                         < corners[2] <= corners[3]):
                raise ValueError(
                    "filter_mhz should be four corners in increasing order, "
                    "with a pass band between the second and third")
            c.taper = _core.Taper(*corners)
        c.gravity = _choice("gravity", self.gravity, _GRAVITY)
        c.attenuation = bool(self.attenuation)
        c.quantity = _choice("quantity", self.quantity, _QUANTITY)
        c.corrections = bool(self.corrections)
        if not 0 < self.tolerance < 1:
            raise ValueError("tolerance should lie between 0 and 1")
        c.solver = _core.SolverOptions(tolerance=self.tolerance)
        return c

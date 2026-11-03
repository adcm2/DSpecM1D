"""The earth model, read from a deck in MINEOS's format."""

import os
from dataclasses import dataclass
from importlib import resources

import numpy as np

from . import _core

# The decks that come with the package, by name.
DECKS = {
    "prem": "prem.200",  # PREM with its ocean
    "prem-no-ocean": "prem.200.no",  # PREM with the ocean taken off
}


@dataclass(frozen=True)
class Layer:
    """A layer of a model, between two discontinuities (km)."""

    index: int
    lower_km: float
    upper_km: float
    fluid: bool


class Model:
    """An earth model: layers of solid and fluid, transversely isotropic and
    anelastic, with the gravity that follows from their density.

        model = yspec.Model("prem.200")
        model = yspec.Model.prem()
    """

    def __init__(self, path):
        self.path = os.fspath(path)
        self._model = _core.Model.read_deck(self.path)

    @classmethod
    def from_text(cls, text):
        """A model from the text of a deck."""
        model = cls.__new__(cls)
        model.path = None
        model._model = _core.Model.read_deck_text(text)
        return model

    @classmethod
    def prem(cls, name="prem"):
        """PREM, from a deck that comes with the package: "prem", with its
        ocean, or "prem-no-ocean"."""
        try:
            deck = DECKS[name]
        except KeyError:
            raise ValueError(f"there is no deck {name!r}: there are "
                             f"{', '.join(DECKS)}") from None
        with resources.as_file(resources.files("yspec") / "models" /
                               deck) as path:
            return cls(path)

    def __repr__(self):
        return (f"Model({self.path!r}, {len(self.layers)} layers, radius "
                f"{self.surface_radius_km:g} km)")

    @property
    def surface_radius_km(self):
        return self._model.surface_radius / 1e3

    @property
    def layers(self):
        """The layers, from the centre up."""
        return [Layer(i, *(r / 1e3 for r in self._model.layer_bounds(i)),
                      self._model.is_fluid(i))
                for i in range(self._model.number_of_layers)]

    def profile(self, radius_km, side="below"):
        """The model at radii (km), in SI units: a dict of arrays of the
        density (rho), the Love moduli (A, C, F, L, N), the Voigt bulk and
        shear moduli (kappa, mu), the isotropic speeds from them (vp, vs), the
        inverse quality factors (q_kappa, q_mu) and gravity (g). At a
        discontinuity, `side` ("below" or "above") says which side."""
        sides = {"below": _core.Side.below, "above": _core.Side.above}
        if side not in sides:
            raise ValueError('side should be "below" or "above"')
        radii = np.atleast_1d(np.asarray(radius_km, float)) * 1e3
        rows = [self._model.material(
                    self._model.layer_index(r, sides[side]), r)
                for r in radii]
        values = {key: np.array([row[key] for row in rows]) for key in rows[0]}
        values["vp"] = np.sqrt((values["kappa"] + 4 * values["mu"] / 3) /
                               values["rho"])
        values["vs"] = np.sqrt(values["mu"] / values["rho"])
        return values

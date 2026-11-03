"""Compact labels for display clipping on spatial panels."""


def clipping_stats(limit, maximum, fraction):
    """Three short lines that fit beside narrow section/map colorbars."""
    return (f"limit ±{limit:.3g}\n"
            f"slice max {maximum:.3g}\n"
            f"saturated {fraction:.2f}%")

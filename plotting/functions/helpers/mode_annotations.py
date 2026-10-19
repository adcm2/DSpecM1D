from pathlib import Path

import numpy as np


def _read_mode_file(filepath):
    path = Path(filepath).expanduser()
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"Cannot read mode catalogue '{path}': {exc}") from exc

    header_idx = next(
        (
            idx
            for idx, line in enumerate(lines)
            if "mode" in line.lower()
            and "phs vel" in line.lower()
            and "w(mhz)" in line.lower()
        ),
        None,
    )
    if header_idx is None:
        raise ValueError(
            f"Mode catalogue '{path}' has no MINEOS table header "
            "('mode', 'phs vel', and 'w(mhz)')."
        )

    modes = []
    for line_number, line in enumerate(lines[header_idx + 2 :], start=header_idx + 3):
        if not line.strip():
            continue
        parts = line.split()
        try:
            radial_order = int(parts[0])
        except ValueError:
            # MINEOS data rows start with an integer radial order. Textual
            # footers and section labels are outside the table body.
            continue

        if len(parts) < 5:
            raise ValueError(
                f"Malformed mode row in '{path}' at line {line_number}: "
                "expected at least 5 columns (mode orders and frequency)."
            )
        try:
            angular_order = int(parts[2])
            frequency = float(parts[4])
        except ValueError as exc:
            raise ValueError(
                f"Malformed mode row in '{path}' at line {line_number}: "
                "expected integer mode orders and a numeric frequency in column 5."
            ) from exc
        if radial_order < 0 or angular_order < 0:
            raise ValueError(
                f"Malformed mode row in '{path}' at line {line_number}: "
                "radial and angular mode orders must be nonnegative."
            )
        if len(parts[1]) != 1 or not parts[1].isalpha():
            raise ValueError(
                f"Malformed mode row in '{path}' at line {line_number}: "
                "mode class must be one alphabetic character."
            )
        if not np.isfinite(frequency):
            raise ValueError(
                f"Malformed mode row in '{path}' at line {line_number}: "
                "frequency must be finite."
            )
        type_mode = parts[1].upper()
        modes.append((parts[0], type_mode, parts[2], frequency))

    return modes


def get_modes_in_range(fmin, fmax, mode_files=None):
    """Return sorted MINEOS mode labels; no catalogues means no annotations."""
    modes = []
    for filepath in mode_files or []:
        for radial_order, type_mode, angular_order, frequency in _read_mode_file(filepath):
            if fmin <= frequency <= fmax:
                name = f"${{_{radial_order}}}${type_mode}${{_{{{angular_order}}}}}$"
                modes.append((name, frequency))

    modes.sort(key=lambda item: item[1])
    return [mode[0] for mode in modes], [mode[1] for mode in modes]


def compute_mode_heights(xvalues, xgrid, y_arrays, window=3, default=0.05):
    heights = []
    for x in xvalues:
        idx = np.searchsorted(xgrid, x)
        if idx >= len(xgrid):
            heights.append(default)
            continue

        start_w = max(0, idx - window)
        end_w = min(len(xgrid), idx + window + 1)
        local_max = max(np.max(y[start_w:end_w]) if len(y[start_w:end_w]) > 0 else 0 for y in y_arrays)
        heights.append(local_max)

    return heights


def annotate_modes(ax, names, xvalues, yvalues, ytop, dy, line_width, alternate=False):
    for i, name in enumerate(names):
        ax.axvline(
            xvalues[i],
            ymin=0,
            ymax=(yvalues[i] + dy) / ytop,
            color="black",
            linestyle="--",
            linewidth=line_width,
        )
        text_y = yvalues[i] + dy
        if alternate and (i % 2 == 1):
            text_y += dy
        ax.text(xvalues[i], text_y, name, fontsize=14, fontweight="bold", ha="center")

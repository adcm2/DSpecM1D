import numpy as np

DEFAULT_MODE_FILES = [
    "/home/adcm2/space/c++/mineos/DEMO/MYEX/lf_prem_R",
    "/home/adcm2/space/c++/mineos/DEMO/MYEX/lf_prem_S",
]


def get_modes_in_range(fmin, fmax, mode_files=None):
    files = mode_files or DEFAULT_MODE_FILES
    modes = []
    for filepath in files:
        try:
            with open(filepath, "r") as f:
                lines = f.readlines()
        except FileNotFoundError:
            continue

        start_idx = 0
        for i, line in enumerate(lines):
            if "mode" in line and "phs vel" in line and "w(mhz)" in line:
                start_idx = i + 2
                break

        if start_idx == 0:
            continue

        for line in lines[start_idx:]:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                n = parts[0]
                type_mode = parts[1].upper()
                l = parts[2]
                freq = float(parts[4])
            except ValueError:
                continue
            if fmin <= freq <= fmax:
                name = f"${{_{n}}}${type_mode}${{_{{{l}}}}}$"
                modes.append((name, freq))

    modes.sort(key=lambda x: x[1])
    return [m[0] for m in modes], [m[1] for m in modes]


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

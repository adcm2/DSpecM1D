import argparse
from pathlib import Path

from helpers.plot_paths import default_build_dir


# User-editable plotting configuration
FREQUENCY_RANGE = (0.0, 5.0)
ZOOM_REGIONS = [(0.28, 0.50), (2.27, 2.55)]


def _arguments():
    parser = argparse.ArgumentParser(
        description="Plot the ex1 DSpecM1D, YSpec and specnm comparison."
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=default_build_dir(),
        help="build tree containing output/ex1_w.out (default: build/release)",
    )
    parser.add_argument(
        "--data",
        type=Path,
        help="comparison file; overrides --build-dir/output/ex1_w.out",
    )
    parser.add_argument(
        "--mode-file",
        type=Path,
        action="append",
        default=[],
        metavar="PATH",
        help="MINEOS mode catalogue to annotate (repeat to add catalogues)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="render without opening a window (use --save to retain the figure)",
    )
    parser.add_argument(
        "--save",
        type=Path,
        metavar="PATH",
        help="save the figure to PATH",
    )
    return parser


def main(argv=None):
    parser = _arguments()
    args = parser.parse_args(argv)
    data_file = args.data if args.data is not None else args.build_dir / "output" / "ex1_w.out"
    data_file = data_file.expanduser().resolve()
    mode_files = [path.expanduser().resolve() for path in args.mode_file]

    if args.no_show:
        import matplotlib

        matplotlib.use("Agg")

    from helpers.comparison_metrics import SpectraComparison
    from helpers.ex1_workflow import output_to_screen_code_diff, plot_ex1_comparison

    try:
        comparison = SpectraComparison(str(data_file))
        fmin, fmax = FREQUENCY_RANGE

        output_to_screen_code_diff(
            comparison=comparison,
            fmin=fmin,
            fmax=fmax,
            zoom_regions=ZOOM_REGIONS,
        )

        fig, _ = plot_ex1_comparison(
            comparison=comparison,
            fmin=fmin,
            fmax=fmax,
            zoom_regions=ZOOM_REGIONS,
            show=not args.no_show,
            mode_files=mode_files,
        )
    except ValueError as exc:
        parser.error(str(exc))

    if args.save is not None:
        save_path = args.save.expanduser().resolve()
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Saved plot: {save_path}")


if __name__ == "__main__":
    main()

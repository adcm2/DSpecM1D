#%%
from helpers.comparison_metrics import SpectraComparison
from helpers.ex1_workflow import output_to_screen_code_diff, plot_ex1_comparison

# User-editable configuration
DATA_FILE = "../outputs/ex1_w.out"
FREQUENCY_RANGE = (0.0, 5.0)
ZOOM_REGIONS = [(0.28, 0.50), (2.27, 2.55)]


def main():
    comparison = SpectraComparison(DATA_FILE)
    fmin, fmax = FREQUENCY_RANGE

    output_to_screen_code_diff(
        comparison=comparison,
        fmin=fmin,
        fmax=fmax,
        zoom_regions=ZOOM_REGIONS,
    )

    plot_ex1_comparison(
        comparison=comparison,
        fmin=fmin,
        fmax=fmax,
        zoom_regions=ZOOM_REGIONS,
        show=True,
    )


if __name__ == "__main__":
    main()
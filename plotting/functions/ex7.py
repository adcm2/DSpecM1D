#%%
from helpers.comparison_metrics import SpectraComparison
from helpers.ex7_workflow import output_to_screen_code_diff_ex7, plot_ex7_comparison

# User-editable configuration
DATA_FILE = "../outputs/ex7_t.out"
LOW_INDEX = 400
UP_INDEX = 1301
TIME_UPPER_LIMIT = None


def main():
    comparison = SpectraComparison(DATA_FILE)
    output_to_screen_code_diff_ex7(
        comparison=comparison,
        lowidx=LOW_INDEX,
        upidx=UP_INDEX,
        tup=TIME_UPPER_LIMIT,
    )
    plot_ex7_comparison(
        comparison=comparison,
        lowidx=LOW_INDEX,
        upidx=UP_INDEX,
        tup=TIME_UPPER_LIMIT,
        show=True,
    )


if __name__ == "__main__":
    main()

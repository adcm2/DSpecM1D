#%%
from helpers.comparison_metrics import SpectraComparison
from helpers.ex2_workflow import output_to_screen_code_diff_ex2, plot_ex2_comparison

# User-editable configuration
DATA_FILE = "../outputs/ex2_t.out"
LOW_INDEX = 0
UP_INDEX = 18001
TIME_UPPER_LIMIT = 18000


def main():
    comparison = SpectraComparison(DATA_FILE)
    output_to_screen_code_diff_ex2(
        comparison=comparison,
        lowidx=LOW_INDEX,
        upidx=UP_INDEX,
        tup=TIME_UPPER_LIMIT,
    )
    plot_ex2_comparison(
        comparison=comparison,
        lowidx=LOW_INDEX,
        upidx=UP_INDEX,
        tup=TIME_UPPER_LIMIT,
        show=True,
    )


if __name__ == "__main__":
    main()

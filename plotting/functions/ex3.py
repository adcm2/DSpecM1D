#%%
from helpers.comparison_metrics import SpectraComparison
from helpers.ex3_workflow import output_to_screen_code_diff_ex3, plot_ex3_comparison

# User-editable configuration
DATA_FILE = "../outputs/ex3_t.out"
LOW_INDEX = 400
UP_INDEX = 1301


def main():
    comparison = SpectraComparison(DATA_FILE)
    output_to_screen_code_diff_ex3(
        comparison=comparison,
        lowidx=LOW_INDEX,
        upidx=UP_INDEX,
    )
    plot_ex3_comparison(
        comparison=comparison,
        lowidx=LOW_INDEX,
        upidx=UP_INDEX,
        show=True,
    )


if __name__ == "__main__":
    main()

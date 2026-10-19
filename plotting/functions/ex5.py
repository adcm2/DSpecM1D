#%%
from helpers.plot_paths import build_output_path
from helpers.ex5_workflow import output_to_screen_code_diff_ex5, plot_ex5

# User-editable configuration
NQ5_FILE = build_output_path("ex5_NQ5_step_error_50.out")
NQ6_FILE = build_output_path("ex5_NQ6_step_error_50.out")
NQ4_FILE = build_output_path("ex5_NQ4_step_error_50.out")


def main():
    output_to_screen_code_diff_ex5(
        path_nq5=NQ5_FILE,
        path_nq6=NQ6_FILE,
        path_nq4=NQ4_FILE,
    )
    plot_ex5(
        path_nq5=NQ5_FILE,
        path_nq6=NQ6_FILE,
        path_nq4=NQ4_FILE,
        show=True,
    )


if __name__ == "__main__":
    main()
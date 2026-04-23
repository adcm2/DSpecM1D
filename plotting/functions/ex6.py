#%%
from helpers.ex6_workflow import output_to_screen_code_diff_ex6, plot_ex6

# User-editable configuration
FILE_MAIN = "../outputs/ex6_radial_response_0.010000.out"
FILE_MAIN_2 = "../outputs/ex6_radial_response_0.001000.out"
FILE_MAIN_3 = "../outputs/ex6_radial_response_0.000100.out"
FILE_FREQS = "../outputs/ex6_w_0.010000.out"
FILE_N2 = "../outputs/ex6_N2_0.010000.out"
MAX_FREQ_PANELS = 6


def main():
    output_to_screen_code_diff_ex6(
        path_mf=FILE_MAIN,
        path_mf2=FILE_MAIN_2,
        path_mf3=FILE_MAIN_3,
        path_freqs=FILE_FREQS,
        path_n2=FILE_N2,
        max_freq_panels=MAX_FREQ_PANELS,
    )
    print("Hello just checking the output above for the workflow context. Now plotting...")
    plot_ex6(
        path_mf=FILE_MAIN,
        path_mf2=FILE_MAIN_2,
        path_mf3=FILE_MAIN_3,
        path_freqs=FILE_FREQS,
        path_n2=FILE_N2,
        max_freq_panels=MAX_FREQ_PANELS,
        show=True,
    )


if __name__ == "__main__":
    main()
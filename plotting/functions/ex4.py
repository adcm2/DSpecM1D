#%%
from helpers.ex4_workflow import output_to_screen_code_diff_ex4, plot_ex4

# User-editable configuration
WAVEFORM_FILE = "../outputs/ex4.out"
TRAVEL_TIME_FILE = "../outputs/travel_times.txt"


def main():
    output_to_screen_code_diff_ex4(
        path_waveforms=WAVEFORM_FILE,
        path_travel_times=TRAVEL_TIME_FILE,
    )
    plot_ex4(
        path_waveforms=WAVEFORM_FILE,
        path_travel_times=TRAVEL_TIME_FILE,
        show=True,
    )


if __name__ == "__main__":
    main()
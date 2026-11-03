# Fresh delta review: post-analysis script extraction

**Result: PASS for the bounded packaging change.** This review compares the current plan with the preserved pre-extraction plan and checks the extracted scripts and argument mapping. No queue was launched, and no numerical run, build, or SSH was performed.

The wait and aggregation inline programs were moved byte-for-byte into `tools/check_wait_for_terminal.py` and `tools/check_post_analysis.py`. Both plan commands now invoke those staged files, and each script is declared as an input. The wait script receives the same dependency-task directory and `{JOB_DIR}/dependency.json` output path. The aggregate script retains the same seven positional arguments and ordering, including tool directory, output directory, roots JSON, global root, layout root, prior summary, and repeat root. For both affected jobs, only `command` and `inputs` differ from the preserved plan; all other job fields are unchanged. The remaining eight jobs are byte-equivalent, and global plan metadata, ten-job ordering, four-thread cap, and environment are unchanged.

Both extracted scripts pass `py_compile`. Every plan job has nonempty inputs. The aggregate argv’s largest token is now 481 bytes, rather than carrying the multi-kilobyte program inline, addressing the reported supervisor file-probe failure without changing its logic or arguments.

## SHA-256

- `configs/post-analysis-plan.json`: `72fc84d5ca243e316c04e6463249b0bd787c7365991bf523044b1be80061af00`
- `review/post-analysis-before-script-extraction.json`: `d54a68ed32389fa6430f9b43b23ce892df7793e9b2b4557e35601ed4e3ed53b2`
- `tools/check_wait_for_terminal.py`: `51fe4b4ee38825618519880a21884f85045238d179922ae6a613bbe18751508b`
- `tools/check_post_analysis.py`: `494bf448495cc8bcc609c7d7e051d939a4d098ca2815417e2f4918831bf9d367`

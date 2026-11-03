# Direct-reference global-kernel wrapper status

## Change log

- Added `check_global_original.py` as a thin wrapper over the unchanged reviewed `check_global_pair.assess()` direct-full-grid comparison. It retains the original physical, native-identity, frequency-audit, q/energy, and `2e-13` channel gates and returns the full helper report.
- Adapts only the final producer's `serial_grid_and_context_setup_wall_s` to the reviewed validator's `serial_context_setup_wall_s` name using the shared in-memory adapter. It asserts every original timing field remains unchanged, restores the reviewed validator in `finally`, and reports adaptation details and the reviewed helper SHA-256.
- Added a completed-receipt end-to-end fixture with 179 direct reference sides, eight selected global sides, the real reviewed comparison path, and the production timing-key adapter.

## Validation

- `python3 -m py_compile check_global_original.py test_check_global_original.py`: passed.
- `python3 -m unittest -v test_check_global_original.py`: 1 end-to-end test passed.
- No build, SSH, or scientific execution was performed.

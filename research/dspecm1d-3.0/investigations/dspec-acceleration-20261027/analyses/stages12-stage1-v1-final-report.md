# DSpecM1D Stage 1 final report

Status: all Stage 1 work completed. All 12 controls and all 10 measurement jobs passed; inputs were unchanged and cleanup was verified. No Stage 2 code has been implemented yet.

Accepted source: `/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/acceleration-inputs/stages12-stage1-v1/candidates/stages12-stage1-profile-v1`. The exact local source is [`investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v1-frozen-inputs/candidate`](investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v1-frozen-inputs/candidate); archive SHA256 `52fccc457361e1fc8b1e865f05738f93a04e933f2725ac24169f59770de3ad63`; tree SHA256 `17609c812f316d91bc22ca058460bc6c19c483e9dbfd8d666000dc746bda5daf`.

## Native L200 repeats

| Run | result.json elapsed_s | supervisor elapsed_s | kernel_s | peak observed RSS GiB | threads |
|---|---:|---:|---:|---:|---:|
| N1-native-L200 | 51.840044503 | 52.480971456 | 40.678408821 | 2.4242 | 9 |
| N2-native-L200 | 51.832922741 | 52.394791400 | 40.663871973 | 2.4272 | 9 |
| N3-native-L200 | 52.140116475 | 52.730584407 | 40.896136167 | 2.4237 | 9 |

Native result elapsed min/median/max: `51.832922741 / 51.840044503 / 52.140116475 s`. Maximum observed native RSS was `2.4272 GiB` and peak thread count was `9`.

P1 profile accounting closed exactly at 52.164586220 s with 150 rows and 100 counters. P1 parity passed to N1 and frozen v9; N2/N3 passed byte parity to N1. The L3 finite-difference relative error was `6.440263106e-09`. The P1 profile packet is [`receipts/stages12-stage1-v3-profile-packet/`](investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v3-profile-packet). Raw rows, counters, and the scope-name correction map are saved as CSV/JSON under `analyses/`.

The archived accepted source predates the optional local cache-counter guard; the guard draft was not built or measured. Two prelaunch path assertions (controls v1/v2) are preserved and no jobs ran in those failed attempts. V3 plans passed all gates. Runtime component label `output_io` remains misclassified for some scope names; use the explicit scope map and raw scope names.

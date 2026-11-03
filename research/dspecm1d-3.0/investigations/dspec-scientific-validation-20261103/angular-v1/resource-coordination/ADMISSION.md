# DSpecM1D angular allocation: resource admission

3 November 2026. ADMITTED after independent review, owner review and live registration probe. This is resource admission only; the DSpec owner retains scientific/source/build/result gates.

Receipt `deployment_receipt.json` SHA256 dc59cc0f82a84409e5338c2f7c73f37a683a82e4aaf6b060e2c0ce35f38503e1. The remote copy has the same hash. The earlier conditional proposal remains preserved as historical evidence.

## Allocation

- Separate `dspecscienceangular20261103.slice`: 40,000,000,000 bytes, 12 tasks, swap 0, CPU quota 800%.
- One `dspec-science-angular-v1-build.service` (32 GiB), OR `dspec-science-angular-v1-run.service` (16 GiB), OR `dspec-science-angular-v1-check.service` (8 GiB). All have TasksMax=12 and MemorySwapMax=0. Never overlap these services; all descendants must remain in the service cgroup.
- Eight outer workers, one inner library thread. Set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, BLIS_NUM_THREADS, VECLIB_MAXIMUM_THREADS and NUMEXPR_NUM_THREADS to 1.
- Existing CCDSM shared 299 GB/90-task slice and its 260 GB/67-task numerical child allocation are unchanged. Pilot/M3D retains a distinct 8 GB/12-task reservation. Combined reservations are 318 GB; combined separate slice caps 339 GB, below the user's 600 GB owned-work ceiling.

## Before every service launch

Read fresh state at `/space/adcm2/research/runs/CCDSM/extended_elastic_20261031_v1/sp12_20261101_v1/resource_coordination_20261101_v1/global_ram_guard_v4/state.json`. Require age between 0 and 5 seconds, schema `tunya-global-memory-observer-v4`, valid active metric, MONITORING, armed, no paused, frozen or unknown owned units, and active guardian service. Verify live guardian source SHA256 `03421ab7cb49f3e28cfa08c3907ff264cae34d3a4b6f9ea170120b1976989d51`.

Require active RAM + 40 GB + all other unconsumed reservations < 825 GB. A conservative sufficient test is active RAM + 318 GB < 825 GB; this overcounts memory already in use. Never substitute cache-inclusive occupancy. Require the exact separate slice caps above, no other live DSpecM1D service, and sufficient task/memory headroom within that slice. Do not adjust existing shared caps or borrow Pilot's reservation.

Launch through `systemd-run --user --slice=dspecscienceangular20261103.slice`, exact service name above, explicit child `MemoryMax`, `MemorySwapMax=0`, `TasksMax=12` and `CPUQuota=800%`, working directory and thread environment. Pass the already-reviewed build/campaign/check command as the service's command, retaining its own immutable source and output contracts. A bare nohup command outside the kernel-owned service is insufficient. Do not relaunch an existing scientific state/output.

Immediately verify MainPID and start time, actual ControlGroup under the registered slice, readable `cgroup.freeze=0`/FreezerState=running, memory and task counters, child/slice limits and zero kernel max/oom/kill events. Inspect the deployed guardian's `owned_units()` from an app.slice service (read-only import, do not invoke main) and confirm the exact live numerical/build service appears with readable freezer/accounting. Observe CPU/log advancement. Continue fresh-state/resource monitoring throughout the service.

## Guardian operation

The SAME `ccdsm-node-memory-guard-v4.service` in app.slice now inventories the original research slice/prefixes and exactly the new sibling slice/prefix. Source and state schema/metric/threshold/pause logic otherwise remain unchanged. No second freezing guardian. Pause thresholds remain 825/875 GB, whole-node active boundary 900 GB. Pause at most the largest ONE owned service; never kill, freeze a whole slice, automatically thaw or freeze a second service. If a guard/metric/freezer/conflict issue occurs, preserve processes and request coordinated operator action.

The guardian restart began a new observed-peak/action epoch. Prior source, unit configuration, state, peak and actions are preserved in the remote registration directory. No historical peak was recalculated. A 64 MB/2-task temporary probe demonstrated exact new-root inventory, readable freezer and kernel accounting. The same SP12 native PID/starttime stayed running with CPU advancement during deployment.

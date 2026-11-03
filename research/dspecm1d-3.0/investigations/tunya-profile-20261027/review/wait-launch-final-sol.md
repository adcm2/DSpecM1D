# Short final Sol review: detached prerequisite waiter

27 October 2026. Read-only static review after `wait-launch-v2-fresh-luna.md` PASS. No helper execution, mock tests, builds, scientific runs or SSH performed.

**PASS for the reviewed helper under the documented frozen, verified input-tree assumption.** It can wait for prerequisites and hand off to the already reviewed staged numerical supervisor. This is orchestration acceptance, not evidence that transfer/build/numerical gates succeeded.

Reviewed hashes:

```text
d83e3689d6c2fa48b26295bf3ac014fd2a84c633543c6a46ebbe14fb456f2fe2  tools/wait_and_launch.py
3af8a351c66aadd1bcc86a2cc04e303f0cdde27a4aba4b634cede871139901a7  tools/qualification_checks.py
e9a0e8eb947895c435c87697f9bcafca159b789b8b4614286bae4bb8080b9d7a  configs-v3/qualification-and-profile-plan.json
```

The helper waits at most four hours for completed build-03 state and a verified full-transfer marker. A partial JSON transfer write remains pending; readable non-success transfer status and unexpected/failed terminal build states fail closed. It rechecks both statuses after waiting. The helper only reads prerequisites and prints status until it replaces itself with the reviewed supervisor; it does not write/delete task evidence or launch individual scientific executables itself.

Before handoff it pins the supplied checker/plan SHA-256 values, permits only `profile-01`, and refuses an existing task directory. The receipt must be completed, its SHA must match `build-result.json` in the completed build-03 supervisor job, and that job must report unchanged inputs. The receipt's build root must match the staged reviewed build config. Exactly seven expected variants are required. Each command receipt must report return code zero, unchanged tracked inputs and the same binary digest as the summary; each current binary must hash to that digest and resolve to its exact configured output path. The hash-pinned plan must use every verified executable, and unverified build-root files in its command/input paths are rejected.

It then uses `exec` to run the reviewed supervisor with the pinned plan path and task ID. The supervisor retains campaign locking, duplicate-task refusal, serial single-thread numerical jobs, the maximum 50-thread userspace watchdog, evidence preservation and the receiver/full-L3/checker-before-high-degree sequence. The checker can stop progression to the 64-angle diagnostics. Parent should preserve the exclusive detached log and verify the detached session, as already planned.

The same residual boundary as fresh Luna applies: checker/plan are hashed during verification, then passed by path; the supervisor and build config are trusted as part of the verified transferred tree. They must remain frozen through handoff. This PASS does not support concurrent modification of that staged tree. No additional implementation change is required under that assumption. Parent-reported build completion/timing/thread counts and upload progress were not independently verified in this static review.

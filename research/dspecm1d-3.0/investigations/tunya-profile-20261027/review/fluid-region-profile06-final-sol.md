# Final short Sol review: descriptive fluid region and profile-06 retry

27 October 2026. Read-only static review after `fluid-region-profile06-fresh-luna.md` PASS. No implementation edits, scientific runs, builds or SSH performed. The parent reports profile-05 stopped in assessment before any high-degree job and the saved radial arrays have zero fluid nodes in the old shell and 22 in the replacement shell; those remote counts/results were not independently recomputed here.

**PASS for this narrow versioned correction and assessment-first profile-06 plan.** This preserves the earlier checker/adapter assessment contract and does not weaken numerical acceptance.

Verified SHA-256:

```text
87605441cb78a938f0461c906b8592092f3bd75e94535eea4bbff4ffe106d7c6  tools-v5/qualification_checks.py
a5fd84a8209e69f657eacb62f6cccde98a4cfd3ec2a68a8e511088a0c220e446  tools-v5/qualification_checks_control.py
1ecfc12aef8c146515082844564f652392e49fb57d868fc91e961e0f5c22c834  configs-v8/qualified-component-profile-plan.json
```

The checker differs from the original only in the descriptive cross-code fluid-shell selection: 4.2–6.2 Mm becomes 2.2–3.2 Mm, with corresponding comment/error/metric-label updates. Both models' fluid flags remain required. Every angular sample at selected radial sides is retained. The existing kappa-versus-A+C+F metric remains `descriptive_only_no_global_accuracy_claim`; there is no numerical accuracy threshold attached to that comparison. The correction removes a descriptive empty-selection failure when the reported saved-node evidence holds. All parity, control, full-grid, radial, energy, q-dot-P, finite-value and provenance checks remain unchanged.

The adapter differs from the reviewed version only by importing `qualification_checks.py` from its own directory. Deployment at I8/tools-v5 for **both** Python files therefore matches the actual import. Existing before/after provenance binding covers the sibling checker used. Original v4/v6 failures remain distinct; control/v6 comparisons still require the unchanged `2e-13`; native/profile exactness and the conditional `PASS_WITH_ANALYTIC_REFERENCE_CONTROL` contract remain intact.

Parsed v7/v8 plan equality check passed after only package/import paths, literal profile-05→profile-06 generated-input paths and the added explanatory diagnostic note. All 26 jobs and their order, binaries, numerical arguments, resource limits and environment are unchanged. Assessment still reads completed profile-04 first; its failure stops the supervisor before later jobs. Subsequent kernels remain 64-angle L128/384/1536 diagnostics, with separate YSpec waveform jobs. No complete high-degree global-kernel or scientific-accuracy claim follows from this source correction.

Parent may retry under a fresh **profile-06** task using the reviewed I6 supervisor after checking exact deployed I8 script/plan hashes and existing binary/evidence receipts. Preserve failed profile-05 and every earlier task; use exclusive detached logging. The maximum 50-thread userspace watchdog, serial numerical jobs, one numerical thread and CPU-0 affinity remain. Actual assessment success is still required before high-degree dispatch, and userspace polling remains distinct from kernel aggregate resource enforcement.

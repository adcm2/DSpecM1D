# DSpec v6 receiver constructor stage map

The v6 source retains v5's three nested receiver-constructor timers. Each scope runs once per receiver and modal block while the frequency workspace is constructed:

| Profile stage | Work measured |
|---|---|
| `receiver_at_pole_closed_form` | Build the closed-form pole receiver basis. |
| `receiver_local_weight_materialization` | Form component and order receiver weights. |
| `receiver_basis_mix_and_unit_adjoint_projection` | Mix weights through the adjoint basis transform and project onto the unit-adjoint basis. |

These names and mappings are also present in `source-v6-frozen/stage-to-math-map.json`. The generic Wigner estimator follow-up does not change the receiver formulas or these timings.

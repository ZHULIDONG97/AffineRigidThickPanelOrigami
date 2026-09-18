# AffineRigidThickPanelOrigami

MATLAB examples of thick-panel folding with affine-metric rigidity constraints,
primarily for the **SAMIO2026 conference**.

## Entry points

| Script | Purpose |
| --- | --- |
| `demo_miura_drag.m` | Four-panel node dragging with LM correction and a hinge-angle limit |
| `FourFaceChirality.m` | Static four-panel folding with chirality constraints |
| `twofacechirality.m` | Static two-panel folding with chirality constraints |

The static examples require Optimization Toolbox (`fmincon`). The drag demo
requires only MATLAB.

## Run

Open this repository in MATLAB and run:

```matlab
demo_miura_drag
```

Adjust the drag node, hinge axis, angular increment, step count, regularization,
and angle limit in the settings section of `demo_miura_drag.m`.

LM uses the same stopping criteria as RigidOrigamiSimulator: scaled RMS
residual `1e-8`, scaled gradient `1e-6`, or relative step `1e-6`.
Each step generates node 18's target coordinates by rotating its **current
corrected position by 1 degree** around hinge 17–24. The difference between the
target and current coordinates is the drag perturbation. The corrected motion
need not be exactly 1 degree.

An angle-limited search bisects this Cartesian perturbation if needed, then
calls LM once. No LM call is made if no admissible positive perturbation is found.
Positive LM exits are accepted; the maximum residual is reported without an additional
`1e-8` acceptance requirement. An LM failure or a corrected trial crossing the
hinge-angle limit stops at the previous accepted state. The final angle may
remain below the limit because corrected states are not bisected or interpolated.

## Reference data

`four_panel_reference.mat` is required by the drag demo. It stores the initial
node coordinates, panel connectivity, master and fixed nodes, and affine
displacement mapping. Keep it beside the script. The two static examples
construct their own models.

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

Adjust the drag node, step length, step count, regularization, and angle limit
in the settings section of `demo_miura_drag.m`.

## Reference data

`four_panel_reference.mat` is required by the drag demo. It stores the initial
node coordinates, panel connectivity, master and fixed nodes, and affine
displacement mapping. Keep it beside the script. The two static examples
construct their own models.

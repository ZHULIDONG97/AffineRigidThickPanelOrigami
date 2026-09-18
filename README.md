# Thick-panel affine-metric examples for SAMIO2026

This folder is primarily intended for the **SAMIO2026 conference**. It contains
MATLAB examples of thick-panel folding with affine-metric rigidity constraints,
including a four-panel Miura node-drag demonstration.

## Entry points

| Script | Purpose | Requires `four_panel_reference.mat` | Requirements |
| --- | --- | --- | --- |
| `demo_miura_drag.m` | Drag a master node along the current panel normal, with LM correction and a hinge-angle limit | Yes | MATLAB |
| `FourFaceChirality.m` | Solve a static four-panel configuration with chirality and fold-branch inequalities | No; constructs its model in the script | MATLAB + Optimization Toolbox (`fmincon`) |
| `twofacechirality.m` | Solve a static two-panel configuration with chirality and fold-branch inequalities | No; constructs its model in the script | MATLAB + Optimization Toolbox (`fmincon`) |

Set the MATLAB current folder to the directory containing these scripts and run:

```matlab
demo_miura_drag
```

The drag demo uses the following settings near the top of the script:

```matlab
dragNode = 18;
dragStepLength = 0.05;
nDragSteps = 130;
regularization = 1e-7;
maximumRotationAngle = deg2rad(180);
```

## Required reference data

**Keep `four_panel_reference.mat` beside `demo_miura_drag.m`.** It is an active
input dependency, not a saved simulation result. The drag script resolves its
path relative to its own location and explicitly loads all eight fields:

| Field | Size | Use |
| --- | --- | --- |
| `x0`, `y0`, `z0` | Each `24 x 1` | Initial coordinates of the 24 nodes |
| `B` | `24 x 8` | Affine displacement basis in the reordered node coordinates |
| `T` | `24 x 24` | Node permutation matrix; `T\B` maps master displacements to the original node order |
| `mainPointIndices` | `1 x 8` | Master nodes `[5 16 17 24 6 18 1 19]` |
| `fix` | `1 x 3` | Fixed nodes `[16 17 24]` |
| `panels` | `4 x 1` cell | Face connectivity of the four solid panels |

In the full `MultiAffineBodyDynamicsCodes` source tree, four tests in
`../tests/test_affine_metric_constraints.m` also load this file
to verify constraint selection, the analytic Jacobian, drag correction, and
the constraint-evaluation pipeline. Removing it breaks both the drag demo and
those tests. The current code does not regenerate or overwrite this file;
changing geometry in `FourFaceChirality.m` does not update the saved reference.

## Drag algorithm and limits

The script starts from `x0,y0,z0`, presets 24 solid-panel metric components,
retains 20 independent `A'*A=I` constraints, recalculates the normal of face
`[17 18 24]` after every accepted correction, and applies regularized LM
without explicit chirality or penetration constraints. If a requested step
is too large, it is automatically halved until LM correction succeeds.
The four physical hinges are checked before every trial state is accepted.
Their signed angles are continuously unwrapped, so limits at or above 180
degrees remain detectable. If a trial crosses `maximumRotationAngle`, the
traction substep is bisected and LM is solved again until the last admissible
angle lies within 0.1 degree below the limit. That partial step is accepted and
the drag simulation then stops; nodal coordinates are never interpolated.
The red arrow marks the traction direction. Coordinate axes and x/y/z history
plots are omitted.

`dragStepLength` is the target perturbation supplied to LM at each displayed
step. Regularization means the corrected node displacement need not equal this
target exactly.

## Verification

The shared six-test suite belongs to the full `MultiAffineBodyDynamicsCodes`
source tree and is not included in this standalone thick-panel repository.
When that source tree is available, run it from the `thick` directory:

```matlab
results = runtests(fullfile('..','tests','test_affine_metric_constraints.m'));
assertSuccess(results);
```

All six tests passed in MATLAB R2025b on 2026-09-18, including the four tests
that load the reference data.

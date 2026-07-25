# A Bistable Model for Durable Remission in Primary Aldosteronism

MATLAB code and generated outputs accompanying:

> L. H. Carter, "A Bistable Model for Durable Remission in Primary Aldosteronism" (preprint).

The model is a reduced-order dynamical system for the renin-aldosterone axis,
augmented by a slow, salt-supported autonomous secretory-capacity state. It
interprets durable off-drug remission in primary aldosteronism as a
basin-to-basin transition in a bistable system. The remitter/non-remitter
boundary is organized by a saddle-node bifurcation in the
suppression-resistant capacity parameter, `z_fixed`.

This repository contains the MATLAB source, numerical result files, and
publication-ready figures needed to reproduce the figures and numerical tables
in manuscript version 6.

---

## Requirements

- MATLAB R2024b or a compatible recent release
- Optimization Toolbox for `fsolve`, used by the equilibrium and continuation
  analyses
- No Simulink dependency

Time integration uses `ode15s`, which is included with MATLAB.

The global robustness analysis evaluates 4,096 parameter sets and can take
substantially longer than the individual trajectory and comparison figures.
The two-parameter fold continuation is also more computationally intensive
than the other figure scripts.

---

## Core model files

These files provide shared model equations and parameter sets. They do not
generate figures on their own.

| File | Role |
| --- | --- |
| `pa_model.m` | Right-hand side of the full four-state model `[r; a; z_p; s]`. Total autonomous capacity is `z = z_fixed + z_p`. |
| `pa_params.m` | Nominal parameter set reported in Table 4 of the manuscript. |
| `pa_params_legacy.m` | Legacy parameter configuration used by the four-scenario mechanism experiment and functional-form analysis. |
| `pa_model_mono.m` | Monostable no-feedback comparator, with reversible and irreversible variants. |
| `pa_params_mono.m` | Parameter configuration for the monostable comparators. |
| `pa_model_turnover.m` | Two-compartment direct-turnover comparator with resistant and drug-sensitive capacity. |

---

## Manuscript figures

The figure numbers below follow manuscript version 6.

| Figure | Script | Generated file(s) | Description |
| --- | --- | --- | --- |
| Figure 1 | `pa_conceptual_summary.m` | `pa_conceptual_summary.png`, `pa_conceptual_summary.pdf` | Biology-first summary of the proposed disease loop, its interruption during baxdrostat treatment, and the possible outcomes after withdrawal. |
| Figure 2 | `pa_experiment.m` | `pa_experiment.png` | Finite drug-pulse experiment under four feedback combinations, showing that durable bistability tracks the slow salt bridge. |
| Figure 3 | `pa_form_sensitivity.m` | `pa_form_sensitivity.png` | Tests six increasing, saturating bridge functions to show that bistability is not specific to a single Hill-function choice. |
| Figure 4 | `pa_global_robustness.m` | `pa_global_robustness.png`, `pa_global_robustness.pdf` | Simultaneous ±20% variation of 16 equilibrium parameters in a reproducible, fixed-seed, 4,096-point Latin-hypercube design. |
| Figure 5 | `pa_bifurcation_structure.m` | `pa_bifurcation_structure.png` | One-parameter equilibrium branches and two-parameter fold structure in `(z_fixed, sigma)` and `(z_fixed, g_z)`. |
| Figure 6 | `pa_phase_map.m` | `pa_phase_map.png` | Remission map over treatment duration and suppression-resistant capacity, including sodium and drug-potency comparisons. |
| Figure 7 | `pa_mono_compare.m` | `pa_mono_compare.png` | Bistable model compared with two monostable no-feedback alternatives. |
| Figure 8 | `pa_turnover_compare.m` | `pa_turnover_compare.png` | Bistable salt-bridge model compared with a direct-turnover competitor, highlighting sharp versus smooth remission boundaries. |

Pre-generated versions of all eight manuscript figures are included.

### Figure 5 prerequisite

`pa_bifurcation_structure.m` reads `pa_two_parameter_fold_data.mat`. Generate
that file first by running:

```matlab
pa_two_parameter_folds
```

This continuation analysis also writes:

- `pa_two_parameter_folds.png`
- `pa_two_parameter_folds.pdf`
- `pa_two_parameter_validation.csv`
- `pa_two_parameter_fold_data.mat`

The first two are diagnostic continuation figures rather than manuscript
figures. The validation CSV records representative equilibrium counts,
stability classifications, and residual checks. The MAT file supplies the
continuation and classification data used to construct Figure 5.

---

## Manuscript tables and numerical results

Tables 1, 3, and 4 are explanatory or parameter tables written directly in the
manuscript. Their model parameter values are defined in `pa_params.m`.

| Manuscript output | Source | Description |
| --- | --- | --- |
| Table 2 | `pa_rigor_4state.m` | Prints the equilibria, full 4-by-4 Jacobian eigenvalues, stability classifications, and characteristic timescales used in the table. It also writes the diagnostic figure `pa_rigor_4state.png`. |
| Table 5 | `pa_global_robustness.m` | Computes the simultaneous-parameter robustness estimates, confidence limits, fold distributions, and sodium comparison reported in the table. |

The global robustness analysis writes the following reproducibility files:

| File | Contents |
| --- | --- |
| `pa_global_robustness_summary.csv` | Summary estimates and uncertainty intervals reported in Table 5. |
| `pa_global_robustness_associations.csv` | Standardized main-effect associations with the fold threshold. |
| `pa_global_robustness_samples.csv` | Parameter multipliers and sample-level results for all 4,096 runs. |
| `pa_global_robustness_data.mat` | Complete MATLAB result structure used by the analysis. |

`pa_bifurcation_sensitivity.m` is retained as a supplementary, local
one-parameter-at-a-time sensitivity check. It is not the source of the v6
simultaneous-parameter robustness figure or Table 5.

---

## Reproducing all manuscript outputs

Start MATLAB in the repository directory, or add the directory to the MATLAB
path. Run the following commands in order:

```matlab
pa_conceptual_summary       % Figure 1
pa_experiment               % Figure 2
pa_form_sensitivity         % Figure 3
pa_global_robustness        % Figure 4 and Table 5 data

pa_two_parameter_folds      % Generate continuation data for Figure 5
pa_bifurcation_structure    % Figure 5

pa_phase_map                % Figure 6
pa_mono_compare             % Figure 7
pa_turnover_compare         % Figure 8

pa_rigor_4state             % Table 2 eigenvalues and diagnostic figure
```

The scripts save their output beside the source files. Existing generated files
with the same names are overwritten.

For the exact fixed-seed robustness results reported in the manuscript, the
repository also includes the generated CSV and MAT files listed above.

---

## Model summary

The full model has four dimensionless states and uses days as its time unit:

- `r`: lumped renin/angiotensin-II drive, a fast variable
- `a`: circulating aldosterone activity, a fast variable and the principal
  clinical observable
- `z_p`: plastic autonomous secretory capacity, a slow variable
- `s`: effective sodium/volume/mineralocorticoid-receptor state, an
  intermediate-timescale variable

Total autonomous capacity is

```text
z = z_fixed + z_p
```

Here, `z_fixed` represents suppression-resistant capacity that is not
appreciably removed by output suppression on the modeled timescale, whereas
`z_p` is plastic capacity sustained by the proposed salt bridge.

The parameter values are nondimensional and are not fitted to clinical data.
They place the system in a regime where bistability is possible; the analyses
then examine the resulting dynamics, bifurcation structure, predictions, and
robustness.

---

## License

Released under the MIT License; see `LICENSE`.

## Citation

If you use this code, please cite the associated manuscript and archived code
release:

L. H. Carter, "A Bistable Model for Durable Remission in Primary
Aldosteronism," Zenodo preprint, 2026.

Manuscript: https://doi.org/10.5281/zenodo.21170037

Code archive: https://doi.org/10.5281/zenodo.2117

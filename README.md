# Durable Remission in Primary Aldosteronism as a Saddle-Node Bifurcation

MATLAB code accompanying the paper:

> L. H. Carter, \\\*A Bistable Model of Primary Aldosteronism: Durable Remission
> as a Saddle-Node Bifurcation\\\* (preprint).

A reduced-order dynamical model of the renin–aldosterone system augmented by a
slow autonomous secretory-capacity state (the "salt bridge"). The model recasts
durable off-drug remission in primary aldosteronism as a basin transition in a
bistable system, with the remitter/non-remitter boundary a saddle-node
bifurcation in the per-patient irreversible-capacity parameter `z\\\_fixed`.

\---

## Requirements

* MATLAB (developed on R2024b; no Simulink required)
* **Optimization Toolbox** — required for `fsolve`, used by the equilibrium /
continuation routine `pa\\\_rigor\\\_4state.m`
* Stiff integration uses `ode15s`, which ships with base MATLAB
* No other toolboxes required

The scripts with no `fsolve` dependency run on base MATLAB alone.

\---

## Core model files

These are shared infrastructure, called by the scripts below; they do not
produce output on their own.

|File|Role|
|-|-|
|`pa\\\_model.m`|Right-hand side of the full four-state model `\\\[r; a; z\\\_p; s]`, Eqs. (1)–(4). Total capacity `z = z\\\_fixed + z\\\_p`. Flags select the dynamic vs. algebraic salt signal, the fast output loop, and the salt bridge.|
|`pa\\\_params.m`|Nominal parameter set (Table 4). Returns the default parameter struct used throughout.|
|`pa\\\_params\\\_legacy.m`|Parameters with post-scaffold features off (single-`z`, algebraic salt signal, `z\\\_fixed = 0`): used by the four-scenario experiment.|
|`pa\\\_model\\\_mono.m`|Monostable no-feedback comparator RHS (Eq. 5), variants M1 (reversible) and M2 (irreversible ratchet).|
|`pa\\\_params\\\_mono.m`|Parameters for the monostable comparator; `mode = 'reversible'` or `'irreversible'`.|
|`pa\\\_model\\\_turnover.m`|Two-compartment direct-turnover comparator RHS (Eq. 6): resistant `z\\\_R` + drug-ablated sensitive `z\\\_S`.|

\---

## Scripts that generate paper figures

Each script writes a `.png` of the same base name and prints its numerical
summary to the console.

|Script|Paper output|Description|
|-|-|-|
|`pa\\\_experiment.m`|**Figure 1** (`pa\\\_experiment.png`)|Finite baxdrostat pulse under the four feedback combinations (both / fast-loop only / slow-bridge only / neither). Confirms a disease state and durable remission appear only with the slow salt bridge.|
|`pa\\\_form\\\_sensitivity.m`|**Figure 2** (`pa\\\_form\\\_sensitivity.png`)|Robustness of bistability to the bridge functional form `H(s)`: six increasing/saturating shapes over a range of bridge strength `g\\\_z ∈ \\\[0.2, 6]`, plus the sustain-vs-decay crossing geometry.|
|`pa\\\_rigor\\\_4state.m`|**Figure 3** (`pa\\\_rigor\\\_4state.png`) **and Table 3**|Full four-state rigor pass. Part 1 prints the 4×4 Jacobian eigenvalues at the three equilibria (**Table 3** data). Part 2 runs predictor–corrector continuation in `z\\\_fixed` and locates the saddle-node at `z\\\_fixed\\\* = 0.367`. Part 3 sweeps sodium `σ` (right panel).|
|`pa\\\_phase\\\_map.m`|**Figure 4** (`pa\\\_phase\\\_map.png`)|Remission map over treatment duration × `z\\\_fixed`, with the sodium and drug-potency comparison panels.|
|`pa\\\_mono\\\_compare.m`|**Figure 5** (`pa\\\_mono\\\_compare.png`)|Bistable model vs. no-feedback comparators M1 and M2, scored on capacity `z`. M1 never durably remits; M2 remits for all `z\\\_fixed`; only the bistable model gives durability + a rarity threshold.|
|`pa\\\_turnover\\\_compare.m`|**Figure 6** (`pa\\\_turnover\\\_compare.png`)|Bistable model vs. the fair two-compartment direct-turnover competitor. Sharp (discontinuous) vs. smooth boundary; boundary-sharpness differs by \~2 orders of magnitude.|

\---

## Script that generates paper table data

|Script|Paper output|Description|
|-|-|-|
|`pa\\\_bifurcation\\\_sensitivity.m`|**Table 5**|Local one-at-a-time sensitivity of the saddle-node threshold `z\\\_fixed\\\*` to ±10% / ±20% perturbations of each parameter. Prints the table; confirms the bifurcation persists in every case and that `τ\\\_s`, `ε` (rate controls) leave the threshold essentially unchanged.|

(The Table 3 eigenvalue data is produced by `pa\\\_rigor\\\_4state.m`, Part 1 — see above.)

\---

## Reproducing the paper figures

From the repository root, in MATLAB:

```matlab
pa\\\_experiment            % Figure 1
pa\\\_form\\\_sensitivity      % Figure 2
pa\\\_rigor\\\_4state          % Figure 3  + Table 3 (console)
pa\\\_phase\\\_map             % Figure 4
pa\\\_mono\\\_compare          % Figure 5
pa\\\_turnover\\\_compare      % Figure 6
pa\\\_bifurcation\\\_sensitivity   % Table 5 (console)
```

Each figure script saves a `.png` beside itself and prints its numerical
summary to the console. Pre-generated `.png` copies of the six figures are
included in the repository.

\---

## Model summary

Four dimensionless states, integrated in days with `ode15s` (stiff, because of
the three-tier timescale separation):

* `r`   — lumped renin / angiotensin II drive (fast, \~hours)
* `a`   — circulating aldosterone activity (fast, \~hours); the clinical observable
* `z\\\_p` — plastic autonomous secretory capacity (slow, \~weeks–months)
* `s`   — effective sodium–volume / mineralocorticoid-receptor state (intermediate, \~weeks)

Total autonomous capacity `z = z\\\_fixed + z\\\_p`, where `z\\\_fixed` is the
irreversible, mutation-locked per-patient parameter that the drug and salt
bridge cannot remove.

Parameters are nondimensional and **not fitted to data**; they place the system
in a regime where bistability is possible, and the analysis then varies them to
determine the dynamical consequences. Nominal values are in Table 4 of the
paper.

\---

## License

Released under the MIT License; see `LICENSE`.

## Citation

If you use this code, please cite the associated manuscript and archived code release:

L. H. Carter, "A Bistable Model of Primary Aldosteronism: Durable Remission
as a Saddle-Node Bifurcation," Zenodo preprint, 2026.
https://doi.org/10.5281/zenodo.21170038

Code archive:
https://doi.org/10.5281/zenodo.21149177

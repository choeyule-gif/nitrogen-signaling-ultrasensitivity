# Quantifying the sources of ultrasensitivity in bacterial nitrogen signaling

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22731768.svg)](https://doi.org/10.5281/zenodo.22731768)

Code and data for

> Choi Y. *Quantifying the sources of ultrasensitivity in bacterial nitrogen signaling.*
> Manuscript prepared for PLOS Computational Biology (2026).

The manuscript analyses a distributive modification cycle in which an effector partitions a single
bifunctional enzyme pool between opposing activities. The shared enzyme cancels from the ratios of
adjacent free-target states, so the local response sensitivity factorises into a ladder factor and
an effector-switch factor; the framework is then applied to published measurements of the
*Escherichia coli* nitrogen-signaling cascade (GlnD–PII and ATase–GS).

Everything numerical in the paper is produced here: the verifications of the derived statements,
the reanalysis of the published titration, the cascade accounting, the values behind every figure
panel, and the Supporting Information package. `tests/test_reproduce.py` checks the output against
the values quoted in the text, and `REPRODUCIBILITY.md` records the comparison statement by
statement.

**Repository:** <https://github.com/choeyule-gif/nitrogen-signaling-ultrasensitivity>

**Archived release:** [v1.0.0 (DOI: 10.5281/zenodo.22731769)](https://doi.org/10.5281/zenodo.22731769). The [concept DOI 10.5281/zenodo.22731768](https://doi.org/10.5281/zenodo.22731768) always resolves to the latest archived version.

## Layout

| path | what it holds |
|---|---|
| `esbm/` | the model as a small Python package: the compartmental matrix and its kernel steady state, the closed forms proved in the paper, and the loaders for `data/` |
| `data/` | inputs only: the twelve digitized points of the published titration and the published constants (`data/README.md` gives the sources and the uncertainty model) |
| `scripts/` | one script per result; each prints what it reproduces and writes `results/<name>.json` |
| `results/` | the generated numbers (`*.json`), the run logs, and `results/expected/`, the same files as produced on the machine the manuscript was written from |
| `figdata/` | the values behind every quantitative panel of main Figs 2–6 and S1–S5 Figs, one CSV per panel, written by `scripts/export_figure_data.py` |
| `matlab/` | the MATLAB scripts that draw the figures from `figdata/` (`matlab/README.md`) |
| `tests/` | comparison of `results/*.json` with the manuscript |
| `extras/` | the Re *F*(iω) figure of the companion stability manuscript, which the Discussion cites |

Generated and not under version control: `figures/` (every rendered figure), `S1_Data/` (the
Supporting Information package, assembled by `scripts/make_s1_data.py`), `results/panels/` and the
run logs.

## Requirements

Python ≥ 3.10 with `numpy`, `scipy`, `sympy`, `mpmath` and `matplotlib`, and MATLAB ≥ R2021b for
the figures. `requirements.txt` pins the versions the manuscript was produced with (CPython
3.14.6, numpy 2.5.2, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0, matplotlib 3.11.1; macOS 26 on
Apple silicon). No compiled code, no GPU, no network access.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Reproducing the paper

```bash
python run_all.py          # every number and every figure data file (~45 min)
python run_all.py --quick  # skips the four long searches (~10 min)
python tests/test_reproduce.py      # or: python -m pytest -q tests
```

Then, in MATLAB:

```matlab
cd matlab
run_all_figures            % one vector PDF and one PNG per panel, into ../figures
```

`make all`, `make quick`, `make test`, `make figdata` and `make s1data` wrap the same commands.
Every script can also be run on its own from any directory (`python scripts/fig2_sites.py`); the
few that depend on another's output say so in their docstring (`fig8_budget.py` and
`sensitivity_major1.py` need `cascade_budget.py`; `fig5_integer.py` needs `fit_titration.py` and
`mcmc_h.py`; `fig6_uniqueness.py` needs `counterexamples.py`; both export scripts need the
`results/` of a full run). All random draws are seeded; the seeds are in the scripts.

## What each script reproduces

Equation numbers refer to the manuscript. Some script and CSV names preserve the numbering of an
earlier draft; the mapping to the final six-figure manuscript is shown below. The right-hand column is the statement in the
text that the script's printed output and `results/<name>.json` are checked against.

| final manuscript figure | analysis/plot identifier |
|---|---|
| Fig 1 | explanatory pathway schematic; no numerical analysis |
| Fig 2 | `Fig1.m`, `Fig1B_theta.csv`, `Fig1C_local_hill.csv`, plus the concentration-dependence panels |
| Fig 3 | `Fig2.m` and `Fig2*.csv` |
| Fig 4 | `Fig3.m` and `Fig3*.csv` |
| Fig 5 | `Fig4.m`, `Fig4*.csv`, plus the direct-glutamine-route panel |
| Fig 6 | `Fig5.m` and `Fig5*.csv` |

### Verifications of the derived statements

| script | manuscript statement |
|---|---|
| `verify_symbolic.py` | computer-algebra checks with `n` symbolic: the chain survives three simultaneous dead ends; the Möbius relation and the quartic reduction with its extreme coefficients; the floored exponent-2 Hill identity of Eq. (8); the two-form chain with variable *l*₁/*l*₂ and the linear/quadratic form of *l*₂ under condition (i); Eq. (S19) on 4,000 random ladders; the sign counts of det *J*; the active bridging complex at *n* = 1 |
| `verify_adair.py` | Eq. (2) against the kernel steady state over 1,500 parameter sets, in double and in 40-digit arithmetic; the two-form network over a thousandfold range of enzyme; *n*_H = *h* exactly for a cooperatively bound effector; ν → 1.999 at *n* = 2 with disparate site constants |
| `verify_nu_bound.py` | 0 < ν ≤ *n* in 50–80-digit arithmetic: no exceedance over 28 decades of the ladder weights, the bound approached to 10⁻²⁰, and ν ≡ 1 for a binomial ladder to 6 × 10⁻⁴³ |
| `verify_deadend_60digit.py` | the invariant in 60-digit arithmetic with a dead-end complex swept over eleven decades, and the drift in exact proportion to *T_E*/*T_L* (Eq. 5) |
| `sequestration_sweeps.py` | the sequestration correction: 0.988 → 0.976 on doubling the enzyme, the affinity sweep, and the target sweep from *K_m*/60 to 100 *K_m* |
| `deficiency_acr.py` | the deficiency and the Shinar–Feinberg structural test at *n* = 1 and in general; absolute concentration robustness in the complex-mediated design; the buffering example |
| `counterexamples.py` | S5 Fig: three positive steady states with condition (i) broken (spread 2.89 in the catalytic ratios, *n* = 4) and with condition (ii) broken (*n* = 4, 5, 6), every root refined in 60-digit arithmetic |
| `mono_adversarial.py` | the adversarial search for a non-monotone total effector: 21 pairs (*n*, *h*), every optimum positive in 50-digit arithmetic |

### Reanalysis of the published data

| script | manuscript statement |
|---|---|
| `fit_titration.py` | the bounded four-parameter fit, the fixed-exponent comparison, Shapiro–Wilk, the profile interval in *h*, the standardised residuals, the floorless and refractory-fraction alternatives, the end-point drops, and the plateaus with the local logit slope |
| `mcmc_h.py` | the posterior of the exponent: acceptance, Gelman–Rubin, effective sample size, median and interval |
| `cascade_budget.py` | the readout factors, composites and residuals at the three PII concentrations, both readings of the reverse input, the supremum over the downstream half-point, and the bootstrap |
| `sensitivity_major1.py` | how much error on the three reported cascade coefficients the intermediate-PII maximum survives |
| `nu2_major2.py` | the cascade composite when the downstream ladder is non-binomial (S4 Fig) |
| `mixture_bias.py` | the composition of the heterotrimer preparation and the bias it can put on the comparison |
| `conversion_factor.py` | the conversion between a reported coefficient and the local ladder factor, on three ladder families and on the full two-coordinate locus |
| `swing_profile.py` | the plateaus, the regulatory swing and its profile lower endpoint, and the implied inhibition of the transferase |
| `direct_input.py` | the additivity of the two log-slopes and the gain a direct glutamine input on the adenylyltransferase can supply (a re-implementation; see the note below) |
| `ratios_table.py` | the enzyme-to-effector ratios of Eq. (5) for each experiment |

### Panel values and figures

`fig1_invariance.py`, `fig2_sites.py`, `fig3_plateau.py`, `fig4_plane.py`, `fig5_integer.py`,
`fig6_uniqueness.py`, `fig7_tests.py` and `fig8_budget.py` compute the underlying quantities,
write them to `results/panels/*.csv`, and render reference versions of the panels in matplotlib.
`scripts/export_figure_data.py` and `scripts/export_extra_figure_data.py` then assemble
`figdata/*.csv`, one file per quantitative panel of the manuscript's main Figs 2–6 and S1–S5 Figs, and
`results/revision_numbers.json`, which holds every quantity the text quotes. The figures in the
paper are drawn from those CSV files in MATLAB; nothing is recomputed at plot time.

## Notes on reproducibility

* Every number in the Results, Methods and figure captions is produced by a script here and
  compared in `tests/test_reproduce.py` at the rounding used in the text. The long searches
  (`mono_adversarial.py`, `counterexamples.py`, `mcmc_h.py`, the profile in `fit_titration.py`)
  are seeded and deterministic on one machine; across numpy/scipy versions their qualitative
  conclusions are what the tests check.
* Steady states are always the kernel of the compartmental matrix (an exact linear solve), never
  a multistart Newton iteration; in this network Newton iteration produces spurious roots from
  cancellations of order *u^n*. In double precision the kernel solve is limited by the condition
  number of the matrix, so `verify_adair.py` reports the wide random family in double precision,
  a well-conditioned family in double precision, and the wide family with the kernel solved in
  40-digit arithmetic, which is the figure that tests the algebra rather than the arithmetic.
* `conversion_factor.py` evaluates the conversion on the three ladder families named in the
  Methods and, separately, along the full two-coordinate locus of three-site ladders that
  reproduce the observed ratio; the locus range is wider than the three-family envelope and is
  reported as such in `results/conversion_factor.json`.
* `direct_input.py` is a re-implementation: the original interactive calculation behind three
  additivity numbers was not preserved. The script states its construction and reproduces the
  exactness of the additivity and the gain at the measured constants.
* The twelve titration points were read from the printed figure; the uncertainty model
  (±0.06 in ordinate, ±6 % in abscissa) is stated in `data/README.md` and used by every
  perturbation analysis.

## License and citation

Code: MIT (`LICENSE`). Please cite the manuscript and the Zenodo record (`CITATION.cff`).

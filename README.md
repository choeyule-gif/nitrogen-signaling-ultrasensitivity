# Discriminating mechanisms of ultrasensitivity in bacterial nitrogen signaling

Code, numerical outputs and manuscript sources for Yule Choi's nitrogen-signaling study. This is a conditional analysis of mechanistic equivalence and measurement design; the calculations do not establish a unique native GlnD–PII–GlnE–GS mechanism.

**Archived release:** [v1.1.0](https://doi.org/10.5281/zenodo.22735182). **Working revision:** adds independent paired-paralog data, assay interval calculations and integrated validation. These additions are not contained in the v1.1.0 archive. The manuscript and submission package are being revised; the older release ZIP is not the current working manuscript.

## Reproduce

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
python tests/run_research_suite.py --regenerate
```

The command regenerates the baseline, matched-model, finite-mechanism, regulatory-state, robust-design and closed-cascade numerical layers and runs their checks, followed by the paired-reporter checks. Run without `--regenerate` to check saved baseline results while freshly solving extension identities. See `REPRODUCIBILITY.md` for scope and independent-data extraction. A passing calculation verifies its stated model, not physiological truth or global parameter optimality.

## Research layout

| Location | Content |
|---|---|
| `esbm/`, `scripts/`, `data/`, `results/`, `figdata/` | Baseline theory, digitized titration, calculations and plot inputs |
| `extras/matched_models/` | Matched titrations, different ladders and downstream calibration |
| `extras/regulatory_states/` | Explicit GlnD/GlnE regulatory-state equivalences |
| `extras/finite_mechanism/` | Finite reaction kinetics and assumption boundaries |
| `extras/closed_cascade/` | Protein-conserving cascades, topology and free/total input bounds |
| `extras/robust_design/` | Nuisance-aware occupancy and binding designs |
| `extras/paired_paralogs/` | Independent published GlnB/GlnK data and calibrated reporter design |
| `extras/revision_audited/`, `matlab/` | Figure calculations and rendering inputs; some filenames retain older numbering |
| `manuscript/` | Working manuscript, appendix and supporting figures |
| `validation/`, `tests/` | Clean-environment evidence and executable checks |

## Current figure mapping

The working manuscript still contains 11 main figures pending compression. Figure 1 is the author's pathway diagram. Figure 2 covers architecture and ladder constraints; Figure 3 sensitivity; Figure 4 titration; Figure 5 matched models; Figure 6 cascade comparison; Figure 7 GlnD states; Figure 8 GlnE states; Figure 9 finite reaction boundaries; Figure 10 closed cascades; Figure 11 measurement design. S1–S10 figures remain separate supporting items. Earlier script names are not submission figure numbers. The paired-paralog analysis currently appears as text and equations, with diagnostic plots retained in its research folder.

Code is MIT licensed. The Gosztolai et al. source workbook and model retain their CC BY 4.0 license and attribution, documented in `extras/paired_paralogs/README.md`. Other published measurements retain their source provenance in `data/README.md`.

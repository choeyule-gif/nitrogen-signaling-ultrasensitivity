# data/

Inputs only. Nothing in this directory is produced by the code.

| file | content | source |
|---|---|---|
| `jiang1998_fig2B_digitised.csv` | twelve points (glutamine, mM; mean UMP per PII trimer, 0–3) read from Fig. 2B of Jiang, Peliska & Ninfa, *Biochemistry* 37:12782 (1998), at PII 0.5 µM and UTase/UR 0.05 µM | digitised by eye from the enlarged printed figure |
| `published_constants.json` | Hill coefficients, half-points, concentrations and affinities quoted from the cited papers | Ventura et al. 2010; Jiang & Ninfa 2011; Jiang et al. 1998a,b |

**Digitisation uncertainty model** (used by every bootstrap in `scripts/`): each ordinate is
perturbed by ±0.06 (12 % of the 0.5 tick spacing) and each abscissa by ±6 %, independently,
uniform in the bootstrap of `cascade_budget.py` and Gaussian in the bootstrap of
`fit_titration.py`; both are stated in the manuscript's Methods.

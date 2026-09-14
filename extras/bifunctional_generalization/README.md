# General graph enclosure and independent AceK examples

This extension separates three model-based analyses from independent-activity checks.

1. `GENERAL_METHOD.md` derives a graph-based cancellation criterion and a
   positive-polynomial enclosure, followed by conservation-aware observation
   linear programs. The matrix-tree theorem and polynomial basis identities
   are existing mathematical foundations. The proposed contribution is the
   executable composition for a declared biochemical class.
2. `acek_kinetic_ambiguity.py` constructs ordered and random initial-rate
   mechanisms using six published AceK kinetic parameter summaries: WT and two
   mutants, kinase and phosphatase. These are source-summary calibrations,
   not six new raw-data fits or independent biological validation datasets.
3. `acek_certificate.py` retains the singly and doubly enzyme-bound target
   complexes of the published nine-target-state AceK–IDH model. It tests
   coefficient cancellation, finite conserved kinetics, and the general
   enclosure independently of the nitrogen cascade implementation.

Numerical work was run on a dedicated Linux computation server. The source data provenance and cross-assay
limitations are in `SOURCE_AUDIT.md`.

## Current verified results

The ordered substrate-binding mechanism realizes the published sequential
initial-rate law exactly. A finite random-binding network approaches that law
as binding becomes fast relative to catalysis. At ratio 1,000, its largest
relative deviation over 625 ATP × protein input pairs per kinetic summary is
0.000979–0.001143, while a low-ATP binding measurement predicts an occupancy
contrast of 0.495–0.497 between mechanisms. This is a constructive interpretation
and prospective design example, not proof that either network is the native
AceK mechanism. The approximation dependence is explicitly calculated at
ratios 1, 10, 100 and 1,000. Twelve independent matrix-exponential solutions
of the finite reaction networks match their stationary distributions.

The full dimeric target network has coefficient rank three, so its free target
distribution does not cancel enzyme concentration exactly. Positive quadratic
weights reproduce the independent stationary-generator solution to 1.8e-14.
Column-scaled bound-state dwell-time families retain the free distribution at
matched free enzyme; they are constructed classes, not measured confidence
regions for the source kinetics. Nine concentration settings and 900 such
conserved realizations lie within the general enclosure. Two independent
conserved mass-action BDF integrations match all species to better than 1e-9.

A second prospective observation-conditioned decision is now available:
at 1 uM IDH dimers and 0.05 uM AceK, the general interval width for the total
unphosphorylated-dimer fraction is 0.056995. If a matched mean phosphate count
lies in [1.435,1.475], the width decreases to 0.036257 and crosses a 0.05
tolerance. Those count limits are declared prospective thresholds, not
existing measurements. Among 1,000 constructed realizations at those totals,
790 satisfy the gate and none violates the conditional interval. Two further
enzyme totals are checked with 1,000 realizations each.

The LaPorte 1985 activity data remain a separately transcribed in vivo dataset.
Its activity units cannot be silently mapped onto the micromolar simulation
grid, and its culture conditions are not Miller's in vitro assay conditions.
It has not yet provided an independent matched-assay validation of the new
concentration certificate. This remaining evidence distinction must be explicit.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extras/bifunctional_generalization/code/acek_certificate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extras/bifunctional_generalization/code/acek_kinetic_ambiguity.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python extras/bifunctional_generalization/code/plot_acek.py
```

## Independent graph and activity checks

`verify_general_graph.py` constructs exact cofactors directly from graph edges and checks four cyclic graphs, a rank-one control, an infeasible observation, and the independently derived AceK coefficients. `acek_holdout.py` adds an ATPase branch without changing the kinase fit, then compares three withheld ATPase halfpoints. Their measured/predicted ratios are 1.353, 1.328 and 0.547; missing exact errors and covariance preclude significance claims. The minimal ATPase embedding is challenged by these observations, not validated as the native mechanism.

Run both scripts after the two model analyses; the integrated suite includes all four.

## Figure caption

**A second bifunctional enzyme exposes both kinetic ambiguity and the boundary
of enzyme cancellation.** (A) WT AceK phosphatase initial-rate law reconstructed
from Miller et al. (1996), Table I, and two finite mechanistic realizations.
The reference line is a published parameter-summary reconstruction, not raw
experimental points. (B) Largest relative departure of finite random binding
from the common sequential rate law, over 625 dose pairs for each of six
published kinetic summaries. The horizontal line is a computational 1%
tolerance, not an experimental confidence limit. (C) Predicted protein-binding
contrast at protein concentration Ka Kmb/Kma, versus ATP/Kma; the random
network uses binding/catalysis ratio 1,000. Lines denote source protein
variants; solid and dashed lines denote kinase and phosphatase assays.
(D) Conservation-aware interval widths in the full dimeric target network,
which includes two-enzyme-bound complexes and fails exact enzyme cancellation.
The reporter is the total unphosphorylated-dimer fraction, not a measured
physiological activity. These widths apply to the declared dwell-time class
and source coefficients. The dashed line marks a prospective 0.05 tolerance.

## Apply the generic algorithm to another graph

Run the following from `extras/bifunctional_generalization/code`. Edge tuples
contain source state, destination state and coefficients in increasing powers
of free enzyme. The example tests an exactly cancelling two-state control;
use the full AceK script for a non-cancelling biological example.

```python
from graph_enclosure import tree_coefficients, bernstein_vertices, observation_interval

edges = [(0, 1, [0, 2]), (1, 0, [0, 3])]
coefficients, exact, generator, enzyme = tree_coefficients(2, edges)
vertices = bernstein_vertices(coefficients, 0, 4)
interval = observation_interval(vertices, epsilon=0.1, reporter=[1, 0])
print(exact.rank(), interval)
```

For an enzyme-bound network, pass the unbound target indices as `free_states`.
All bound mass must be aggregated into the same target categories as the free
vertices. Supply measurements as `(weights, lower, upper)` tuples through
`observations`. An incompatible observation returns `None`; it is not a
zero-width certificate. The caller must justify concentration units, the
free-enzyme interval, the bound-fraction allowance and observation coverage.

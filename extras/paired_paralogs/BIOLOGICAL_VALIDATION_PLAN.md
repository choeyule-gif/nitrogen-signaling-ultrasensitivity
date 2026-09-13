# Evidence boundary and a prospective validation sequence

## Existing direct binding measurements

Cong X, Liu Y, Liu W, Liang X, Laganowsky A. Allosteric modulation of protein-protein interactions by individual lipid binding events. Nature Communications 8, 2203 (2017). https://doi.org/10.1038/s41467-017-02397-0

The primary article reports AmtB–GlnK binding by native MS and SPR. At 50 micromolar ADP, native MS gives a dissociation constant of 1.12 ± 0.41 micromolar; SPR gives an average of 0.65 micromolar under the corresponding buffer/detergent conditions. ADP titration changes the binding constant substantially, and individual lipid binding events also modulate interaction. These results establish direct binding assay feasibility and show why nucleotide and membrane conditions matter. They do not provide a matched GlnB/GlnK uridylylation calibration or validate occupancy reconstruction in the Gosztolai cellular time series. A reported binding constant cannot be transplanted into those cells without concentration, input and environment matching.

## Required evidence and what the present analysis can replace

| Concern | Evidence required for validation | Extension available now | Remaining limitation |
|---|---|---|---|
| Stationarity | A time course showing stable modification and partition under fixed inputs, with a predeclared tolerance | Mixture bounds themselves apply at each time; an independently calibrated dynamic response may replace a stationary curve | Existing transient peptide data do not establish either calibration |
| Free GlnB reporter | Direct free/bound reporter measurements | A quantitative upper bound on reporter sequestration can be propagated | An unmeasured bound cannot be assumed small |
| Constant relative specificity | Matched no-binding controls across assay inputs | A calibrated response curve replaces a single fixed odds ratio | Transfer of that curve still needs testing |
| Homotrimers | Composition-resolved native measurement or defined reconstitution | Protein-subunit mixture bounds survive mixed oligomers | Calibration can still change with composition; particle counts require composition |
| Bound-GlnK selectivity | Measure modification in separated bound material | Bound-modification intervals replace the zero-modification assumption | Those intervals must be measured or justified independently |
| Actual sequestration | A held-out direct binding measurement | Compare direct measurement with the calculated interval | Numerical recovery in synthetic data is not this validation |

## A test that can reject the observation model

First determine a no-binding calibration curve over the intended input range. Use defined protein composition and measure any residual reporter partition. Then measure both total paralog modification fractions and bound-pool GlnK modification under binding-enabled conditions. Fix the uncertainty model before examining a separate direct occupancy endpoint. The sharp interval in ROBUST_CALIBRATION.md predicts an allowed bound protein fraction. Compare its interval with the independently measured fraction, including the direct assay uncertainty.

Disjoint intervals reject the jointly specified calibration/observation assumptions. Overlap establishes compatibility at those conditions, not universal mechanistic correctness. Multiple held-out input settings are needed to examine transferability. Calibration controls, selectivity measurements and validation occupancy measurements have distinct roles and should not be reused as if independent evidence.

No new wet-laboratory measurements are present in this repository. The current result is a better specified, less restrictive testable assay and a documented route to biological validation, not completed biological validation.

## Published perturbation evidence for the observation limitation

Radchenko MV, Thornton J, Merrick M (2014). Association and dissociation of the GlnK–AmtB complex in response to cellular nitrogen status can occur in the absence of GlnK post-translational modification. Frontiers in Microbiology 5:731. https://doi.org/10.3389/fmicb.2014.00731

Figure 3 reports changing membrane association of unmodifiable GlnK Y51A after ammonium shock, despite unchanged modification state. Native PAGE and membrane Western blots were repeated three times. The methods specify 5–10 micrograms of membrane protein per lane. These observations provide a biological example of modification state failing to determine partition. They do not supply calibrated absolute occupancy or a matched GlnB reporter curve, so they cannot validate the quantitative reconstruction. No band intensities have been converted into absolute sequestration fractions here.

Coutts et al. (2002), *Membrane sequestration of the signal transduction protein GlnK by the ammonium transporter AmtB*, https://pmc.ncbi.nlm.nih.gov/articles/PMC125854/, additionally reports membrane-associated modified GlnK under some preparation conditions and GlnB association when GlnK is absent. This cautions against treating either strict modification selectivity or a universally free GlnB reporter as a condition-independent fact. The bounded-departure formulation is appropriate only when those bounds have experimental justification.

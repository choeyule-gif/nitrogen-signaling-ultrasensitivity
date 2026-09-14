# Response to the supplied internal reviewer assessment

This is an author revision responding to the supplied assessment of `PLOS_CB_reviewer_revision3 2`. It is not a response to an actual journal decision. The original input folder is preserved. The manuscript remains a theoretical and methodological study without new wet-laboratory measurements.

## 1. Complete a decision that changes because of the bound

The revised main Figure 3 and Results give a prospective concentration case, rather than only a generic inequality. A concurrently justified upper bound on total PII modification constrains the common free zero-UMP fraction. That constrains the concentration coefficient in the GlnE drive, without requiring independent PII sites or a known AR gain. At the declared totals, reducing GlnE from 0.10 to 0.02 micromolar gives a count-conditioned upper bound of 0.0412 against a 0.05 tolerance. The general twelve-site bound remains 0.370. Thus the extra observation changes what can be certified.

The certificate is conditional on a concurrent count gate and a supported coefficient envelope. The historical rounded point 0.65 does not establish the chosen 0.67 gate as a confidence limit. The bound concerns pairwise model differences, not error relative to biological truth. Failed certification does not imply detectability. The example is not a measured optimization of an actual assay.

## 2. Connect biochemical constraints to candidates and contrast

S2 D4 now compares the 1998 UT/UR, 2007 GlnE and 2011 coupled-system assay conditions. The 1998 UT and UR assays use 25 and 50 mM MgCl2, while the coupled system uses 10 mM. Cross-assay transport is explicitly a stress hypothesis.

Relaxing constant catalytic capacity while retaining both transported GlnD specificity constraints reduces RMSE from the earlier 0.540 stress result to about 0.0409, but the binding scale reaches an extreme boundary. This is not presented as a plausible affinity estimate. A fixed-K1 profile in Fig4D shows that equivalence persists while its maximum occupancy contrast can shrink from 0.242 to 0.0108 as fit RMSE improves from 0.0493 to 0.0411. The original 0.244 occupancy gap is therefore not generalized to all constrained candidates.

The GlnE case uses the primary supplementary EP/EG/EGP coefficients 7.5/175/1315. The main-text K_P unit and conflicting supplementary unit are both reported and evaluated. The chosen K_G span is source motivated, not a confidence interval. Candidate and nuisance values have not become a validated native parameter region.

## 3. Test the full observation and estimation process

Fig4A–C and D5 replace the main known-variance budget with finite-depletion binding measurements, three glutamine conditions, a shared baseline, estimated enzyme total, jointly fitted background/gain and halfpoints, and a variance estimated across independent complete titration sets. The 75 main and twelve local-effect configurations give 8,700 Monte Carlo trials. Every plotted interval quantifies Monte Carlo sampling uncertainty.

At the tested reference boundary, rejection counts out of 100 are 6 for ideal calibration, 3 with stock/enzyme errors, 96 with systematic dose bias and 0 with independently corrected bias. The upper Wilson limit for the last count is 0.037. Small-effect rejection varies strongly across the stated K_G range. These are declared noise scenarios, not measured assay precision or a universal type-I error guarantee.

A hidden two-affinity-class generator produces 100/100 single-class reference rejections with mean residual about 0.0279, close to correctly specified scenarios. This is a model-misspecification demonstration: rejection cannot identify a specific change in the single-class regulatory coefficient. Concentration controls correct the tested bias but do not solve hidden heterogeneity.

The text also specifies why GlnD occupancy resolution is a separate, unvalidated assay-development requirement; no successful GlnD native-MS experiment is claimed.

## 4. Show when the bound remains conservative

Fig3A–C compares the general and count-conditioned bounds, the observation/unit dependence and actual differences from 120 paired finite-reaction realizations. The chosen numerical differences can remain many orders of magnitude below both bounds. This does not invalidate a uniform guarantee, but prevents interpreting the bound as a predicted effect size. The conditional refinement is proved in A11 and is not applied unchanged to the productive EPU extension in Fig3D.

## 5. Adjust the paired-paralog role and synchronize the source

Fig5 retains the observed peptide fractions, deterministic synthetic intervals and normalization-dependent upper limits. It is explicitly a prospective application. Its general calibration and mixture results remain in A13–A14, while the main text is shortened. No calibration-transfer validation or direct native occupancy estimate is added by inference from unrelated source measurements.

The revised source, figures, results, manuscript and verification entry point were synchronized as v1.3.2 (comparison-clarification patch following v1.3.1). Version 1.3.3 retains those scientific corrections and completes submission preparation with a shorter initial-submission cover letter. The current S1 Data carries a per-file manifest. Earlier releases are not represented as containing the new work. A new Zenodo DOI is not asserted or substituted until its record is verified; the author manages Zenodo.

## Remaining limits

The revision adds a conditional methodological decision and observation-level performance evidence. It does not resolve native biochemical validation. Applying the certificate requires matched measurements of the count gate and coefficient bounds, plus stationarity at the lower enzyme total. Applying the binding design requires readout/concentration calibration and tests for heterogeneous states. The paired-paralog application still needs matched calibration transfer and an independent bound-fraction endpoint. These limits should remain visible to editors and reviewers.

The subsequent pre-submission review corrections are recorded in `validation/MINOR_REVISION_V1_3_2.md`. The GlnD RMSE change follows joint relaxation of fixed-C, intermediate-state and transported UR-half-range constraints; it is not the isolated effect of variable state capacities.

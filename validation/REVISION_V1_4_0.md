# v1.4.0 evidence and scope audit

## Completed research changes

- Concentration claims now state sufficient conditions for a future matched assay. The nitrogen count gate, coefficient envelope, stationarity, topology and the conflicting source units remain explicit; no native prediction-error guarantee is claimed.
- Expanded raw-assay calibration includes 3,000 calibration trials per null cell and 2,000 independent evaluation trials per cell. Both signal fits undergo nonlinear halfpoint estimation. Finite-grid calibration uncertainty, optimizer failures, active bounds, continuous mixtures and complete nuisance grids are supplied.
- A reusable exact graph/cofactor algorithm tests enzyme cancellation and constructs Bernstein/linear-program enclosures when it fails. Four independent graphs, a rank-one control and an infeasible-observation control pass.
- The independently implemented AceK–IDH application retains singly and doubly bound dimers, full conservation and finite kinetics. Its conditional count gate narrows the example interval width from 0.056995 to 0.036257. The gate is prospective.
- Six published WT/mutant kinetic summaries calibrate finite binding realizations; three withheld ATPase halfpoints test an additional biochemical implication. The disagreement constrains adoption of a minimal mechanism. Published cellular activity data are retained with their distinct units and assay conditions, rather than treated as a matched reconstitution assay.
- Prediction profile likelihood and model-discrimination design are discussed as complementary prior methods. No blanket superiority or claim that prior methods require identified parameters is made.
- The original paired-paralog Fig. 5 and its interpretation are preserved as S7 Fig and S3 Appendix E7. The current Fig. 5 presents AceK. S8 preserves raw-titration/GlnD-profile diagnostics; S9 exposes all null and small-effect nuisance cells.

## Verification

The current integrated check completed 13 stages with zero failures on the designated computation server. Its report is `integrated_checks_v1.4.0.json`; dependencies are in `computation_environment_v1.4.0.json`. The check independently recounts all 840,000 stored evaluation decisions, replays selected complete trials and checks 50 independent raw-curve refits. Fresh AceK, exact-cofactor, conservation and initial-value checks are included. This verification run did not regenerate the entire 8,016,000-fit Monte Carlo study a second time. Generation records retain their actual parameters and scope.

Computational verification is not matched biological validation. Source-summary reconstructions are labeled as such, and uncertainty is not invented where the original publication does not report exact errors or covariance.

## Archive policy

v1.4.0 is a new source snapshot. The old v1.3.4 tag and Zenodo record are not rewritten. Current source metadata, manuscript and S1 Data consistently cite v1.4.0 and the explicitly version-independent concept DOI. The final external submission manifest identifies the newly verified version DOI and hashes. Keeping that post-publication record external avoids falsely claiming an archived commit contains changes added after its release.

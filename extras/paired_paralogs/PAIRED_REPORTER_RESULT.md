# Conditional recovery of GlnK sequestration from paired paralogs

## Biological setting and assumptions

GlnB and GlnK share GlnD. In vitro isolation and mass spectrometry support binding of fully deuridylylated GlnK to AmtB (Durand and Merrick, 2006; see the sources below). We consider count-state trimers with independent modification sites in the free pool, no modification of bound GlnK, and a constant relative UT/UR specificity ratio between paralogs. GlnB is an unsequestered reporter in this construction. These are model conditions, not established properties of every in vivo assay.

Let q_B be the free GlnB modified-site fraction, q_K the total GlnK modified-site fraction, s the fraction of GlnK trimers bound to AmtB, and rho the ratio of the free GlnK and GlnB modification odds. At fixed inputs and steady state,

    F_rho(q_B) = rho q_B / (1 - q_B + rho q_B)
    q_K = (1 - s) F_rho(q_B)
    s = 1 - q_K / F_rho(q_B).

The common GlnD activity cancels, including the constructed regulatory-state transformations. Thus a second paralog can act as an internal reporter for the free-pool driving ratio. The result does not require knowing the glutamine switch law or the common enzyme concentration.

## What must be calibrated

The two means alone do not identify s when rho is unconstrained. A control without AmtB binding gives rho = odds(q_K_control)/odds(q_B_control), provided relative catalytic specificities and free-site ladders are unchanged. This may be a defined reconstitution rather than a cellular knockout, because deleting amtB can also alter input metabolism. With a bounded rho and interval observations, monotonicity propagates uncertainty to an s interval. With completely unknown rho the upper envelope collapses to s <= 1 - q_K, so the extra reporter supplies no unique sequestration estimate.

A time series is not automatically a set of stationary observations. At changing protein totals, synthesis, loss and modification derivatives contribute to the mean balance. The count-state kinetic model in native_trimers.py preserves positivity and explicitly accounts for an assumed unmodified synthesis input. Its initial fits do not establish that the 2017 observations are stationary or that the chosen input/calibration model explains the knockouts.

## Verification and independent-data role

`paired_reporter.py` independently solves 2,000 five-state generators, verifies reconstruction of known bound fractions, and checks calibration at a different input. Maximum error is below 3e-14. This proves the implementation of the specified construction, not physiological validity or priority over all earlier theory.

The 2017 data provide measured paired GlnB/GlnK quantities. Their ratio-derived fractions can be mapped across rho scenarios, but cannot yet be reported as measured sequestration. At late WT starvation both paralogs approach high modification; following upshift their mean trajectories separate. How much of this difference represents binding depends on relative specificity, dynamics and peptide calibration. The public workbook reports SE and n but not replicate vectors or cross-assay covariance.

The mathematical identity is elementary. The prospective contribution is using a biologically available second paralog and a binding-disabled control to remove the unknown GlnD switch from sequestration inference, with explicit failure conditions and observation error. Additional literature comparison is required before claiming novelty for this assay construction.

## Primary sources

- Gosztolai et al. (2017), https://doi.org/10.1016/j.bpj.2017.04.012; measured data and original model https://doi.org/10.6084/m9.figshare.4880003.
- In vitro analysis of the E. coli AmtB-GlnK complex (Durand and Merrick, 2006), https://doi.org/10.1074/jbc.M602477200; the article reports fully deuridylylated GlnK in the isolated complex.
- Control of AmtB-GlnK complex formation by ATP, ADP and 2-oxoglutarate (2010), https://pmc.ncbi.nlm.nih.gov/articles/PMC2945594/; supports a 1:1 trimeric complex and additional ligand controls. The present model does not identify those controls from uridylylation alone.

## What the independent measurements can bound without stationary kinetics

If only fully unmodified GlnK binds AmtB, the bound fraction is at most the unmodified-trimer fraction p0, and p0 <= 1 - qK for any count-state distribution. No site independence or stationarity is needed. At WT -15 min, paired-peptide normalization gives qK = 0.9675734 and a point ceiling 0.0324266. The disclosed component-mean t-interval construction gives a ceiling 0.1615882; separate-total normalization gives 0.6680689. Both observation models and their uncertainty are retained in main Fig5D. These are conditional upper limits, not measured occupancy. The biological message is that calibration determines how strongly the existing high-modification observation constrains sequestration. The sharper paired-reporter reconstruction requires the additional calibration and stationary assumptions above.

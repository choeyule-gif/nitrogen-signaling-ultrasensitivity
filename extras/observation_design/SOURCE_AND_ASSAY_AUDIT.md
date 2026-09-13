# Source and assay constraints for the revised benchmark

## Primary sources inspected

Local full texts: Jiang et al. 1998, doi:10.1021/bi980667m; Jiang et al. 2007, doi:10.1021/bi0620510; Jiang and Ninfa 2011, doi:10.1021/bi201410x. The 2007 public supplement was independently retrieved from ACS Figshare, doi:10.1021/bi0620510.s001, file bi0620510_si_001.pdf. Its metadata and source hash are retained locally. The s002 file has a different figure sequence and lacks the S10 fit used here; it is not substituted for s001.

## Kinetic transfer is not automatic

The 1998 Mg-UT assay uses 25 mM MgCl2 and the Mg-UR assay 50 mM; the 2011 reconstituted system uses 10 mM. Their pH and temperature agree (7.5,30 C), but Tris concentration, nucleotide and effector conditions also differ. The 1998 Mn reactions omit regulators present in the Mg system. Initial-rate apparent kcat/Km ratios are therefore not automatically steady-state microscopic coefficients in the later coupled system. The strict-transfer fit is a stress hypothesis, not the default biochemical constraint.

The 2011 source verifies stationary values with multiple time points. This supports its reconstituted steady-state data, not stationarity of the independent 2017 cellular time series. In the current benchmark the actual observation being bounded must be checked in its own conditions.

## A source-anchored GlnE example

The 2007 supplement S10 uses GlnE .04 uM, GS 2.5 uM, ATP .5 mM and 2-OG .05 mM. Its fit gives alpha_1=.17 and EP/EG/EGP activities 7.5/175/1315. The main article states K_G=15.6 mM and K_P=.18 uM, with a stated sixfold reduction from the apparent direct scale. The supplement prints .18 mM for K_P, inconsistent with that statement. The benchmark follows the main-text unit and explicitly tests the literal alternative. Neither interpretation is silently chosen to maximize performance.

For the observation benchmark, E=.04 uM uses the source enzyme scale. K_G in [7.8,15.6] spans the main-text direct-scale comparison and joint-fit example; it is not a measured confidence interval. Intermediate points are sensitivity calculations. GlnE is held at the source reference parameter family in the finite certificate checks; arbitrary sigma changes are not claimed to preserve the original absolute-rate measurements.

## Constraints change a GlnD assay's contrast

Transporting both the unliganded specificity ratio 38.9012 and UR fold 6.75248 while allowing nonconstant state catalytic capacities can reduce the strict-transfer RMSE to about .0409. This result reaches an extreme K1 boundary (~148 mM); it is not evidence of a plausible binding affinity. A profile with K1=.3,1,10,100 mM gives RMSEs ~.0493,.0430,.0411,.0409, while maximum midpoint singly-liganded separation falls from .242 to .0954,.0108,.00109. The equivalence family persists, but its proposed occupancy contrast is not invariant to biochemical constraints. This directly prevents extrapolating the original .244 gap or five-replicate calculation to every biochemically restricted candidate.

## What an experimental observable would require

The GlnE benchmark assumes a solution readout proportional to the fraction of GlnE bound by PII, with free and saturating-ligand controls. A labeled-GlnE binding titration could target that quantity; label effects, glutamine-dependent signal changes and concentration standards must be checked. The simulated error levels describe a sensitivity benchmark, not established instrument performance. Every independent experiment contains all three glutamine conditions, one shared zero-glutamine baseline and two technical readings at eleven ligand doses. Eight experiments means eight complete three-condition titration sets.

GlnD singly liganded occupancy requires separating zero-, one- and two-glutamine enzyme states, not merely measuring mean binding. Native intact-protein MS is a candidate pilot: compare apo and +1/+2 glutamine masses, validate assignments with ligand titration and ACT controls, and test ligand loss/adducts across source conditions. No GlnD precision or successful state resolution has been demonstrated here; GlnK native-MS success does not establish it. The revised central design therefore does not rest on an assumed GlnD state-resolving assay.

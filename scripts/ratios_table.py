"""Table 'tab:ratios': the ratio T_E/T_L that eq (sandwich) controls, for each experiment.
First cycle: 0.05 uM UTase/UR against 0.55 mM glutamine at the half-point of the titration
reanalysed here (9e-5); the largest value across the first-cycle titrations of Jiang & Ninfa
2011 is 3.6e-3 (highest enzyme, lowest half-point).  Second cycle: 0.1 uM ATase against PII
at 36, 5 and 0.5 uM (0.003, 0.02, 0.2).  Also eq (sandwich) at h = 2 for the refit: 2e-4."""
import _path  # noqa: F401
from esbm import save_json
from esbm.data import constants
c = constants()["jiangninfa2011"]
first = c["UTase_uM_titration"]/(c["first_cycle_halfpoint_mM"]*1000)
second = [c["ATase_uM"]/p for p in c["cascade_PII_uM"]]
# the h = 2 refit: the relative gap between free and total effector is at most T_E/T_L at the half-point
refit = 2*c["UTase_uM_titration"]/(0.5595*1000)          # h = 2: each liganded enzyme holds two effector molecules
print("first cycle, titration reanalysed: T_E/T_L = %.1e" % first)
print("first cycle, refit half-point 0.5595 mM, h = 2 (factor h): %.1e" % refit)
print("second cycle at PII 36, 5, 0.5 uM: " + ", ".join("%.3g" % x for x in second))
save_json("ratios_table", dict(first_cycle=first, refit=refit, second_cycle=second))

#!/usr/bin/env python3
"""Reproduce every number and figure of the manuscript, in dependency order.

    python run_all.py            # everything (about 40 minutes on an Apple-silicon laptop)
    python run_all.py --quick    # skip the four long searches (adversarial, counterexamples, MCMC, profile)

Each script prints the quantities it is responsible for and writes results/<name>.json;
tests/test_reproduce.py then compares those files with the values quoted in the manuscript.
"""
import os, sys, subprocess, time
ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ROOT, "scripts")
QUICK = "--quick" in sys.argv
LONG = {"mono_adversarial.py", "counterexamples.py", "mcmc_h.py", "fit_titration.py"}
ORDER = [
    # verifications of the theory (independent of the data)
    "verify_symbolic.py", "verify_adair.py", "verify_nu_bound.py", "verify_deadend_60digit.py",
    "sequestration_sweeps.py", "deficiency_acr.py", "counterexamples.py", "mono_adversarial.py",
    # reanalysis of the published data
    "fit_titration.py", "mcmc_h.py", "cascade_budget.py", "sensitivity_major1.py", "nu2_major2.py",
    "mixture_bias.py", "conversion_factor.py", "swing_profile.py", "direct_input.py", "ratios_table.py",
    # figures
    "fig1_invariance.py", "fig2_sites.py", "fig3_plateau.py", "fig4_plane.py",
    "fig5_integer.py", "fig6_uniqueness.py", "fig7_tests.py", "fig8_budget.py",
]
failed = []
for s in ORDER:
    if QUICK and s in LONG:
        print("skip %s (--quick)" % s); continue
    t0 = time.time(); print("=" * 78); print("running", s, flush=True)
    log = open(os.path.join(ROOT, "results", s.replace(".py", ".log")), "w")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, s)], stdout=log, stderr=subprocess.STDOUT, cwd=SCRIPTS)
    log.close()
    print("   %s in %.0f s -> results/%s" % ("ok" if r.returncode == 0 else "FAILED", time.time()-t0, s.replace(".py", ".log")), flush=True)
    if r.returncode: failed.append(s)
print("=" * 78)
print("failed:", failed if failed else "none")
sys.exit(1 if failed else 0)

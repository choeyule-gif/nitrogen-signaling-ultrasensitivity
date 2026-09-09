"""Compare results/*.json (written by the scripts) with the values quoted in the manuscript.

Run `python run_all.py` (or `make all`) first, then `python -m pytest -q tests`.
Each assertion names the manuscript statement it checks.  Tolerances are those of the rounding
used in the text; the long stochastic searches are checked on their qualitative conclusion.
"""
import os, json
try:
    import pytest
except ImportError:           # pytest is optional: python tests/test_reproduce.py runs the same checks
    class _Skip(Exception): pass
    class pytest:  # noqa: N801
        @staticmethod
        def skip(msg): raise _Skip(msg)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
def R(name):
    p = os.path.join(ROOT, "results", name + ".json")
    if not os.path.exists(p): pytest.skip("results/%s.json missing: run run_all.py" % name)
    return json.load(open(p))
def close(x, y, rel=0.0, abs_=0.0): return abs(x-y) <= abs_ + rel*abs(y)

# ---------------------------------------------------------------- Fig 1 (eq. chain, eq. nogo)
def test_fig1_enzyme_cancels():
    d = R("fig1_invariance")
    assert d["enzyme_fold_range"] > 5000                       # "over a 5,000-fold range"
    assert d["max_abs_dev_from_Ki_over_l"] < 1e-13             # "agree with K_i/l to 5e-14"
    assert d["deadend_ratio_drift"] < 1e-9                     # "drift by 5e-10 at T_E/T_L = 1.1e-9"
def test_deadend_60_digit():
    d = R("verify_deadend_60digit")
    assert d["1e-9"]["chain_rel_dev"] < 1e-28                  # "1.3e-29 in 60-digit arithmetic"
    assert close(d["1e-4"]["ratio_drift"], 5.4e-5, rel=0.05)   # "5.4e-5 at T_E/T_L = 1.1e-4"
    assert close(d["1e-9"]["ratio_drift"], 5.4e-10, rel=0.05)  # "5.4e-10 at 1.1e-9, in exact proportion"
# ---------------------------------------------------------------- Fig 2, eq. hill1, eq. zeroorder
def test_fig2_sites():
    d = R("fig2_sites")
    assert d["max_abs_nH_minus_1"] < 1e-7                      # "equal to 1 to 4e-8"
    assert close(d["nH_TE_over_TS_0p1"], 0.988, abs_=0.001)    # "0.988 rather than 1"
    assert close(d["nH_TE_over_TS_0p2"], 0.976, abs_=0.001)    # "doubling the enzyme moves it to 0.976"
def test_plateau():
    d = R("fig3_plateau")
    assert close(d["interior_min"], 0.99745, abs_=2e-5)        # "between 0.99745 and 1"
    assert close(d["interior_width"], 2.5e-3, rel=0.05)        # "a plateau of width 2.5e-3"
def test_sequestration_sweeps():
    d = R("sequestration_sweeps")
    assert close(d["target_sweep"]["nH_min"], 0.981, abs_=0.001)   # "leaves n_H between 0.981 and 0.9996"
    assert close(d["target_sweep"]["nH_max"], 0.9996, abs_=0.0002)
# ---------------------------------------------------------------- eq. adair, eq. factor, eq. nubound, eq. integer
def test_adair_and_twoform():
    d = R("verify_adair")
    assert d["adair_sets"] == 1500 and d["adair_worst_rel_dev"] < 1e-7        # eq. adair over 1,500 sets (double precision, conditioning-limited)
    assert d["adair_40digit_worst"] < 1e-20 and d["twoform_40digit_drift"] < 1e-30   # the same sets with the kernel in 40-digit arithmetic
    assert d["twoform_drift_1000fold"] < 1e-8                                 # enzyme independence of the two-form ratios (see README note)
    for h in ("1", "2", "3"):
        assert all(close(v, float(h), abs_=1e-8) for v in d["integer_hill"][h])   # eq. integer: n_H = h exactly
    assert close(d["nu_sup_n2_K1_1e_minus3_K2_1e3"], 1.999, abs_=5e-4)        # "nu rises to 1.999"
def test_nu_bound():
    d = R("verify_nu_bound")
    assert d["exceedances"] == 0                                              # "no evaluation of nu exceeded n"
    assert d["gap_to_n"]["n2_eps1e-20"] < 1e-19                               # "approached to 1e-20"
    assert d["binomial_max_dev_overall"] < 6e-43                              # "below 6e-43"
# ---------------------------------------------------------------- eq. quartic, eq. cond, symbolic identities
def test_symbolic():
    d = R("verify_symbolic")
    assert d["A_chain_free_of_e"] and d["B"]["mobius_ok"] and d["B"]["degree"] == 4
    assert d["B"]["c4_ok"] and d["B"]["c0_ok"] and d["B"]["quartic_form_ok"]  # eq. quarticco, eq. quartic
    assert d["C_hillid_residual_zero"] and d["D_twoform_chain_ok"]
    assert d["E_twoform_l2_degrees"] == [1, 2]                                # "linear to a quadratic polynomial in e"
    assert d["F_omcov_violations"] == 0
    assert d["G_detJ_signs"]["1"]["negative"] == 0 and d["G_detJ_signs"]["1"]["positive"] == 18
    assert d["G_detJ_signs"]["2"]["positive"] == 44 and d["G_detJ_signs"]["2"]["negative"] == 7   # "44 positive against 7 negative at n=2"
    assert d["H"]["free_of_e"] and d["H"]["carries_s0"] and d["H"]["reduces_to_chain"]
def test_counterexamples():
    d = R("counterexamples")
    assert d["cond_i"]["n_roots"] == 3 and d["cond_i"]["all_positive"]        # Fig 6b: three positive steady states
    assert close(d["cond_i"]["spread"], 2.89, abs_=0.01)                      # "a spread of 2.89 in the catalytic ratios"
    assert len(d["cond_ii"]) >= 1 and all(x["n_roots"] == 3 and x["all_positive"] and x["max_elasticity"] > 1 for x in d["cond_ii"])
def test_deficiency():
    d = R("deficiency_acr")
    assert d["n1"]["deficiency"] == 1 and d["n1"]["class_deficiencies"] == [0, 0, 0] and d["n1"]["sf_pairs"] == 0
    for n in range(1, 7):
        g = d["general"][str(n)]; assert g["complexes"] == 4*n+4 and g["linkage_classes"] == 3 and g["rank"] == 3*n+1 and g["deficiency"] == n
    b = d["buffering"]
    assert close(b["e_x1"], 0.1128, abs_=1e-4) and close(b["e_x16"], 0.0326, abs_=1e-4)     # "e 0.1128 -> 0.0326"
    assert close(b["s0_x1"], 0.2669, abs_=1e-4) and close(b["s0_x16"], 4.7647, abs_=1e-3)   # "s0 0.2669 -> 4.7647"
def test_adversarial():
    d = R("mono_adversarial")
    assert d["n_pairs"] == 21 and d["violations"] == 0                        # "no violation ... 21 pairs (n,h)"
    assert 5e-15 < d["verified_min"] < 2e-14 and d["verified_max"] < 2e-10    # "between 9e-15 and 1e-10"
# ---------------------------------------------------------------- R2: the titration (Methods, tab:ic)
def test_fit_titration():
    d = R("fit_titration"); t = d["table"]
    assert close(t["h=2"]["SSE"], 0.0149, abs_=1e-4) and close(t["h=1"]["SSE"], 0.3658, abs_=1e-4)
    assert close(d["dAIC_h1_vs_h2"], 38.38, abs_=0.05) and close(d["dAIC_h3_vs_h2"], 24.93, abs_=0.05)
    assert close(d["dAIC_h2_vs_free"], 1.95, abs_=0.02)
    assert close(d["h_free"], 2.018, abs_=0.002)                              # "n_H^fit = 2.02"
    assert close(d["shapiro_W"], 0.964, abs_=0.001) and close(d["shapiro_p"], 0.83, abs_=0.01)
    assert close(d["profile"]["lo"], 1.86, abs_=0.011) and close(d["profile"]["hi"], 2.20, abs_=0.011)
    assert close(d["profile"]["stat_h1"], 38.4, abs_=0.1) and close(d["profile"]["stat_h3"], 25.0, abs_=0.1)
    r = d["standardised_residuals"]
    assert close(r["1"]["max_abs"], 6.9, abs_=0.05) and close(r["3"]["max_abs"], 4.7, abs_=0.05) and r["2"]["max_abs"] < 1.9
    assert close(d["sigma_hat_h2"], 0.0407, abs_=1e-4)
    assert close(d["floorless"]["two_dL"], 31.9, abs_=0.1) and close(d["floorless"]["h"], 1.42, abs_=0.01)
    assert close(d["refractory"]["q"], 0.15, abs_=0.005) and close(d["refractory"]["h"], 1.96, abs_=0.01) and close(d["refractory"]["two_dL"], 1.30, abs_=0.02)
    e = d["endpoints"]; assert close(e["drop_highest"], 2.06, abs_=0.006) and close(e["drop_lowest"], 2.03, abs_=0.006) and close(e["drop_last_two"], 2.04, abs_=0.006)
    assert close(d["free_fit"]["Umin"], 0.45, abs_=0.006) and close(d["free_fit"]["se"][1], 0.04, abs_=0.005)
    assert close(d["hillid"]["eta_at_theta_half"], 1.63, abs_=0.006) and close(d["hillid"]["eta_at_S"], 1.71, abs_=0.006)
def test_mcmc():
    d = R("mcmc_h")
    assert close(d["acceptance_mean"], 0.245, abs_=0.005) and d["R_hat"] < 1.001 and d["ESS"] > 1e4
    assert close(d["median"], 2.025, abs_=0.003) and close(d["ci95"][0], 1.81, abs_=0.006) and close(d["ci95"][1], 2.28, abs_=0.006)
    assert d["P_h_lt_1p5"] < 1e-4 and d["P_h_lt_1"] == 0
# ---------------------------------------------------------------- R1, R3, R4, R5 and the referee analyses
def test_mixture_bias():
    d = R("mixture_bias")
    assert close(d["di_tri_trimer_percent"], 15.0, abs_=0.1) and close(d["di_tri_subunit_percent"], 27.0, abs_=0.6)
    row = [r for r in d["rows"] if abs(r["nu3"]-1.40) < 0.01][0]; assert close(row["bias"], 1.028, abs_=0.001)   # "raises the coefficient by 2.8%"
def test_conversion_factor():
    d = R("conversion_factor")
    lo = min(d["fit_range"][0], d["range_range"][0]); hi = max(d["fit_range"][1], d["range_range"][1])
    assert 1.02 <= lo and hi <= 1.06                                          # the three named families: "1.04 +- 0.02, an envelope over both"
    assert d["locus"]["fit"][2] > 10 and d["locus"]["range"][2] > 10          # the full locus is computed and reported (wider; see README)
    vals = list(d["h_independence"].values()); assert max(vals)-min(vals) < 1e-6   # "exactly independent of h"
def test_swing():
    d = R("swing_profile")
    assert close(d["theta0_free"], 0.989, abs_=0.001) and close(d["theta_inf_free"], 0.151, abs_=0.001)
    assert close(d["w0_over_Khat"], 0.011, abs_=0.0006) and 400 < d["swing_ML"] < 600     # "about 5e2"
    assert close(d["lower95"], 175, rel=0.02) and not d["closed_above"]        # "[175, inf) at 95%"
    assert close(d["stat"]["6.8"], 54, abs_=1) and close(d["min_inhibition_fold"], 26, abs_=0.5)
def test_cascade_budget():
    d = R("cascade_budget")
    def near(xs, ys, tol=0.006): return all(close(x, y, abs_=tol) for x, y in zip(xs, ys))
    assert near(d["lin"]["readout"], [2.44, 1.87, 1.76]) and near(d["lin"]["composite"], [4.915, 3.78, 3.55])   # text: 4.92 (computed 4.9145)
    assert near(d["compl"]["readout"], [2.32, 1.67, 1.57]) and near(d["compl"]["composite"], [4.67, 3.37, 3.17])
    assert close(d["lin"]["factor_sup"], 2.4399, abs_=2e-4) and close(d["matched"]["factor"], 1.25935, abs_=1e-5)
    assert close(d["matched"]["predicted_halfpoint_mM"], 1.08, abs_=0.005)
    assert min(d["compl"]["factor_min_24dec"], d["lin"]["factor_min_24dec"]) < 1.16   # "runs from 1.15"
    assert len(set(round(v, 5) for v in d["sup_vs_h"].values())) == 1         # "2.43990 to five decimals for h between 1 and 4"
    b = d["bootstrap"]; assert b["n_usable"] == 536 and b["peak_middle"] == 536
    for p, med, lo, hi in (("36", 0.92, 0.83, 1.03), ("5", 1.70, 1.57, 1.84), ("0.5", 1.47, 1.35, 1.58)):
        r = b["PII_"+p]["residual"]; assert close(r[0], med, abs_=0.006) and close(r[1], lo, abs_=0.006) and close(r[2], hi, abs_=0.006)
def test_major1():
    d = R("sensitivity_major1")
    assert close(d["drop_percent"], 13.9, abs_=0.1) and close(d["middle_to_vanish"], 5.57, abs_=0.01)
    assert close(d["P_peak_vs_sigma"]["0.05"][0], 0.98, abs_=0.01) and close(d["P_peak_vs_sigma"]["0.1"][0], 0.85, abs_=0.01) and close(d["P_peak_vs_sigma"]["0.15"][0], 0.76, abs_=0.01)
def test_major2():
    d = R("nu2_major2")
    r134 = [r for r in d["rows"] if abs(r["nu2"]-1.34) < 0.01][0]; assert all(close(x, y, abs_=0.006) for x, y in zip(r134["residual"], [0.735, 1.38, 1.20]))
    r223 = [r for r in d["rows"] if abs(r["nu2"]-2.23) < 0.01][0]; assert all(close(x, y, abs_=0.006) for x, y in zip(r223["residual"], [0.51, 0.955, 0.83]))
    assert d["peak_always"]
def test_ratios():
    d = R("ratios_table")
    assert close(d["first_cycle"], 9e-5, rel=0.03) and close(d["refit"], 2e-4, rel=0.12)
    assert [round(x, 3) for x in d["second_cycle"]] == [0.003, 0.02, 0.2]
def test_fig7():
    d = R("fig7_tests")
    assert all(close(v, 1.000, abs_=5e-4) for v in d["b_ours"])              # "1.000 in every condition"
    assert close(d["c_bridge"][0], 1.08, abs_=0.006) and close(d["c_bridge"][2], 1.26, abs_=0.006)   # "1.08 against 1.26"
def test_fig4():
    d = R("fig4_plane")
    assert close(d["nu_eta_product"], 2.02, abs_=0.001)

if __name__ == "__main__":            # plain runner, so that pytest is optional
    import sys, traceback
    fails = 0; skips = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print("PASS %s" % name)
            except Exception as ex:      # pytest.skip raises an Exception subclass as well
                if type(ex).__name__.lower().startswith("skip"): skips += 1; print("SKIP %s: %s" % (name, ex))
                else: fails += 1; print("FAIL %s: %s" % (name, traceback.format_exc().strip().splitlines()[-1]))
    print("%d failed, %d skipped" % (fails, skips)); sys.exit(1 if fails else 0)

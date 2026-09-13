"""Numbers behind the five figures of the revised manuscript, written as CSV for the MATLAB
figure scripts in matlab/ (nothing is recomputed at plot time), plus revision_numbers.json
with every quantity the revised text quotes.

The analysis follows the revised manuscript's conventions: a bounded four-parameter Hill fit
(0 <= U_min <= 1.5, 1.5 <= U_max <= 3, -8 <= log S <= 5, 0.1 <= h <= 8) to the twelve digitised
points of Jiang & Ninfa (2011) Fig 2B; fixed-exponent fits under the same bounds; corrected
AIC with the error variance counted; a profile likelihood in h and in the regulatory swing
W = f_max/w_0; 2,000 Gaussian coordinate perturbations (log-input sd 0.06, ordinate sd 0.06,
seed 20260909); and the cascade composite Y = k p_0/(k p_0 + theta) scored by the range
coefficient after calibrating k to each reported cascade midpoint.
"""
import os, sys, json, csv
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.special import gammaln, logsumexp
from esbm.data import digitised_titration, constants
from esbm import RESULTS, ROOT
FD = os.path.join(ROOT, "figdata")
os.makedirs(FD, exist_ok=True)
x, y = digitised_titration(); N = len(x)
cc = constants()["jiangninfa2011"]
BOUNDS = ([0, 1.5, -8, .1], [1.5, 3, 5, 8])
def csvw(name, cols, comment):
    keys = list(cols); rows = zip(*[np.asarray(cols[k]).ravel() for k in keys])
    with open(os.path.join(FD, name), "w", newline="") as f:
        f.write("# " + comment + "\n"); w = csv.writer(f); w.writerow(keys); w.writerows(rows)
def hill(xx, p):
    lo, hi, ls, h = p
    return lo + (hi-lo)/(1+(np.asarray(xx)/np.exp(ls))**h)
def fit(xx, yy, h=None, start=None):
    p0 = np.array([.45, 2.97, np.log(.56), 2.]) if start is None else np.array(start)
    if h is None:
        r = least_squares(lambda p: hill(xx, p)-yy, p0, bounds=BOUNDS, xtol=1e-10, ftol=1e-10, gtol=1e-10)
        return r.x, float(r.fun@r.fun)
    r = least_squares(lambda q: hill(xx, [*q, h])-yy, p0[:3], bounds=(np.array(BOUNDS[0])[:3], np.array(BOUNDS[1])[:3]), xtol=1e-10, ftol=1e-10, gtol=1e-10)
    return np.r_[r.x, h], float(r.fun@r.fun)
out = {}
p, sse = fit(x, y); S = np.exp(p[2]); tlo, thi = p[:2]/3
out['fit'] = dict(Umin=p[0], Umax=p[1], S_mM=S, h=p[3], SSE=sse, theta0=thi, theta_inf=tlo)   # theta0 is the low-glutamine plateau
# NB: p[0] = U_min (high-glutamine floor) -> theta_inf ; p[1] = U_max -> theta_0
out['fit']['theta0'] = p[1]/3; out['fit']['theta_inf'] = p[0]/3
print("fit: Umin %.6f Umax %.6f S %.6f h %.6f SSE %.8f" % (p[0], p[1], S, p[3], sse))
# ---- model comparison (bounded)
fits = {}; rows = []
for h in [1, 2, 3, None]:
    pp, ss = fit(x, y, h); k = 5 if h is None else 4
    aic = N*np.log(ss/N)+2*k; aicc = aic+2*k*(k+1)/(N-k-1)
    fits[str(h)] = (pp, ss); rows.append([0 if h is None else h, ss, k, aic, aicc])
rows = np.array(rows); rows[:, 3] -= rows[:, 3].min(); rows[:, 4] -= rows[:, 4].min()
out['model_comparison'] = {("h=%d" % r[0] if r[0] else "h free"): dict(SSE=r[1], k=int(r[2]), dAIC=r[3], dAICc=r[4]) for r in rows}
r1 = least_squares(lambda q: q[0]+(q[1]-q[0])/(1+x/np.exp(q[2]))-y, [.4, 3.3, -.5], xtol=1e-12, ftol=1e-12, gtol=1e-12)
out['h1_unbounded'] = dict(Umax=r1.x[1], SSE=float(r1.fun@r1.fun))
for k_, v in out['model_comparison'].items(): print("  %-7s SSE %.6f dAIC %.3f dAICc %.3f" % (k_, v['SSE'], v['dAIC'], v['dAICc']))
# ---- profile in h
def hprof(h): return N*np.log(fit(x, y, h)[1]/sse)
hci = [brentq(lambda h: hprof(h)-3.841458820694, 1, p[3]), brentq(lambda h: hprof(h)-3.841458820694, p[3], 3)]
out['profile_h'] = hci; print("profile h interval [%.3f, %.3f]" % tuple(hci))
# ---- swing profile
def swingfit(w):
    def fun(q):
        low, ls, h = q; t = low/3; hi = 3*w*t/(1-t+w*t)
        return hill(x, [low, hi, ls, h])-y
    r = least_squares(fun, [p[0], p[2], p[3]], bounds=([1e-8, -8, .1], [1.5, 5, 8]), xtol=1e-10, ftol=1e-10, gtol=1e-10)
    return float(r.fun@r.fun)
swing = (thi/(1-thi))/(tlo/(1-tlo))
wlo = np.exp(brentq(lambda z: N*np.log(swingfit(np.exp(z))/sse)-3.841458820694, np.log(10), np.log(swing)))
Rrev = 7.93/1.17
out['swing'] = dict(point=swing, lower95=wlo, R_rev=Rrev, R_fwd_lower=wlo/Rrev, w0_over_Khat=1/out['fit']['theta0']-1)
print("swing %.2f lower %.2f -> R_fwd >= %.2f" % (swing, wlo, wlo/Rrev))
ww = np.logspace(1.5, 5, 110); loss = np.array([N*np.log(swingfit(w)/sse) for w in ww])
csvw("Fig3D_swing_profile.csv", dict(swing=ww, loss=loss), "Fig 4D: profile likelihood loss N log(SSE_W/SSE_free) against the regulatory swing W = f_max/w_0; threshold 3.841; lower endpoint %.2f" % wlo)
csvw("Fig3D_markers.csv", dict(quantity=["threshold", "lower95", "point_estimate"], value=[3.841458820694, wlo, swing]), "Fig 4D reference values")
# refractory alternative
rr = least_squares(lambda q: hill(x, [q[0], 3, q[1], q[2]])-y, [p[0], p[2], p[3]], bounds=([0, -8, .1], [1.5, 5, 8]))
out['refractory'] = dict(q=rr.x[0]/3, h=rr.x[2], two_dL=N*np.log(float(rr.fun@rr.fun)/sse))
print("refractory: q %.3f h %.3f 2dL %.2f" % (rr.x[0]/3, rr.x[2], out['refractory']['two_dL']))
# ---- cascade composite
def response(l, pp, k, reverse='linear'):
    t = hill(l, pp)/3; fw = (1-t)**3; rv = t if reverse == 'linear' else 1-fw
    return k*fw/(k*fw+rv)
def limits(pp, k, reverse): return response(0, pp, k, reverse), response(1e15, pp, k, reverse)
def calibration(pp, smid, reverse='linear'):
    def fun(logk):
        k = np.exp(logk); a, b = limits(pp, k, reverse)
        return response(smid, pp, k, reverse)-(a+b)/2
    return np.exp(brentq(fun, -28, 28))
def score(pp, k, reverse='linear'):
    tl, th = pp[:2]/3; a, b = limits(pp, k, reverse); ts = []
    for q in [.1, .9]:
        target = a+q*(b-a)
        def fun(t):
            fw = (1-t)**3; rv = t if reverse == 'linear' else 1-fw
            return k*fw/(k*fw+rv)-target
        t = brentq(fun, tl+1e-13, th-1e-13)
        ts.append(pp[2]+np.log((th-t)/(t-tl))/pp[3])
    return np.log(81)/(ts[1]-ts[0])
pii = np.array([.5, 5., 36.]); mid = np.array([.600, .530, .225]); observed = np.array([5.23, 6.46, 4.48])
ks = np.array([calibration(p, m) for m in mid]); pred = np.array([score(p, k) for k in ks]); res = observed/pred
predc = np.array([score(p, calibration(p, m, 'compl'), 'compl') for m in mid])
out['cascade'] = dict(PII_uM=pii.tolist(), midpoint_mM=mid.tolist(), reported=observed.tolist(), kappa=ks.tolist(), calculated=pred.tolist(), ratio=res.tolist(), calculated_complementary=predc.tolist(), ratio_complementary=(observed/predc).tolist())
print("cascade: calculated", np.round(pred, 3), "ratio", np.round(res, 3), "| complementary", np.round(predc, 3))
# ---- 2,000 Gaussian coordinate perturbations, seed 20260909
rng = np.random.default_rng(20260909); B = 2000; ps = []; preds = []; sw = []; failed = 0
for i in range(B):
    xx = x*np.exp(rng.normal(0, .06, N)); yy = y+rng.normal(0, .06, N)
    try:
        pp, ss = fit(xx, yy, start=p)
        vals = np.array([score(pp, calibration(pp, m)) for m in mid])
        ps.append(pp); preds.append(vals)
        l_, h_ = pp[:2]/3; sw.append((h_/(1-h_))/(l_/(1-l_)) if h_ < 1 else np.inf)
    except (ValueError, RuntimeError, FloatingPointError): failed += 1
ps = np.array(ps); preds = np.array(preds); sw = np.array(sw)
q_h = np.quantile(ps[:, 3], [.025, .5, .975]); q_res = np.quantile(observed/preds, [.025, .5, .975], axis=0)
out['perturbations'] = dict(B=B, usable=len(ps), failed=failed, h_quantiles=q_h.tolist(), ratio_quantiles=q_res.tolist(), swing_p025=float(np.quantile(sw, .025)), swing_at_ceiling=int(np.sum(ps[:, 1] > 3-1e-5)))
print("perturbations: h %s ; ratio quantiles %s ; swing p2.5 %.1f ; U_max at ceiling %d" % (np.round(q_h, 3), np.round(q_res, 3).tolist(), np.quantile(sw, .025), out['perturbations']['swing_at_ceiling']))
# Major-1 style sensitivity with these draws: lognormal error on the three reported coefficients
rng2 = np.random.default_rng(1); rdraw = observed/preds; sens = {}
for sig in (0.05, 0.10, 0.15):
    M = 200000; k_ = rng2.integers(0, len(rdraw), M); m = observed*np.exp(sig*rng2.standard_normal((M, 3)))
    r = m/preds[k_]           # ratios at PII 0.5, 5, 36
    sens[str(sig)] = float(np.mean((r[:, 1] > r[:, 0]) & (r[:, 1] > r[:, 2])))
out['peak_sensitivity'] = dict(P_peak_middle=sens, middle_to_vanish=float(pred[1]*max(res[0], res[2])), drop_percent=float(100*(1-pred[1]*max(res[0], res[2])/observed[1])))
print("peak survives with probability", sens, "; 6.46 would have to fall to %.2f (%.1f%%)" % (out['peak_sensitivity']['middle_to_vanish'], out['peak_sensitivity']['drop_percent']))
# ---- local slope of the physical mean fraction
def nlocal(l, pp):
    t = hill(l, pp)/3; q = (l/np.exp(pp[2]))**pp[3]
    return (pp[1]-pp[0])/3*pp[3]*q/(1+q)**2/(t*(1-t))
out['local_slope_at_S'] = float(nlocal(S, p))
# ---- state distributions
n = 12; ii = np.arange(n+1); logbin = gammaln(n+1)-gammaln(ii+1)-gammaln(n-ii+1)
def dist(j):
    z = logbin+j*(ii-n/2)**2; return np.exp(z-logsumexp(z))
J = brentq(lambda j: 4*np.sum(dist(j)*(ii-n/2)**2)/n-2.2, 0, .3); out['nu2_example_J'] = J

# ======================= figure data
import shutil
PANELS = os.path.join(RESULTS, "panels")        # written by the Python figure scripts
for src, dst in (("Fig2a_theta.csv", "Fig1B_theta.csv"), ("Fig2b_local_hill.csv", "Fig1C_local_hill.csv"), ("Fig7_titrations.csv", "Fig5BC_titrations.csv"),
                 ("Fig1a_enzyme_sweep.csv", "FigS1A_enzyme_sweep.csv"), ("Fig1b_deadend_sweep.csv", "FigS1B_deadend_sweep.csv"),
                 ("Fig6a_residual.csv", "FigS5A_residual.csv"), ("Fig6b_residual.csv", "FigS5B_residual.csv"),
                 ("Fig3_plateau_surface.csv", "FigS3_plateau_surface.csv")):
    shutil.copy(os.path.join(PANELS, src), os.path.join(FD, dst))
nu = np.linspace(.55, 3, 300)
csvw("Fig2A_isocurves.csv", dict(nu=nu, **{"eta_c%g" % c: c/nu for c in (1, 2, 3, 4, 6)}), "Fig 3A: curves of constant local coefficient nu*eta = c; reference lines nu = 1 (independent ladder) and nu = 3 (trimer ceiling)")
rat = np.logspace(0, 4, 300)
csvw("Fig2B_eta_max.csv", dict(affinity_ratio=rat, eta_max=2*np.sqrt(rat)/(1+np.sqrt(rat))), "Fig 3B: maximum switch elasticity of eq (6) against A_UT/A_UR")
csvw("Fig2B_markers.csv", dict(ion=["Mg", "Mn"], affinity_ratio=[80/70, 700/150], eta_max=[2*np.sqrt(80/70)/(1+np.sqrt(80/70)), 2*np.sqrt(700/150)/(1+np.sqrt(700/150))]), "Fig 3B markers from the reported kinetic constants")
xx = np.logspace(-2, 1.3, 400)
csvw("Fig2C_local_slope.csv", dict(glutamine_mM=xx, nH_loc=nlocal(xx, p)), "Fig 3C: local logit slope of the fitted physical occupancy; fitted exponent %.4f; value at S = %.4f" % (p[3], out['local_slope_at_S']))
csvw("Fig2C_markers.csv", dict(quantity=["S_mM", "nH_loc_at_S", "nH_fit"], value=[S, out['local_slope_at_S'], p[3]]), "Fig 3C markers")
xx = np.logspace(-2, 1.25, 350); curves = np.array([hill(xx, pp) for pp in ps]); band = np.quantile(curves, [.025, .975], axis=0)
csvw("Fig3A_points.csv", dict(glutamine_mM=x, UMP_per_trimer=y), "Fig 4A: the twelve digitised points (data/jiangninfa2011_fig2B_digitised.csv)")
csvw("Fig3A_fit.csv", dict(glutamine_mM=xx, fit=hill(xx, p), band_lo=band[0], band_hi=band[1]), "Fig 4A: free-exponent fit and pointwise 2.5-97.5 percentiles over 2,000 coordinate perturbations")
csvw("Fig3A_fixed_h_fits.csv", dict(glutamine_mM=xx, **{"fit_h%d" % h: hill(xx, fits[str(h)][0]) for h in (1, 2, 3)}),
     "Fig 4A: the bounded fixed-exponent fits h = 1, 2, 3 as curves, on the grid of Fig3A_fit.csv")
sigma = np.sqrt(fits['2'][1]/(N-3))
csvw("Fig3B_residuals.csv", dict(glutamine_mM=x, **{"r_h%d" % h: (y-hill(x, fits[str(h)][0]))/sigma for h in (1, 2, 3)}), "Fig 4B: residuals of the bounded fixed-exponent fits divided by sigma_hat = %.4f (h = 2 fit, 9 dof)" % sigma)
csvw("Fig3C_exponent_draws.csv", dict(h=ps[:, 3]), "Fig 4C: fitted exponent in each of the 2,000 perturbation draws; quantiles 2.5/50/97.5 = %.4f/%.4f/%.4f" % tuple(q_h))
csvw("Fig3C_markers.csv", dict(quantity=["h_fit", "p025", "p975"], value=[p[3], q_h[0], q_h[2]]), "Fig 4C markers")
xx = np.logspace(-2, 1, 400); cols = dict(glutamine_mM=xx)
for t, k in zip(pii, ks):
    aa, bb = limits(p, k, 'linear'); cols["Y_norm_PII_%g" % t] = (response(xx, p, k)-aa)/(bb-aa)
csvw("Fig4A_composites.csv", cols, "Fig 5A: calibrated composite responses normalised between their plateaus; midpoints at the reported cascade half-points")
csvw("Fig4A_midpoints.csv", dict(PII_uM=pii, midpoint_mM=mid), "Fig 5A midpoints")
csvw("Fig4B_coefficients.csv", dict(PII_uM=pii, calculated=pred, reported=observed), "Fig 5B: calculated range coefficients (linear reading of the reverse input) and reported cascade coefficients")
csvw("Fig4C_ratios.csv", dict(PII_uM=pii, ratio=res, p025=q_res[0], p975=q_res[2]), "Fig 5C: reported / calculated with central 95%% ranges over the 2,000 upstream coordinate perturbations")
xx = np.logspace(-2, 2, 300)
csvw("Fig5A_binding.csv", dict(scaled_glutamine=xx, one_ligand=xx/(1+xx), two_ligands=2*xx**2/(1+xx**2)), "Fig 6A: ligands bound per enzyme for a one-site law and a concerted two-site law")
probs = np.c_[dist(0), dist(J)]
csvw("Fig5D_distributions.csv", dict(modified_sites=ii, independent_nu1=probs[:, 0], example_nu2p2=probs[:, 1]), "Fig 6D: twelve-site distributions with mean 1/2, weights binomial(12,i) exp[J (i-6)^2], J = 0 and J = %.6f" % J)
csvw("S1_perturbation_draws.csv", dict(Umin=ps[:, 0], Umax=ps[:, 1], S_mM=np.exp(ps[:, 2]), h=ps[:, 3], composite_PII_0p5=preds[:, 0], composite_PII_5=preds[:, 1], composite_PII_36=preds[:, 2], regulatory_swing=sw),
     "S1 Data: the 2,000 coordinate-perturbation draws (seed 20260909; log-input sd 0.06, ordinate sd 0.06): refitted parameters, calibrated composite coefficients at PII 0.5, 5, 36 uM, and the regulatory swing")
nu2 = json.load(open(os.path.join(ROOT, "results", "nu2_major2.json")))
csvw("S1_downstream_ladder_factor.csv", dict(cooperativity_w=[r['w'] for r in nu2['rows']], nu2_at_halfpoint=[r['nu2'] for r in nu2['rows']],
     composite_PII_36=[r['composite'][0] for r in nu2['rows']], composite_PII_5=[r['composite'][1] for r in nu2['rows']], composite_PII_0p5=[r['composite'][2] for r in nu2['rows']],
     ratio_PII_36=[r['residual'][0] for r in nu2['rows']], ratio_PII_5=[r['residual'][1] for r in nu2['rows']], ratio_PII_0p5=[r['residual'][2] for r in nu2['rows']]),
     "S1 Data: the cascade composite when the twelve-site downstream ladder is non-binomial (K_i = (n-i+1)/i Khat2 w^(i-1)); ratio = reported/calculated; from scripts/nu2_major2.py")
with open(os.path.join(ROOT, "results", "revision_numbers.json"), "w") as f: json.dump(out, f, indent=1, default=float)
print("written", len(os.listdir(FD)), "csv files and revision_numbers.json")

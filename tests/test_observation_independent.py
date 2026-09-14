"""Independent audit of printed guarantees, source-rate dynamics and decisions.

This complements saved-result checks: high-precision interval calculations,
linear-program sharpness, direct mass-action integration from unmodified pools,
and a refit using a free-ligand root instead of the production quadratic.
"""
from pathlib import Path
import csv
import json
import sys
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares, linprog
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / 'extras/observation_design/code'
OUT = ROOT / 'extras/observation_design/results'
sys.path.insert(0, str(CODE))
from model_source import Cascade


def main():
    report = {}
    # The improved GlnD fit releases the transported UR half-range as well
    # as the fixed-C/middle-state restrictions. Audit both stored fits.
    earlier = json.loads((ROOT / 'extras/finite_mechanism/results/summary.json').read_text())['strict_cross_assay_transfer']
    fitted = json.loads((OUT / 'constrained_glnd.json').read_text())
    q0, q2, middle, logs, logc = earlier['parameters']
    q = np.array([q0, q2 + middle * (q0 - q2), q2])
    scale = np.array([1., (1-q0)/(1-q[1]), fitted['UR_fold']*(1-q0)/(1-q2)])
    old_a, old_b = scale*q, scale*(1-q)
    s, c = np.exp(logs), np.exp(logc)
    old_k1, old_k2 = s*np.sqrt(c)*scale[1], s*scale[2]/(np.sqrt(c)*scale[1])
    comparison = []
    for name, aa, bb, k1, k2, rmse in [
        ('fixed_C', old_a, old_b, old_k1, old_k2, earlier['RMSE']),
        ('two_endpoint_ratios', np.array(fitted['a']), np.array(fitted['b']),
         fitted['K1_mM'], fitted['K2_mM'], fitted['RMSE'])
    ]:
        midpoint = (bb[0]+bb[2])/2
        def capacity(g):
            weights = np.array([1., g/k1, g*g/(k1*k2)])
            return weights@bb/weights.sum()
        half = brentq(lambda g: capacity(g)-midpoint, 1e-10, 100, xtol=1e-14)
        coef = [(bb[2]-midpoint)/(k1*k2), (bb[1]-midpoint)/k1, bb[0]-midpoint]
        positive_root = max(np.roots(coef))
        assert abs(half-positive_root) < 1e-10
        assert abs(aa[0]/bb[0]-fitted['ratio_unliganded']) < 1e-10
        assert abs(bb[2]/bb[0]-fitted['UR_fold']) < 1e-10
        comparison.append(dict(model=name, RMSE=rmse, UR_half_range_mM=float(half),
                               state_capacities=(aa+bb).tolist(),
                               endpoint_ratios=[float(aa[0]/bb[0]), float(bb[2]/bb[0])]))
    assert abs(comparison[0]['UR_half_range_mM']-.080) < 1e-10
    assert abs(comparison[1]['UR_half_range_mM']-.2566242923150519) < 1e-7
    assert not np.allclose(scale, scale[0])
    report['GlnD_constraint_comparison'] = dict(cases=comparison,
        interpretation='The two endpoint ratios are retained; fixed C, the middle-state constraint and the UR half-range are jointly relaxed. RMSE improvement is not attributable to one relaxation.')
    mp.mp.dps = 65
    d = mp.mpf
    eps, e, pt = d('.108'), d('.04'), d('5')
    lower, upper = pt * (1 - eps), pt
    p0 = 1 - d('.67') / (1 - eps)
    g = d('2') / d('15.6')
    ell = (d(175) / d('7.5')) * g * d('.18') / (
        p0 * (1 + (d(1315) / d('7.5')) * g / d('.17')))
    gs_bound = lambda x: e + mp.tanh(3 * mp.log((1 + x / lower) / (1 + x / upper)))
    threshold = mp.findroot(lambda x: gs_bound(x) - d('.05'), (d('.1'), d('.2')))
    assert gs_bound(d('.1417')) < d('.05')
    assert gs_bound(d('.141799')) < d('.05')
    assert gs_bound(d('.142')) > d('.05')  # Previously rounded upward.
    pool_limits = {}
    for n, printed in [(1, '.03536'), (12, '.002995')]:
        exact = 1 - mp.exp(-4 * mp.atanh(d('.009')) / n)
        assert d(printed) < exact
        assert d('.001') + mp.tanh(-n * mp.log(1 - d(printed)) / 4) < d('.01')
        pool_limits[str(n)] = {'exact': str(exact), 'printed': printed}
    main_tex = ROOT / 'manuscript/manuscript_submission.tex'
    if main_tex.exists():
        assert r'\ell\leq0.1417\,\mu' in main_tex.read_text()
        assert '0.002995 for the general twelve-site ladder' in (ROOT / 'manuscript/S1_Appendix.tex').read_text()
    report['high_precision'] = dict(nominal_bound=str(gs_bound(ell)), exact_ell_threshold=str(threshold),
                                  printed_ell_bound=str(gs_bound(d('.1417'))), pool_limits=pool_limits)

    # Positive ladders approaching a two-endpoint distribution attain the
    # tanh(n*d/4) constant. Evaluate from high-precision normalized weights.
    worst = d(0)
    for n in [1, 3, 12]:
        for separation in [d('.0001'), d('.2'), d('2')]:
            coef = [d(1) if i in [0, n] else d('1e-50') for i in range(n + 1)]
            def mean(logu):
                w = [coef[i] * mp.exp(i * logu) for i in range(n + 1)]
                return sum(i * w[i] for i in range(n + 1)) / (n * sum(w))
            error = abs(mean(separation / 2) - mean(-separation / 2) - mp.tanh(n * separation / 4))
            assert error < d('1e-45')
            worst = max(worst, error)
    report['ladder_extremizers'] = {'cases': 9, 'max_error': str(worst)}

    # Charnes-Cooper variables: normalized free fractions y, bound/free z,
    # and t=1/free fraction. This independently minimizes free pi_0.
    lp_error = 0.
    for capacity in [0., .108, .8]:
        for observed in [0., .1, .67, 1., 2.9]:
            ae = np.zeros((2, 9)); ae[0, :4] = 1; ae[1, 4:8] = 1; ae[1, 8] = -1
            au = np.r_[np.arange(4), np.arange(4), -observed][None, :]
            fit = linprog(np.r_[1., np.zeros(8)], A_ub=au, b_ub=[0], A_eq=ae, b_eq=[1, -1],
                          bounds=[(0, None)] * 8 + [(1, 1 / (1 - capacity))], method='highs')
            assert fit.success
            error = abs(fit.fun - max(0, 1 - observed / (1 - capacity)))
            assert error < 1e-9
            lp_error = max(lp_error, error)
    report['count_gate_linear_programs'] = {'cases': 15, 'max_error': lp_error}

    # No stationary species are used to initialize these trajectories.
    cases = []
    for lam in [1, 4]:
        model = Cascade(G=2000, Pt=5, St=.5, Dt=.5, Et=.02, beta=175/7.5,
                        beta_prime=1315/7.5, gain=.001, lam=lam)
        n = len(model.names); totals = np.array([model.Pt, model.St, model.Dt, model.Et])
        y = np.zeros(n); y[[model.P[0], model.S[0], model.D[0], model.E[0]]] = totals
        reference = model.steady(); ii = np.arange(len(model.ks)); mask = model.r2 < n
        def jac(time, state):
            yy = np.append(state, 1.); v = np.zeros((len(model.ks), n))
            np.add.at(v, (ii, model.r1), model.ks * yy[model.r2])
            np.add.at(v, (ii[mask], model.r2[mask]), model.ks[mask] * yy[model.r1[mask]])
            return model.N @ v
        t0 = 0.; minimum = 0.; conservation = 0.
        for t1 in [1e3, 1e4, 1e5, 1e6, 1e7]:
            sol = solve_ivp(model.rhs, (t0, t1), y, method='BDF', jac=jac, rtol=2e-9, atol=2e-12)
            assert sol.success
            conservation = max(conservation, float(np.max(abs(model.cons @ sol.y - totals[:, None]))))
            minimum = min(minimum, float(sol.y.min()))
            y = sol.y[:, -1]; t0 = t1
            error = float(np.max(abs(y - reference)))
            residual = float(np.max(abs(model.rhs(t1, y))))
            if error < 1e-7 and residual < 1e-10:
                break
        assert error < 1e-7 and residual < 1e-10 and minimum > -1e-10 and conservation < 1e-8
        cases.append(dict(lam=lam, time=t1, max_species_error=error, max_rhs=residual,
                          conservation_error=conservation, min_concentration=minimum,
                          GS_total=model.obs(y)['GStotal']))
    report['source_rate_ODE'] = dict(parameters=dict(G=2000, Pt=5, St=.5, Dt=.5, Et=.02,
                                                   beta=175/7.5, beta_prime=1315/7.5, gain=.001), cases=cases)

    # Independently refit the archived raw example with a free-ligand root,
    # a different starting point and unconstrained Levenberg-Marquardt.
    example = json.loads((OUT / 'binding_observation_example.json').read_text())
    doses = np.array(example['P_nominal_uM']); obs = np.array(example['observations']); enz = example['E_measured_uM']
    def residual(z):
        h = np.exp(z[:3]); pred = np.empty((3, len(doses)))
        for j in range(3):
            for i, dose in enumerate(doses):
                free = 0 if dose == 0 else brentq(lambda p: p + enz*p/(h[j]+p) - dose, 0, dose, xtol=1e-14)
                pred[j, i] = z[3] + z[4] * free / (h[j] + free)
        return (pred[:, :, None] - obs).ravel()
    fit = least_squares(residual, np.r_[np.log([.1, .08, .05]), .02, .95], method='lm',
                        ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=1000)
    assert fit.success
    h = np.exp(fit.x[:3]); error = float(max(abs(h - example['H_estimated_uM'])))
    assert error < 1e-6
    report['independent_raw_refit'] = dict(max_halfpoint_error_uM=error,
                                         rmse=float(np.sqrt(np.mean(fit.fun**2))))

    null = sum(mp.log((1 + x) / (1 + x / d('.17'))) for x in [d(1), d(2)]) / 2
    counted = 0
    for stem, expected in [('binding_observation', 7500), ('binding_local_effects', 1200)]:
        rows = list(csv.DictReader((OUT / (stem + '_trials.csv')).open()))
        assert len(rows) == expected
        for row in rows:
            statistic = (float(row['mean_score']) - float(null)) / float(row['score_se'])
            decision = bool(int(row['valid_fit'])) and t.sf(statistic, df=7) < .05
            assert decision == bool(int(row['rejected']))
        counted += len(rows)
    report['independently_recomputed_trial_decisions'] = counted
    report['passed'] = True
    report['scope'] = 'Conditional computational validation; two integrations do not prove global convergence or native kinetics.'
    (OUT / 'independent_verification.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

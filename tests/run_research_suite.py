"""One entry point for all existing research layers.
--regenerate recreates numerical results before assertions; the default mode reads
baseline saved outputs but freshly solves extension identities and reaction models.
"""
from pathlib import Path
import subprocess,sys,json,time,argparse
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--regenerate',action='store_true');a=p.parse_args()
commands=[]
if a.regenerate:
 (R/'results').mkdir(exist_ok=True)
 commands=[['run_all.py'],['extras/matched_models/run.py'],['extras/finite_mechanism/run.py'],['extras/regulatory_states/code/analyze.py'],['extras/robust_design/code/analyze.py']]+[[f'extras/closed_cascade/code/{n}.py'] for n in ['analyze','topology','total_input','verify_algebra']]
if a.regenerate:
 commands += [[f'extras/paired_paralogs/code/{n}.py'] for n in ['extract_data','state_bounds']]
commands += [['tests/test_reproduce.py'],['extras/matched_models/code/verify_final.py'],['tests/test_extensions.py']]
commands += [[f'extras/paired_paralogs/code/{n}.py'] for n in ['paired_reporter','reporter_intervals','robust_calibration']]
if a.regenerate:
 commands += [[f'extras/observation_design/code/{n}.py'] for n in ['observation_certificate','verify_observation_certificate','constrained_glnd','glnd_constraint_profile','binding_observation','binding_local_effects']]
commands += [['extras/observation_design/code/verify_revision.py']]
commands += [['tests/test_observation_independent.py']]
(R/'results').mkdir(exist_ok=True)
report=[]
for args in commands:
 start=time.time();print('RUN',*args,flush=True)
 status=subprocess.run([sys.executable,*args],cwd=R).returncode
 report.append(dict(command=args,exit_code=status,seconds=time.time()-start))
 (R/'results/research_suite.json').write_text(json.dumps({'regenerated':a.regenerate,'steps':report,'completed':len(report)==len(commands),'passed':len(report)==len(commands) and all(s['exit_code']==0 for s in report)},indent=2))
 if status:sys.exit(status)
print('All research layers completed',flush=True)

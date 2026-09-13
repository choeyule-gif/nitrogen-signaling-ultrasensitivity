"""Reproduce the new analyses and figures in dependency order."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for script in ['analyze.py','extend_and_verify.py','joint_ladder.py','asymmetric_gs.py','verify_final.py','figures.py']:
    print('\nRunning',script,flush=True)
    subprocess.run([sys.executable,str(root/'code'/script)],check=True)

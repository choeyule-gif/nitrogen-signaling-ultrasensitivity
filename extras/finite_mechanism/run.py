"""Reproduce finite-mechanism results and Figure 7."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for name in ['analyze.py','verify_network.py']:
    subprocess.run([sys.executable,str(root/'code'/name)],check=True)

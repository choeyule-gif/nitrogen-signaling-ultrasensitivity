"""Repository-relative paths.  No absolute paths appear anywhere else in the code."""
import os
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
S1DATA = os.path.join(ROOT, "S1_Data")
for _d in (RESULTS, FIGURES, S1DATA):
    os.makedirs(_d, exist_ok=True)

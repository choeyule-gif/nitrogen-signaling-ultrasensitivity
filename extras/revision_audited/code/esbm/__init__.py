"""esbm -- effector-switched bifunctional modification.

Shared model code for the manuscript "Two separable sources of ultrasensitivity in
bifunctional modification cycles, and how a bacterial nitrogen cascade uses them".
Every script in scripts/ imports from here; nothing in here is specific to one figure.
"""
from .paths import ROOT, DATA, RESULTS, FIGURES, S1DATA
from .report import save_json, load_json, save_csv
__all__ = ["ROOT", "DATA", "RESULTS", "FIGURES", "S1DATA", "save_json", "load_json", "save_csv"]

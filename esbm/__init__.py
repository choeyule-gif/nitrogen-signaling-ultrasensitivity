"""esbm -- effector-switched bifunctional modification.

Shared model code for "Mechanistic equivalence and measurement design in models
of bacterial nitrogen signaling".
Every script in scripts/ imports from here; nothing in here is specific to one figure.
"""
from .paths import ROOT, DATA, RESULTS, FIGURES, S1DATA
from .report import save_json, load_json, save_csv
__all__ = ["ROOT", "DATA", "RESULTS", "FIGURES", "S1DATA", "save_json", "load_json", "save_csv"]

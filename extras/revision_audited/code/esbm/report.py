"""Write the numbers a script reports to results/<name>.json, so that tests/ can
compare them with the values quoted in the manuscript."""
import json, os, numpy as np
from .paths import RESULTS

def _clean(x):
    if isinstance(x, dict): return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [_clean(v) for v in x]
    if isinstance(x, np.ndarray): return _clean(x.tolist())
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    return x

def save_json(name, obj):
    path = os.path.join(RESULTS, name + ".json")
    with open(path, "w") as f:
        json.dump(_clean(obj), f, indent=1)
    return path

def load_json(name):
    with open(os.path.join(RESULTS, name + ".json")) as f:
        return json.load(f)

def save_csv(name, columns, header_comment=None):
    """Write results/panels/<name>.csv from a dict of equal-length columns (the numbers behind one
    panel of a Python figure).  scripts/export_figure_data.py draws on these for the figdata/ files
    the MATLAB figures read, and scripts/make_s1_data.py assembles the S1 Data package."""
    import csv
    from .paths import RESULTS
    S1DATA = os.path.join(RESULTS, "panels"); os.makedirs(S1DATA, exist_ok=True)
    keys = list(columns); rows = zip(*[np.asarray(columns[k]).ravel() for k in keys])
    path = os.path.join(S1DATA, name + ".csv")
    with open(path, "w", newline="") as f:
        if header_comment: f.write("# " + header_comment + "\n")
        w = csv.writer(f); w.writerow(keys); w.writerows(rows)
    return path

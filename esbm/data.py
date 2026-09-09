"""Published data used by the reanalysis, read from data/.

digitised_titration()  -- the twelve points of Fig. 2B of Jiang, Peliska & Ninfa (1998),
                          Biochemistry 37:12782 (PII uridylylation state against glutamine)
constants()            -- published Hill coefficients, half-points and affinities
                          (Ventura et al. 2010; Jiang & Ninfa 2011; Jiang et al. 1998)
"""
import csv, json, os
import numpy as np
from .paths import DATA

def digitised_titration():
    g, U = [], []
    with open(os.path.join(DATA, "jiang1998_fig2B_digitised.csv")) as f:
        for row in csv.DictReader(f):
            g.append(float(row["glutamine_mM"])); U.append(float(row["uridylylation_state"]))
    return np.array(g), np.array(U)

def constants():
    with open(os.path.join(DATA, "published_constants.json")) as f:
        return json.load(f)

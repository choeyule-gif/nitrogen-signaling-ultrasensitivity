"""Create compact main-text panels from checked supplementary calculations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIGDATA = ROOT / "figdata"
AUDIT = ROOT / "data" / "audit"
OUT = ROOT / "figures" / "promoted"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
BLUE = "#4d88cf"
ROSE = "#d87978"


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, comments="#", skip_header=1)


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


# Fig. 1C: where a total-target measurement departs from the free-target invariant.
t = read_csv(FIGDATA / "FigS3_plateau_surface.csv")
xvals = np.unique(t["TE_over_TS"])
yvals = np.unique(t["TS_over_Km"])
z = np.full((len(yvals), len(xvals)), np.nan)
xi = {v: i for i, v in enumerate(xvals)}
yi = {v: i for i, v in enumerate(yvals)}
for row in t:
    z[yi[row["TS_over_Km"]], xi[row["TE_over_TS"]]] = row["nH_local"]
fig, ax = plt.subplots(figsize=(3.3, 2.55), layout="constrained")
mesh = ax.pcolormesh(xvals, yvals, z, shading="auto", cmap="viridis", vmin=0.65, vmax=1.0)
ax.set(xscale="log", yscale="log", xlabel=r"$T_E/T_S$", ylabel=r"$T_S/K_m$")
ax.axvline(0.01, color="white", ls="--", lw=1.2)
cb = fig.colorbar(mesh, ax=ax, pad=0.02)
cb.set_label(r"Total-readout $n_{H,\mathrm{local}}$")
save(fig, "Fig1C_total_heatmap")


# Fig. 1D: representative slices make the concentration trend readable at panel size.
fig, ax = plt.subplots(figsize=(3.3, 2.55), layout="constrained")
for target, color in zip([1e-4, 1e-2, 1e-1], [BLUE, "#33a99e", ROSE]):
    xv = xvals[np.argmin(np.abs(np.log10(xvals) - np.log10(target)))]
    rows = t[t["TE_over_TS"] == xv]
    order = np.argsort(rows["TS_over_Km"])
    ax.semilogx(rows["TS_over_Km"][order], rows["nH_local"][order], color=color,
                lw=1.8, label=rf"$T_E/T_S={xv:g}$")
ax.axhline(1, color="0.55", ls="--", lw=1)
ax.set(xlabel=r"$T_S/K_m$", ylabel=r"Total-readout $n_{H,\mathrm{local}}$", ylim=(0.65, 1.015))
ax.legend(fontsize=7, loc="lower left")
save(fig, "Fig1D_total_slices")


# Fig. 2D: a two-dimensional projection emphasizes estimator-dependent non-identifiability.
s = read_csv(FIGDATA / "FigS2_nu_sup_surface.csv")
avals = np.unique(s["log10_a"])
bvals = np.unique(s["log10_b"])
surface = np.full((len(bvals), len(avals)), np.nan)
ai = {v: i for i, v in enumerate(avals)}
bi = {v: i for i, v in enumerate(bvals)}
for row in s:
    surface[bi[row["log10_b"]], ai[row["log10_a"]]] = row["nu_sup"]
fig, ax = plt.subplots(figsize=(3.3, 2.55), layout="constrained")
levels = np.linspace(1, 3, 17)
mesh = ax.contourf(avals, bvals, surface, levels=levels, cmap="viridis", extend="max")
for suffix, color, label in [("fit", BLUE, "Fitted exponent"), ("range", ROSE, "Response range")]:
    locus = read_csv(FIGDATA / f"FigS2_locus_{suffix}.csv")
    ax.plot(locus["log10_a"], locus["log10_b"], color=color, lw=2, label=label)
ax.scatter([0], [0], facecolor="white", edgecolor="black", s=32, zorder=4)
ax.set(xlabel=r"$\log_{10}a$", ylabel=r"$\log_{10}b$", xlim=(-3, 3), ylim=(-3, 3))
cb = fig.colorbar(mesh, ax=ax, pad=0.02)
cb.set_label(r"$\sup_u\nu$")
ax.legend(fontsize=7, loc="lower right")
save(fig, "Fig2D_ladder_constraints")


# Fig. 4D: the checked direct-route calculation, normalized by the fitted midpoint S.
d = np.genfromtxt(AUDIT / "direct_route_checked.csv", delimiter=",", names=True)
S = json.loads((ROOT / "data" / "analysis_summary.json").read_text())["parameters"]["S_mM"]
fig, ax = plt.subplots(figsize=(3.3, 2.55), layout="constrained")
ax.semilogx(d["KG"] / S, d["gain"], color=BLUE, lw=2)
ax.axhline(1, color="0.55", ls="--", lw=1)
ax.set(xlabel=r"$K_G/S$", ylabel="Direct-route gain", ylim=(0.93, 1.19))
ax.text(0.96, 0.94, r"Illustrative $R=3.7$", transform=ax.transAxes,
        ha="right", va="top", fontsize=8)
save(fig, "Fig4D_direct_route")

print("Four promoted main-text panels prepared.")

"""Matplotlib styling for the figures: MATLAB-like colours, inward ticks, no minor ticks,
no prose inside the axes (numbers, symbols and legends only)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator, MaxNLocator, LogLocator

MC = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F']
GRAY = '#BFBFBF'

COMPACT = {'font.size': 9, 'axes.linewidth': 0.8, 'lines.linewidth': 1.6,
           'xtick.direction': 'in', 'ytick.direction': 'in',
           'xtick.top': True, 'ytick.right': True,
           'xtick.major.size': 4, 'ytick.major.size': 4,
           'axes.labelsize': 10, 'legend.fontsize': 8,
           'legend.frameon': True, 'legend.edgecolor': '0', 'legend.framealpha': 1,
           'pdf.fonttype': 42}

def use_compact(**over):
    """Style of Figs 1, 2, 4, 7, 8: explicit axes placement, legends boxed."""
    d = dict(COMPACT); d.update(over); plt.rcParams.update(d)

def use():
    """Style of Figs 5 and 6 (subplot grid with tidy())."""
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans', 'font.size': 9, 'axes.labelsize': 10,
        'axes.linewidth': 0.9, 'axes.edgecolor': 'black', 'axes.facecolor': 'white',
        'axes.grid': False, 'axes.axisbelow': True,
        'axes.spines.top': True, 'axes.spines.right': True,
        'axes.prop_cycle': mpl.cycler(color=MC),
        'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
        'xtick.major.size': 4.0, 'ytick.major.size': 4.0, 'xtick.minor.size': 0.0, 'ytick.minor.size': 0.0,
        'xtick.major.width': 0.9, 'ytick.major.width': 0.9, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
        'lines.linewidth': 1.6, 'lines.markersize': 5,
        'legend.frameon': True, 'legend.edgecolor': 'black', 'legend.fancybox': False,
        'legend.framealpha': 1.0, 'legend.fontsize': 9, 'legend.borderpad': 0.35,
        'legend.handlelength': 1.7, 'legend.labelspacing': 0.28,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
        'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02, 'pdf.fonttype': 42, 'ps.fonttype': 42})

def tidy(ax, logx=False, logy=False, nx=4, ny=4):
    ax.xaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=nx) if logx else MaxNLocator(nbins=nx))
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=ny) if logy else MaxNLocator(nbins=ny))
    for s in ax.spines.values(): s.set_zorder(10)
    return ax

def panel(ax, letter, dx=-0.155, dy=1.045):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11, fontweight='bold', va='top', ha='left')

def no_minor(fig):
    for a in fig.axes:
        a.xaxis.set_minor_locator(NullLocator()); a.yaxis.set_minor_locator(NullLocator())

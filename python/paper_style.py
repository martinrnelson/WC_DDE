"""Colours and sizes shared by the Python figures (same as tools/paper_style.m).

Solid line = stable, dashed = unstable; colour = type of object
(Okabe-Ito palette, colour-blind safe).
"""
import matplotlib.pyplot as plt

PO    = '#0072B2'   # periodic orbits                blue
PO2   = '#E69F00'   # period-doubled orbits          orange
HOPF  = '#D55E00'   # Hopf points and curves         vermillion
TORUS = '#CC79A7'   # torus points and curves        reddish purple
PD    = '#009E73'   # period-doubling points/curves  bluish green
FOLD  = '#56B4E9'   # folds of periodic orbits       sky blue
GREY  = '0.55'
CM = 1 / 2.54       # inches per cm (Matplotlib sizes are in inches)


def use():
    """Matplotlib settings close to the MATLAB figures: serif text, LaTeX-style maths."""
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm',
                         'font.size': 9, 'axes.linewidth': 0.6,
                         'xtick.direction': 'out', 'ytick.direction': 'out',
                         'legend.frameon': False})

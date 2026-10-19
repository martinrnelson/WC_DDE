"""Hopf curves of the rich parameter set in the (b, tau) plane (Section 3.1).

The analytic part of Figure 5: no continuation is needed. For each
cross-coupling b = c, the symmetric equilibrium solves u* = f((a+b)u* + P);
the in-phase and anti-phase modes have mu = k(a+b) and mu = k(a-b), with
k = beta u*(1 - u*), and each gives a family of Hopf delays. Where an
in-phase curve crosses an anti-phase curve there is a double Hopf point.
The torus and period-doubling curves of Figure 5 need DDE-BIFTOOL (MATLAB
script compute/c4_two_parameter.m).

Produces: pdf/hopf_curves.pdf
Usage:    python hopf_curves.py
"""
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import wcmodel as wc
import paper_style as ps

here = os.path.dirname(os.path.abspath(__file__))
p = wc.parameters('rich')
f = lambda z: 1 / (1 + math.exp(-p['beta'] * z))

b_values = np.linspace(0.05, 5.5, 500)
tau_in = np.full((2, b_values.size), np.nan)     # first 2 in-phase Hopf delays
tau_anti = np.full((4, b_values.size), np.nan)   # first 4 anti-phase Hopf delays
for i, b in enumerate(b_values):
    ustar = brentq(lambda u: u - f((p['a'] + b) * u + p['P']), 0, 1)
    k = p['beta'] * ustar * (1 - ustar)
    tau_in[:, i] = wc.hopf_delays_from_mu(k * (p['a'] + b), 2)[0]
    tau_anti[:, i] = wc.hopf_delays_from_mu(k * (p['a'] - b), 4)[0]

# double Hopf points: sign changes of (in-phase - anti-phase), interpolated
double_hopf = []
for ti in tau_in:
    for ta in tau_anti:
        d = ti - ta
        for q in np.where(d[:-1] * d[1:] < 0)[0]:
            s = d[q] / (d[q] - d[q + 1])
            double_hopf.append((b_values[q] + s * (b_values[q + 1] - b_values[q]),
                                ti[q] + s * (ti[q + 1] - ti[q])))
double_hopf = np.array(sorted(double_hopf, key=lambda z: z[1]))
for b, tau in double_hopf[double_hopf[:, 1] < 0.35]:
    print(f'double Hopf point at b = {b:.3f}, tau = {tau:.3f}')

ps.use()
fig, ax = plt.subplots(figsize=(10 * ps.CM, 8 * ps.CM))
ax.fill_between(b_values, 0, tau_anti[0], color='0.92', linewidth=0,
                label='stable equilibrium')
ax.plot(b_values, tau_anti.T, '-', color=ps.HOPF, linewidth=1.2)
ax.plot(b_values, tau_in.T, '-.', color=ps.HOPF, linewidth=1.2)
ax.plot(double_hopf[:, 0], double_hopf[:, 1], 'ko', markersize=3, label='double Hopf')
ax.plot([], [], '-', color=ps.HOPF, label='Hopf (anti-phase)')
ax.plot([], [], '-.', color=ps.HOPF, label='Hopf (in-phase)')
ax.axvline(p['b'], color='0.5', linestyle=':', linewidth=0.8)
ax.set_xlim(b_values[0], b_values[-1])
ax.set_ylim(0, 0.35)
ax.set_xlabel('$b = c$')
ax.set_ylabel(r'$\tau$')
fig.legend(loc='lower center', ncol=2)
fig.tight_layout(rect=(0, 0.14, 1, 1))
out = os.path.join(here, 'pdf', 'hopf_curves.pdf')
fig.savefig(out)
print('Saved', out)

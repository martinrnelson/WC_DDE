"""Figure 6 (Python version): route to chaos through torus breakdown.

(a) orbit diagram: local maxima of u after transients, by simulation
(b) largest Lyapunov exponent
Dashed lines mark the Hopf point, the torus point and the onset of chaos.
Data: data/route_to_chaos.npz, made by route_to_chaos.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import wcmodel as wc
import paper_style as ps

here = os.path.dirname(os.path.abspath(__file__))
D = np.load(os.path.join(here, 'data', 'route_to_chaos.npz'))
ps.use()

tau_H = wc.hopf_delays(wc.parameters('rich'))['tau_c']   # exact (Section 3.1)
tau_T = 0.0639             # torus point from DDE-BIFTOOL (MATLAB script c2)
marks = [(tau_H, ps.HOPF, r'$\tau_H$'), (tau_T, ps.TORUS, r'$\tau_T$'),
         (float(D['tau_chaos']), '0.35', r'$\tau_C$')]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(14 * ps.CM, 10 * ps.CM))
ax1.plot(D['points_tau'], D['points_umax'], 'k.', markersize=1)
ax1.set_ylabel('local maxima of $u$')
ax1.set_title('(a)')
ax2.axhline(0, color='0.5', linestyle=':', linewidth=0.8)
ax2.plot(D['lyap_tau'], D['lyap_lambda'], '.-', color=ps.PO, linewidth=1, markersize=3)
ax2.set_xlabel(r'$\tau$')
ax2.set_ylabel(r'$\lambda_{\max}$')
ax2.set_title('(b)')
for ax in (ax1, ax2):
    for x, colour, _ in marks:
        ax.axvline(x, color=colour, linestyle='--', linewidth=0.8)
    ax.set_xlim(D['diagram_tau'][0], D['diagram_tau'][-1])
for x, colour, label in marks:
    ax1.annotate(label, (x, 1), xycoords=('data', 'axes fraction'),
                 xytext=(2, -2), textcoords='offset points', va='top', color=colour)
fig.tight_layout()
out = os.path.join(here, 'pdf', 'fig6_route_to_chaos.pdf')
fig.savefig(out)
print('Saved', out)

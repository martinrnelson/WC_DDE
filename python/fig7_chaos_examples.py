"""Figure 7 (Python version): solutions along the route to chaos.

One column per solution: periodic, quasi-periodic, mode-locked, chaotic.
Top:    u(t) (blue) and v(t) (grey) over a short time window
Middle: the solution in the (u, v) plane
Bottom: return map of successive maxima of u: a point (periodic), a closed
        curve (quasi-periodic), a finite set (mode-locked), a folded set (chaotic)
Data: data/route_to_chaos.npz, made by route_to_chaos.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import paper_style as ps

here = os.path.dirname(os.path.abspath(__file__))
D = np.load(os.path.join(here, 'data', 'route_to_chaos.npz'))
ps.use()
n = len(D['example_tau'])

fig, axes = plt.subplots(3, n, figsize=(17 * ps.CM, 12 * ps.CM))
for k in range(n):
    ax = axes[0, k]
    ax.plot(D[f't{k}'], D[f'u{k}'], color=ps.PO, linewidth=1, label='$u$')
    ax.plot(D[f't{k}'], D[f'v{k}'], color='0.6', linewidth=1, label='$v$')
    ax.set_xlim(0, D[f't{k}'].max())
    ax.set_xlabel('$t$')
    ax.set_title(f"({'abcd'[k]}) {D['example_label'][k]}\n"
                 rf"$\tau = {D['example_tau'][k]:.4g}$")

    ax = axes[1, k]
    ax.plot(D[f'pu{k}'], D[f'pv{k}'], color=ps.PO, linewidth=0.4)
    ax.set_xlabel('$u$')

    ax = axes[2, k]
    m = D[f'max{k}']
    ax.plot(m[:-1], m[1:], 'k.', markersize=2)
    lo, hi = m.min(), m.max()
    pad = max(0.05 * (hi - lo), 1e-4)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$u_{\max}(n)$')
    for a in axes[:, k]:
        a.tick_params(labelsize=7)
        a.locator_params(nbins=3)

axes[0, 0].set_ylabel('activity')
axes[1, 0].set_ylabel('$v$')
axes[2, 0].set_ylabel(r'$u_{\max}(n+1)$')
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=2)
fig.tight_layout(rect=(0, 0.05, 1, 1))
out = os.path.join(here, 'pdf', 'fig7_chaos_examples.pdf')
fig.savefig(out)
print('Saved', out)

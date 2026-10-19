"""Python version of compute/c5_route_to_chaos.m.

Rich parameter set: orbit diagram, largest Lyapunov exponent and four
example solutions, all by direct simulation with the Euler map.

Paper:    Section 6; Figures 6 and 7
Produces: data/route_to_chaos.npz, used by fig6_route_to_chaos.py and
          fig7_chaos_examples.py
Run time: about 15 minutes (most of it in step 2)

Usage:    python route_to_chaos.py
"""
import os
import time
import numpy as np
import wcmodel as wc

here = os.path.dirname(os.path.abspath(__file__))
p = wc.parameters('rich')
history = (0.3, 0.1)      # u and v must differ, or the solution stays on u = v
                          # and never shows the anti-phase dynamics

# 1. Orbit diagram: local maxima of u after a transient, for each tau.
#    A periodic solution gives a few distinct maxima, a quasi-periodic one
#    fills an interval, and a chaotic one a scattered set.
diagram_tau = np.arange(0.040, 0.130 + 1e-9, 5e-4)
T_transient, T_record = 200, 50
points_tau, points_umax = [], []
start = time.time()
for tau in diagram_tau:
    t, u, v = wc.simulate(p, tau, T_transient + T_record, history)
    peaks, _ = wc.local_maxima(u[t > T_transient])
    if len(peaks) == 0:              # a steady state has no maxima:
        peaks = u[-1:]               # record its value instead
    points_tau.extend([tau] * len(peaks))
    points_umax.extend(peaks)
print(f'Orbit diagram: {len(diagram_tau)} delays in {time.time() - start:.0f} s')

# 2. Largest Lyapunov exponent for each tau
lyap_tau = np.arange(0.040, 0.130 + 1e-9, 5e-4)
start = time.time()
lyap_lambda = np.array([wc.lyapunov_exponent(p, tau, 200, 200, history)
                        for tau in lyap_tau])
print(f'Lyapunov exponents: {len(lyap_tau)} delays in {time.time() - start:.0f} s')
tau_chaos = lyap_tau[np.argmax(lyap_lambda > 0.1)]   # first clearly positive
print(f'lambda first exceeds 0.1 at tau = {tau_chaos:.4f}')

# 3. Example solutions for Figure 7
example_tau = np.array([0.060, 0.068, 0.0735, 0.100])
example_label = np.array(['periodic', 'quasi-periodic', 'mode-locked', 'chaotic'])
T = 400
examples = {}
for k, tau in enumerate(example_tau):
    t, u, v = wc.simulate(p, tau, T, history)
    w = t > T - 2                          # short window for the time series
    examples[f't{k}'], examples[f'u{k}'], examples[f'v{k}'] = t[w] - (T - 2), u[w], v[w]
    w = t > T - 40                         # longer window for the (u,v) plane
    examples[f'pu{k}'], examples[f'pv{k}'] = u[w], v[w]
    examples[f'max{k}'], _ = wc.local_maxima(u[t > 100])   # for the return map

# 4. Save
np.savez(os.path.join(here, 'data', 'route_to_chaos.npz'),
         diagram_tau=diagram_tau, points_tau=np.array(points_tau),
         points_umax=np.array(points_umax), lyap_tau=lyap_tau,
         lyap_lambda=lyap_lambda, tau_chaos=tau_chaos,
         example_tau=example_tau, example_label=example_label, **examples)
print('Saved data/route_to_chaos.npz')

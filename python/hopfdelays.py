"""Listing 1 (Python version): Hopf delays of the symmetric Wilson-Cowan equilibrium.

For a = d, b = c, P = Q the linearisation splits into an in-phase and an
anti-phase mode, each satisfying lambda + 1 = mu exp(-lambda tau).
"""
import math
import numpy as np
from scipy.optimize import brentq

a, b, P, beta = -6.0, 2.5, 0.2, 60.0
f = lambda z: 1 / (1 + math.exp(-beta * z))

ustar = brentq(lambda u: u - f((a + b) * u + P), 0, 1)   # equilibrium u* = v*
k = beta * ustar * (1 - ustar)                           # slope f'(z*)
mu = np.array([k * (a + b), k * (a - b)])                # in-phase, anti-phase
w = np.sqrt(mu ** 2 - 1)                                 # Hopf frequencies (|mu| > 1)
j = np.arange(4)
tau = (np.angle(mu)[:, None] - np.arctan(w)[:, None] + 2 * np.pi * j) / w[:, None]

print(f'u* = {ustar:.5f},  k = {k:.4f}')
print(f'in-phase:   mu = {mu[0]:7.3f}, omega = {w[0]:6.3f}, tau_j =', np.round(tau[0], 4))
print(f'anti-phase: mu = {mu[1]:7.3f}, omega = {w[1]:6.3f}, tau_j =', np.round(tau[1], 4))

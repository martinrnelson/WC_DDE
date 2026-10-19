"""Listing 6 (Python version): largest Lyapunov exponent of the Wilson-Cowan DDE.

Two nearby solutions are advanced with the Euler map. The state of a DDE is
the whole history on [t - tau, t], so the separation is measured over all
m+1 stored values and renormalised once per delay interval.
"""
import numpy as np

a, b, P, beta, tau = -6, 2.5, 0.2, 60, 0.1      # also try tau = 0.068
W = np.array([[a, b], [b, a]])
f = lambda z: 1 / (1 + np.exp(-beta * z))
F = lambda x, xd: -x + f(W @ xd + P)            # rhs for x(t), x(t - tau)
m = 200; h = tau / m

X = np.tile([[0.3], [0.1]], (1, m + 1))   # circular buffer of the last m+1 states
i = m                                     # column holding the current state
for n in range(round(100 / h)):           # discard a transient of length 100
    j = (i + 1) % (m + 1)                 # oldest column = state at t - tau
    X[:, j] = X[:, i] + h * F(X[:, i], X[:, j])
    i = j

d0 = 1e-8; Y = X + d0 / np.sqrt(X.size)   # perturbed copy, distance d0
S = 0.0; N = round(200 / h)
for n in range(1, N + 1):
    j = (i + 1) % (m + 1)
    X[:, j] = X[:, i] + h * F(X[:, i], X[:, j])
    Y[:, j] = Y[:, i] + h * F(Y[:, i], Y[:, j])
    i = j
    if n % m == 0:                        # once per delay interval
        D = Y - X; s = np.linalg.norm(D)
        S += np.log(s / d0); Y = X + (d0 / s) * D
print('lambda =', S / (N * h))

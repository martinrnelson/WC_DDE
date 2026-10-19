"""Listing 2 (Python version): the Euler map (2.9) for the Wilson-Cowan model.

The step is h = tau/m, so the delayed state x(t - tau) is always exactly m
steps back and no interpolation is needed. Constant history on [-tau, 0].
"""
import numpy as np
import matplotlib.pyplot as plt

a, b, c, d, P, Q, beta = 10, -10, 10, 2, -2, -4, 1
W = np.array([[a, b], [c, d]])
I = np.array([P, Q])
f = lambda z: 1 / (1 + np.exp(-beta * z))
F = lambda x, xd: -x + f(W @ xd + I)        # x = [u, v] at t and t - tau

tau, T = 0.5, 400                           # also try tau = 0.05 and 0.12
m = 100; h = tau / m; N = round(T / h)      # m steps per delay interval
x = np.tile([[0.19644], [0.13377]], (1, N + m + 1))   # columns 0..m = history
for n in range(m, N + m):
    x[:, n + 1] = x[:, n] + h * F(x[:, n], x[:, n - m])
t = (np.arange(N + m + 1) - m) * h

plt.plot(t, x[0], 'b-'); plt.xlabel('t'); plt.ylabel('u(t)'); plt.show()
late = x[0, t > T - 50]
print('amp =', late.max() - late.min())    # oscillation amplitude of u

"""The single-delay Wilson-Cowan model: Python versions of the MATLAB files
in ../model and ../tools that do not need DDE-BIFTOOL.

    u'(t) = -u(t) + f(a u(t-tau) + b v(t-tau) + P)
    v'(t) = -v(t) + f(c u(t-tau) + d v(t-tau) + Q),   f(z) = 1/(1 + exp(-beta z))

Functions
    parameters(name)          the 'gentle' or 'rich' parameter set (Table 1)
    equilibrium(p)            (u*, v*), by Newton's method
    hopf_delays_from_mu(mu)   Hopf delays of lambda + 1 = mu exp(-lambda tau)
    hopf_delays(p)            Hopf delays of the model's equilibrium
    simulate(p, tau, T, history)    the Euler map (2.9)
    local_maxima(u)           peaks of a sampled signal
    lyapunov_exponent(p, tau, ...)  largest Lyapunov exponent

The simulation loops use plain Python numbers (not NumPy arrays) because,
for a system of two equations, that is both simpler and faster.
"""
import math
import numpy as np


def parameters(name):
    """The two parameter sets of the paper (Coombes and Laing, 2009).

    'gentle': beta = 1, a single delay-induced Hopf bifurcation (Figure 1)
    'rich':   beta = 60, Hopf, torus, period doubling and chaos (Figures 2-7)
    Returns a dict with keys a, b, c, d, P, Q, beta and x0 (a guess for the
    equilibrium).
    """
    if name == 'gentle':
        return dict(a=10.0, b=-10.0, c=10.0, d=2.0, P=-2.0, Q=-4.0, beta=1.0,
                    x0=np.array([0.19, 0.13]))
    if name == 'rich':
        return dict(a=-6.0, b=2.5, c=2.5, d=-6.0, P=0.2, Q=0.2, beta=60.0,
                    x0=np.array([0.07, 0.07]))
    raise ValueError("use 'gentle' or 'rich'")


def equilibrium(p):
    """The equilibrium [u*, v*].

    At an equilibrium u(t - tau) = u(t), so the delay drops out and we solve
    G(x) = x - f(W x + I) = 0 by Newton's method, using f' = beta f (1 - f).
    """
    W = np.array([[p['a'], p['b']], [p['c'], p['d']]])
    I = np.array([p['P'], p['Q']])
    x = p['x0'].copy()
    for _ in range(50):
        F = 1 / (1 + np.exp(-p['beta'] * (W @ x + I)))
        G = x - F
        DG = np.eye(2) - np.diag(p['beta'] * F * (1 - F)) @ W
        dx = -np.linalg.solve(DG, G)
        x = x + dx
        if np.linalg.norm(dx) < 1e-13:
            break
    return x


def hopf_delays_from_mu(mu, nj=1):
    """Hopf delays tau_0 < tau_1 < ... of lambda + 1 = mu exp(-lambda tau).

    Section 3.1: omega = sqrt(|mu|^2 - 1) (modulus) and
    omega tau = arg(mu) - arctan(omega) + 2 pi j (argument).
    Returns (tau, omega); NaNs if |mu| <= 1, when there are no Hopf points.
    """
    if abs(mu) <= 1:
        return np.full(nj, np.nan), np.nan
    omega = math.sqrt(abs(mu) ** 2 - 1)
    j = np.arange(nj)
    tau = (np.mod(np.angle(mu) - math.atan(omega), 2 * np.pi) + 2 * np.pi * j) / omega
    return tau, omega


def hopf_delays(p, nj=1):
    """Hopf delays of the equilibrium, for every eigenvalue mu of M.

    The linearisation is delta' = -delta + M delta(t - tau) with
    M = diag(f'(z*)) W. Returns a dict with the equilibrium, M, mu, omega,
    tau (one row per eigenvalue), and tau_c, omega_c: the first Hopf point.
    """
    W = np.array([[p['a'], p['b']], [p['c'], p['d']]])
    I = np.array([p['P'], p['Q']])
    x = equilibrium(p)
    F = 1 / (1 + np.exp(-p['beta'] * (W @ x + I)))
    M = np.diag(p['beta'] * F * (1 - F)) @ W
    mu = np.linalg.eigvals(M)
    tau = np.full((len(mu), nj), np.nan)
    omega = np.full(len(mu), np.nan)
    for k, m in enumerate(mu):
        tau[k], omega[k] = hopf_delays_from_mu(m, nj)
    kc = np.nanargmin(tau[:, 0])
    return dict(xstar=x, M=M, mu=mu, omega=omega, tau=tau,
                tau_c=tau[kc, 0], omega_c=omega[kc])


def simulate(p, tau, T, history, m=None):
    """Solve the model with the Euler map (2.9), step h = tau/m.

    history = (u0, v0) is the constant history on [-tau, 0]. Because
    h = tau/m, the delayed state is exactly m steps back. By default m is
    chosen so that h is close to 5e-4. Returns arrays t, u, v starting at
    t = -tau.
    """
    if m is None:
        m = max(1, round(tau / 5e-4))
    h = tau / m
    N = round(T / h)
    a, b, c, d, P, Q, beta = (p[k] for k in ('a', 'b', 'c', 'd', 'P', 'Q', 'beta'))
    u = [float(history[0])] * (m + 1)       # entries 0..m: the history
    v = [float(history[1])] * (m + 1)
    for n in range(m, N + m):
        ud, vd = u[n - m], v[n - m]          # state at t - tau
        u.append(u[n] + h * (-u[n] + 1 / (1 + math.exp(-beta * (a * ud + b * vd + P)))))
        v.append(v[n] + h * (-v[n] + 1 / (1 + math.exp(-beta * (c * ud + d * vd + Q)))))
    t = (np.arange(N + m + 1) - m) * h
    return t, np.array(u), np.array(v)


def local_maxima(u):
    """Values and indices of the local maxima of a sampled signal."""
    u = np.asarray(u)
    idx = np.where((u[1:-1] > u[:-2]) & (u[1:-1] >= u[2:]))[0] + 1
    return u[idx], idx


def lyapunov_exponent(p, tau, T_transient, T_average, history, m=None):
    """Largest Lyapunov exponent, from two nearby solutions of the Euler map.

    The state of a DDE is the whole history on [t - tau, t], so the distance
    between the two solutions is measured over all m+1 stored values. Once
    per delay interval we record the distance s and pull the second solution
    back to distance d0; then lambda = sum(log(s/d0)) / (total time).

    Only the last m+1 states are kept, in a circular buffer: index i holds
    the current state, and index j = (i+1) mod (m+1) the state at t - tau,
    which is overwritten by the new state once it has been used.
    """
    if m is None:
        m = max(1, round(tau / 5e-4))
    h = tau / m
    a, b, c, d, P, Q, beta = (p[k] for k in ('a', 'b', 'c', 'd', 'P', 'Q', 'beta'))
    L = m + 1

    def step(U, V, i, j):
        ud, vd = U[j], V[j]
        U[j] = U[i] + h * (-U[i] + 1 / (1 + math.exp(-beta * (a * ud + b * vd + P))))
        V[j] = V[i] + h * (-V[i] + 1 / (1 + math.exp(-beta * (c * ud + d * vd + Q))))

    # 1. transient
    U = [float(history[0])] * L
    V = [float(history[1])] * L
    i = m
    for _ in range(round(T_transient / h)):
        j = (i + 1) % L
        step(U, V, i, j)
        i = j

    # 2. two nearby solutions, renormalised once per delay interval
    d0 = 1e-8
    shift = d0 / math.sqrt(2 * L)
    U2 = [x + shift for x in U]
    V2 = [x + shift for x in V]
    S = 0.0
    N = round(T_average / h)
    for n in range(1, N + 1):
        j = (i + 1) % L
        step(U, V, i, j)
        step(U2, V2, i, j)
        i = j
        if n % m == 0:
            s = math.sqrt(sum((x2 - x) ** 2 for x, x2 in zip(U, U2)) +
                          sum((y2 - y) ** 2 for y, y2 in zip(V, V2)))
            S += math.log(s / d0)
            U2 = [x + (x2 - x) * d0 / s for x, x2 in zip(U, U2)]
            V2 = [y + (y2 - y) * d0 / s for y, y2 in zip(V, V2)]
    return S / (N * h)

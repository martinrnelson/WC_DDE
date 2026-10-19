function lambda = lyapunov_exponent(p, tau, Ttrans, Tavg, history, m)
%LYAPUNOV_EXPONENT  Largest Lyapunov exponent of the Wilson-Cowan DDE.
%
%   lambda = lyapunov_exponent(p, tau, Ttrans, Tavg, history, m)
%
%   p        parameter structure from wc_parameters
%   tau      the delay
%   Ttrans   length of the transient that is discarded
%   Tavg     length of the time over which the exponent is averaged
%   history  [u; v], the constant history on [-tau, 0]
%   m        steps per delay interval (default: h close to 5e-4)
%
%   Method (Section 6, Listing 6). Two nearby solutions are advanced with
%   the Euler map (2.9). The state of a DDE at time t is the whole function
%   on [t - tau, t], so their distance is measured over all m+1 stored
%   values, not just the current one. Once per delay interval the distance
%   s is recorded and the second solution is pulled back to distance d0
%   along the same direction. Then
%
%       lambda = (sum of log(s/d0)) / (total time).
%
%   lambda > 0 means nearby solutions separate exponentially (chaos);
%   lambda = 0 is typical of periodic and quasi-periodic solutions.
%
%   Only the last m+1 states are needed, so they are kept in a circular
%   buffer X: column i holds the current state and column j = i+1 (mod
%   m+1) holds the state at t - tau, which is overwritten by the new state
%   as soon as it has been used.

if nargin < 6
    m = max(1, round(tau / 5e-4));
end
h = tau / m;
W = [p.a, p.b; p.c, p.d];
I = [p.P; p.Q];
f = @(z) 1 ./ (1 + exp(-p.beta*z));
F = @(x, xd) -x + f(W*xd + I);          % right-hand side for x(t), x(t-tau)

% 1. Transient
X = repmat(history(:), 1, m+1);
i = m + 1;
for n = 1:round(Ttrans / h)
    j = mod(i, m+1) + 1;
    X(:, j) = X(:, i) + h * F(X(:, i), X(:, j));
    i = j;
end

% 2. Follow two nearby solutions, renormalising once per delay interval
d0 = 1e-8;
Y  = X + d0 / sqrt(numel(X));           % perturbed copy at distance d0
S  = 0;
N  = round(Tavg / h);
for n = 1:N
    j = mod(i, m+1) + 1;
    X(:, j) = X(:, i) + h * F(X(:, i), X(:, j));
    Y(:, j) = Y(:, i) + h * F(Y(:, i), Y(:, j));
    i = j;
    if mod(n, m) == 0
        D = Y - X;
        s = norm(D(:));
        S = S + log(s / d0);
        Y = X + (d0 / s) * D;
    end
end
lambda = S / (N * h);
end

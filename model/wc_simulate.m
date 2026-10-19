function [t, u, v] = wc_simulate(p, tau, T, history, m)
%WC_SIMULATE  Solve the Wilson-Cowan DDE with the Euler map (2.9).
%
%   [t, u, v] = wc_simulate(p, tau, T, history)
%   [t, u, v] = wc_simulate(p, tau, T, history, m)
%
%   p        parameter structure from wc_parameters
%   tau      the delay
%   T        final time
%   history  [u; v], the constant value of the solution on [-tau, 0]
%   m        number of steps per delay interval, so the step is h = tau/m
%            (default: m chosen so that h is close to 5e-4)
%
%   The outputs t, u, v are row vectors, starting at t = -tau.
%
%   Because h = tau/m, the delayed state x(t - tau) is always exactly m
%   steps back, x(:, n-m), and no interpolation is needed. The first m+1
%   columns of x hold the history function.
%
%   The Euler map is first-order accurate: halving h roughly halves the
%   error. For the rich parameter set (beta = 60) the solutions change
%   quickly, so h must be small; h = 5e-4 reproduces the bifurcation
%   values found with DDE-BIFTOOL.

if nargin < 5
    m = max(1, round(tau / 5e-4));
end
h = tau / m;
N = round(T / h);

W = [p.a, p.b; p.c, p.d];
I = [p.P; p.Q];
f = @(z) 1 ./ (1 + exp(-p.beta*z));

x = repmat(history(:), 1, N + m + 1);       % columns 1..m+1: the history
for n = m+1 : N+m
    x(:, n+1) = x(:, n) + h * ( -x(:, n) + f(W*x(:, n-m) + I) );
end

t = ((1:N+m+1) - m - 1) * h;
u = x(1, :);
v = x(2, :);
end

%LYAPDDE  Largest Lyapunov exponent of the Wilson-Cowan DDE
%
% Two nearby solutions are advanced with the Euler map. The state of a
% DDE is the whole history on [t-tau, t], so the separation is measured
% over all m+1 stored values and renormalised once per delay interval.

a = -6; b = 2.5; P = 0.2; beta = 60; tau = 0.1;   % also try tau = 0.068
f = @(z) 1./(1 + exp(-beta*z));
F = @(x,xd) -x + f([a b; b a]*xd + P);            % rhs for x(t), x(t-tau)
m = 200; h = tau/m;

X = repmat([0.3; 0.1], 1, m+1);   % circular buffer of the last m+1 states
i = m+1;                          % column holding the current state
for n = 1:round(100/h)            % discard a transient of length 100
    j = mod(i, m+1) + 1;          % oldest column = state at t - tau
    X(:,j) = X(:,i) + h*F(X(:,i), X(:,j));
    i = j;
end

d0 = 1e-8;  Y = X + d0/sqrt(numel(X));   % perturbed copy, distance d0
S = 0;  N = round(200/h);
for n = 1:N
    j = mod(i, m+1) + 1;
    X(:,j) = X(:,i) + h*F(X(:,i), X(:,j));
    Y(:,j) = Y(:,i) + h*F(Y(:,i), Y(:,j));
    i = j;
    if mod(n, m) == 0             % once per delay interval
        D = Y - X;  s = norm(D(:));
        S = S + log(s/d0);  Y = X + (d0/s)*D;
    end
end
lambda = S/(N*h)

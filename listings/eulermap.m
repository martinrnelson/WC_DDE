%EULERMAP  Euler map (2.9) for the Wilson-Cowan model (3.1)
%
% The step is h = tau/m, so the delayed state x(t-tau) is always exactly
% m steps back and no interpolation is needed. Constant history on [-tau,0].

a = 10; b = -10; c = 10; d = 2; P = -2; Q = -4; beta = 1;
f = @(z) 1./(1 + exp(-beta*z));
F = @(x,xd) -x + f([a b; c d]*xd + [P; Q]);  % x = [u; v] at t and t - tau

tau = 0.5; T = 400;                     % also try tau = 0.05 and 0.12
m = 100; h = tau/m; N = round(T/h);     % m steps per delay interval
x = repmat([0.19644; 0.13377], 1, N+m+1);   % columns 1:m+1 = history
for n = m+1:N+m
    x(:,n+1) = x(:,n) + h*F(x(:,n), x(:,n-m));
end
t = ((1:N+m+1) - m - 1)*h;

plot(t, x(1,:), 'b-'), xlabel('t'), ylabel('u(t)')
late = x(1, t > T-50);
amp  = max(late) - min(late)            % oscillation amplitude of u

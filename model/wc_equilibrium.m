function xstar = wc_equilibrium(p)
%WC_EQUILIBRIUM  The equilibrium (u*, v*) of the Wilson-Cowan model.
%
%   xstar = wc_equilibrium(p) returns xstar = [u*; v*] for the parameter
%   structure p from wc_parameters.
%
%   At an equilibrium the solution is constant, so u(t - tau) = u(t) and
%   the delay drops out. The equilibrium therefore does not depend on tau
%   and solves
%
%       u = f(a u + b v + P),    v = f(c u + d v + Q).
%
%   Writing this as G(x) = x - f(W x + I) = 0, with W = [a b; c d] and
%   I = [P; Q], we solve it by Newton's method, x <- x - DG(x) \ G(x),
%   starting from the guess p.x0. The Jacobian uses f'(z) = beta f (1 - f).

W = [p.a, p.b; p.c, p.d];
I = [p.P; p.Q];
f = @(z) 1 ./ (1 + exp(-p.beta*z));

x = p.x0;
for iteration = 1:50
    F  = f(W*x + I);
    G  = x - F;                                    % we want G = 0
    DG = eye(2) - diag(p.beta * F .* (1 - F)) * W; % Jacobian of G
    dx = -DG \ G;
    x  = x + dx;
    if norm(dx) < 1e-13
        break
    end
end
xstar = x;
end

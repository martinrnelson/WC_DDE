function y = wc_rhs(xx, p, tied)
%WC_RHS  Single-delay Wilson-Cowan right-hand side for DDE-BIFTOOL.
%   xx : 2 x 2 x nvec array, xx(:,1,:) = state at t, xx(:,2,:) = state at t-tau
%   p  : parameters [a b c d P Q beta tau]
%   tied (optional, default false): if true, c is set equal to b, so that
%        parameter 2 controls both cross-couplings (used in part 3).
%   Returns 2 x 1 x nvec.
if nargin < 3, tied = false; end
a = p(1); b = p(2); c = p(3); d = p(4); P = p(5); Q = p(6); be = p(7);
if tied, c = b; end
ud = xx(1,2,:);  vd = xx(2,2,:);
z1 = a*ud + b*vd + P;
z2 = c*ud + d*vd + Q;
y = [ -xx(1,1,:) + 1./(1 + exp(-be*z1));
      -xx(2,1,:) + 1./(1 + exp(-be*z2)) ];
end

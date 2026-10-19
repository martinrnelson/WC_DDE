function J = wc_deri(xx, p, nx, np, v, tied)
%WC_DERI  Analytic derivatives of wc_rhs in DDE-BIFTOOL 'sys_deri' format.
%
%   nx : delay indices (0 = current state, 1 = delayed state)
%   np : parameter index (1..8 = a b c d P Q beta tau)
%   v  : vector for second derivatives
%
%   Cases (as required by DDE-BIFTOOL):
%     nx = k,      np = []  : df/dx_k                      (2 x 2)
%     nx = [],     np = j   : df/dp_j                      (2 x 1)
%     nx = k,      np = j   : d/dp_j (df/dx_k)             (2 x 2)
%     nx = [k l],  np = []  : d/dx_k (df/dx_l * v)         (2 x 2)
%   With a third (vector) dimension nvec when xx is 2 x 2 x nvec.
%
%   Uses f(z) = 1/(1+exp(-beta z)),  f' = beta f(1-f),
%   f'' = beta^2 f(1-f)(1-2f),  df/dbeta = z f(1-f).
%   The delay tau (parameter 8) does not appear explicitly in the
%   right-hand side, so all derivatives with respect to it are zero.
if nargin < 6, tied = false; end
a = p(1); b = p(2); c = p(3); d = p(4); P = p(5); Q = p(6); be = p(7);
if tied, c = b; end
W = [a b; c d];

nvec = size(xx, 3);
ud = reshape(xx(1,2,:), 1, nvec);
vd = reshape(xx(2,2,:), 1, nvec);
z  = [a*ud + b*vd + P; c*ud + d*vd + Q];      % 2 x nvec
F  = 1./(1 + exp(-be*z));
g  = be*F.*(1 - F);                           % f'(z)
h  = be^2*F.*(1 - F).*(1 - 2*F);              % f''(z)

if numel(nx) == 1 && isempty(np)
    %% first derivative w.r.t. state
    if nx == 0
        J = repmat(-eye(2), 1, 1, nvec);
    else
        J = rowscale(W, g, nvec);
    end

elseif isempty(nx) && numel(np) == 1
    %% first derivative w.r.t. a parameter
    if np == 7
        dfdp = z.*F.*(1 - F);
    elseif np >= 8
        dfdp = zeros(2, nvec);
    else
        dz = dpar(np, ud, vd, tied, nvec);
        dfdp = g.*dz;
    end
    J = reshape(dfdp, 2, 1, nvec);

elseif numel(nx) == 1 && numel(np) == 1
    %% mixed derivative d/dp (df/dx)
    if nx == 0 || np >= 8
        J = zeros(2, 2, nvec);
    else
        if np == 7
            dg = F.*(1 - F).*(1 + be*z.*(1 - 2*F));
            dW = zeros(2);
        else
            [dz, dW] = dpar(np, ud, vd, tied, nvec);
            dg = h.*dz;
        end
        J = rowscale(W, dg, nvec) + rowscale(dW, g, nvec);
    end

elseif numel(nx) == 2 && isempty(np)
    %% second derivative w.r.t. state, applied to v
    if all(nx == 1)
        v = reshape(v, 2, []);
        if size(v, 2) == 1, v = repmat(v, 1, nvec); end
        Wv = W*v;                                % 2 x nvec
        J = rowscale(W, h.*Wv, nvec);
    else
        J = zeros(2, 2, nvec);
    end

else
    error('wc_deri: requested derivative not implemented');
end
end

%% ------------------------------------------------------------------------
function J = rowscale(M, s, nvec)
% J(i,j,k) = M(i,j) * s(i,k)
J = zeros(2, 2, nvec);
for i = 1:2
    J(i,:,:) = reshape(M(i,:)' * s(i,:), 1, 2, nvec);
end
end

function [dz, dW] = dpar(k, ud, vd, tied, nvec)
% Derivatives of the sigmoid arguments z (2 x nvec) and of the weight
% matrix W with respect to parameter k (k = 1..6).
dz = zeros(2, nvec);  dW = zeros(2);
switch k
    case 1                                   % a
        dz(1,:) = ud;  dW(1,1) = 1;
    case 2                                   % b (and c, if tied)
        dz(1,:) = vd;  dW(1,2) = 1;
        if tied, dz(2,:) = ud;  dW(2,1) = 1; end
    case 3                                   % c (unused if tied)
        if ~tied, dz(2,:) = ud;  dW(2,1) = 1; end
    case 4                                   % d
        dz(2,:) = vd;  dW(2,2) = 1;
    case 5                                   % P
        dz(1,:) = 1;
    case 6                                   % Q
        dz(2,:) = 1;
end
end

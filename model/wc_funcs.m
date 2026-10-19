function funcs = wc_funcs(tie_bc)
%WC_FUNCS  The Wilson-Cowan model in the form DDE-BIFTOOL needs.
%
%   funcs = wc_funcs()       the model with parameters [a b c d P Q beta tau]
%   funcs = wc_funcs(true)   the same model, but with c set equal to b, so
%                            that changing parameter 2 changes both
%                            cross-couplings at once. This keeps the two
%                            populations identical in the two-parameter
%                            diagram (Figure 5).
%
%   DDE-BIFTOOL needs three things:
%     sys_rhs   the right-hand side. It is called as f(x, par), where
%               x(:,1) is the state at time t and x(:,2) the state at
%               time t - tau (one column per delay).
%     sys_tau   which entry of par is the delay.
%     sys_deri  derivatives of the right-hand side. If these are not
%               supplied, DDE-BIFTOOL approximates them by finite
%               differences; the exact formulas in wc_deri are faster and
%               make the Floquet multipliers more accurate.
%   'x_vectorized' tells DDE-BIFTOOL that wc_rhs and wc_deri accept many
%   states at once, which speeds up the computations.

if nargin < 1
    tie_bc = false;
end
ix = wc_index();
funcs = set_funcs( ...
    'sys_rhs',  @(x, par) wc_rhs(x, par, tie_bc), ...
    'sys_deri', @(x, par, nx, np, v) wc_deri(x, par, nx, np, v, tie_bc), ...
    'sys_tau',  @() ix.tau, ...
    'x_vectorized', true);
end

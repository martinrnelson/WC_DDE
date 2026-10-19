function p = wc_parameters(name)
%WC_PARAMETERS  The two parameter sets used in the paper (Table 1).
%
%   p = wc_parameters('gentle')   beta = 1: a single delay-induced Hopf
%                                 bifurcation (Figure 1)
%   p = wc_parameters('rich')     beta = 60: Hopf, torus, period doubling,
%                                 folds and chaos (Figures 2-7)
%
%   Both sets are taken from Coombes and Laing (2009). The model is
%
%     u'(t) = -u(t) + f( a u(t-tau) + b v(t-tau) + P )
%     v'(t) = -v(t) + f( c u(t-tau) + d v(t-tau) + Q ),
%
%   with f(z) = 1/(1 + exp(-beta z)). The output p is a structure with
%   fields a, b, c, d (connection weights), P, Q (external inputs), beta
%   (steepness of the sigmoid) and x0 (a starting guess for the
%   equilibrium, used by wc_equilibrium).
%
%   The delay tau is not stored here, because it is the parameter we vary.
%   DDE-BIFTOOL wants all parameters, including tau, in one vector: see
%   wc_parvec.

switch lower(name)
    case 'gentle'
        p.a = 10;   p.b = -10;  p.c = 10;   p.d = 2;
        p.P = -2;   p.Q = -4;   p.beta = 1;
        p.x0 = [0.19; 0.13];
    case 'rich'
        p.a = -6;   p.b = 2.5;  p.c = 2.5;  p.d = -6;
        p.P = 0.2;  p.Q = 0.2;  p.beta = 60;
        p.x0 = [0.07; 0.07];
    otherwise
        error('wc_parameters: unknown set ''%s''; use ''gentle'' or ''rich''.', name)
end
p.name = lower(name);
end

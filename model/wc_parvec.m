function par = wc_parvec(p, tau)
%WC_PARVEC  The parameters as the single vector DDE-BIFTOOL expects.
%
%   par = wc_parvec(p, tau) returns  par = [a b c d P Q beta tau],
%   where p comes from wc_parameters. The positions of tau and b in this
%   vector are given by wc_index.

par = [p.a, p.b, p.c, p.d, p.P, p.Q, p.beta, tau];
end

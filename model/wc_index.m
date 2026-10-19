function ix = wc_index()
%WC_INDEX  Positions of the parameters we vary, in the vector from wc_parvec.
%
%   ix.tau = 8   the delay (varied in every one-parameter computation)
%   ix.b   = 2   the cross-coupling b (varied, together with tau, in the
%                two-parameter diagram of Figure 5)
%
%   DDE-BIFTOOL refers to parameters by their position, so these numbers
%   appear in calls such as  'contpar', ix.tau.

ix.tau = 8;
ix.b   = 2;
end

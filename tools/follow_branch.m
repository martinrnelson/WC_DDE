function [idx, ifold] = follow_branch(tau, nfolds)
%FOLLOW_BRANCH  Follow a branch from its first point through nfolds folds.
%
%   [idx, ifold] = follow_branch(tau, nfolds)
%
%   tau     the delay at each point of a branch, in the order computed
%   nfolds  how many folds (turning points in tau) to pass
%
%   idx     indices of the points from the start of the branch up to, but
%           not including, the next fold after the first nfolds
%   ifold   indices (into tau) of the folds that were passed
%
%   A fold is where tau stops decreasing and starts increasing, or vice
%   versa: the branch turns back on itself. We use this to show only part
%   of a branch that folds many times in a narrow range of tau.

tau   = tau(:)';
dtau  = diff(tau);
turn  = find(sign(dtau(1:end-1)) ~= sign(dtau(2:end))) + 1;   % turning points
ifold = turn(1:min(nfolds, numel(turn)));
if numel(turn) > nfolds
    idx = 1:turn(nfolds+1) - 1;
else
    idx = 1:numel(tau);
end
end

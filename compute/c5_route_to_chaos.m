%% C5  Rich parameter set: orbit diagram, Lyapunov exponents, example solutions
%
% Paper:    Section 6; Figures 6 and 7
% Produces: data/route_to_chaos.mat, used by fig6_route_to_chaos.m and
%           fig7_chaos_examples.m
% Needs:    nothing beyond this project (no DDE-BIFTOOL)
% Run time: about 15 minutes
%
% Everything here is direct simulation with the Euler map (wc_simulate),
% so it shows what a solution actually does, including quasi-periodic and
% chaotic behaviour that continuation cannot follow.
%
% Steps
%   1. Orbit diagram: for each tau, simulate, discard the transient and
%      record the local maxima of u.
%   2. Largest Lyapunov exponent for each tau.
%   3. Four example solutions for Figure 7.
%   4. Save.

clear; close all
p       = wc_parameters('rich');
history = [0.3; 0.1];      % u and v must differ, or the solution stays on u = v
                           % and never shows the anti-phase dynamics

%% 1. Orbit diagram
% A periodic solution gives a few distinct maxima, a quasi-periodic one
% fills an interval, and a chaotic one a scattered set.
diagram.tau = 0.040 : 5e-4 : 0.130;
diagram.T_transient = 200;   diagram.T_record = 50;
diagram.points_tau = [];  diagram.points_umax = [];
tic
for k = 1:numel(diagram.tau)
    tau = diagram.tau(k);
    [t, u] = wc_simulate(p, tau, diagram.T_transient + diagram.T_record, history);
    peaks = local_maxima(u(t > diagram.T_transient));
    if isempty(peaks)                    % a steady state has no maxima:
        peaks = u(end);                  % record its value instead
    end
    diagram.points_tau  = [diagram.points_tau,  tau*ones(1, numel(peaks))];
    diagram.points_umax = [diagram.points_umax, peaks];
end
fprintf('Orbit diagram: %d delays in %.0f s\n', numel(diagram.tau), toc);

%% 2. Largest Lyapunov exponent
lyap.tau    = 0.040 : 5e-4 : 0.130;   % fine enough to locate the onset of chaos
lyap.lambda = zeros(size(lyap.tau));
tic
for k = 1:numel(lyap.tau)
    lyap.lambda(k) = lyapunov_exponent(p, lyap.tau(k), 200, 200, history);
end
fprintf('Lyapunov exponents: %d delays in %.0f s\n', numel(lyap.tau), toc);

% where the exponent first becomes clearly positive
lyap.tau_chaos = lyap.tau(find(lyap.lambda > 0.1, 1));
fprintf('lambda first exceeds 0.1 at tau = %.3f\n', lyap.tau_chaos);

%% 3. Example solutions for Figure 7
examples.tau   = [0.060, 0.068, 0.0735, 0.100];
examples.label = {'periodic', 'quasi-periodic', 'mode-locked', 'chaotic'};
T = 400;
for k = 1:numel(examples.tau)
    [t, u, v] = wc_simulate(p, examples.tau(k), T, history);
    window = t > T - 2;                      % short window for the time series
    examples.series(k) = struct('t', t(window) - (T - 2), 'u', u(window), 'v', v(window));
    window = t > T - 40;                     % longer window for the (u,v) plane
    examples.phase(k)  = struct('u', u(window), 'v', v(window));
    examples.maxima{k} = local_maxima(u(t > 100));   % for the return map
end

%% 4. Save
save(wc_file('data', 'route_to_chaos.mat'), 'p', 'history', 'diagram', ...
     'lyap', 'examples');
fprintf('Saved data/route_to_chaos.mat\n');

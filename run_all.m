%% RUN_ALL  Recompute every result and redraw every figure in the paper.
%
% Total run time: roughly 45 minutes, almost all of it in c3, c4 and c5.
% To redraw the figures from results you already have (a few seconds),
% run only the second block. The figure scripts do not need DDE-BIFTOOL.
%
% (Each script starts with 'clear', so this file calls them one by one
% rather than in a loop.)

wc_setup

%% 1. Computations: each saves its results in data/
c1_gentle_hopf              % Figure 1                       ~1 min
c2_rich_hopf_torus          % Figures 2(a), 3(a), 4(a)       a few min
c3_rich_period_doubling     % Figures 2(b,c), 3(b), 4(b,c)   ~10 min
c4_two_parameter            % Figure 5 (needs c2 and c3)     ~15 min
c5_route_to_chaos           % Figures 6 and 7 (no DDE-BIFTOOL)  ~15 min

%% 2. Figures: each loads from data/ and saves a PDF in pdf/
fig1_gentle_hopf
fig2_rich_bifurcations
fig3_floquet_multipliers
fig4_orbit_profiles
fig5_two_parameter
fig6_route_to_chaos
fig7_chaos_examples

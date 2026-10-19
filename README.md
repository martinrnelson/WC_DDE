# Code for the SIAM Review paper on delay differential equations in neuroscience

Code for every listing and figure in the paper by I. Pinder, M. R. Nelson and J. J. Crofts. It studies the Wilson–Cowan model with one delay τ,

  u′(t) = −u(t) + f(a u(t−τ) + b v(t−τ) + P),
  v′(t) = −v(t) + f(c u(t−τ) + d v(t−τ) + Q),    f(z) = 1/(1 + e^(−βz)),

and shows how a delay creates oscillations, and then torus bifurcations, period doubling and chaos.

The bifurcation analysis (Figures 1–5) is done in MATLAB with the continuation package **DDE-BIFTOOL**. The simulations behind Figures 6–7 use the Euler map and need only MATLAB, and are also provided in **Python** (see "Python versions" below).

The code is written to be read. Every script starts with a short description of what it does, which part of the paper it belongs to, what it needs and how long it takes, and is then divided into numbered steps that follow the paper.

## What you need

- MATLAB R2020a or later (for `exportgraphics`). No toolboxes are needed.
- [DDE-BIFTOOL](https://sourceforge.net/projects/ddebiftool/) version 3.x, free to download. It does the continuation in `c1`–`c4` and Listings 4–5. The simulations (`c5`, Listings 1–3 and 6) do not need it, and nor do the figure scripts, which only plot results already computed.

## Getting started

1. Open `wc_setup.m` and change the line `ddebiftool_folder = ...` to the folder where you put DDE-BIFTOOL.
2. In MATLAB, go to this folder and type `wc_setup`.
3. Run a listing, for example `hopfdelays`, or a figure script, for example `fig1_gentle_hopf`.

`run_all` recomputes everything and redraws every figure.

## How the code is organised

| Folder | What is in it |
|---|---|
| `listings/` | The short programs printed in the paper (Listings 1–6). Each one runs on its own. |
| `model/` | The model: parameter sets, right-hand side, equilibrium, exact Hopf delays, Euler-map solver. |
| `compute/` | The computations behind the figures. Each saves its results in `data/`. |
| `figures/` | One script per figure. Each loads from `data/` and saves a PDF in `pdf/`. |
| `tools/` | Small helpers: plotting styles, peak finding, Lyapunov exponent. |

The computations and the plotting are kept apart on purpose. The computations can take minutes, while a figure can then be redrawn, restyled or changed in seconds. The saved data contain only ordinary arrays, so you can plot them without DDE-BIFTOOL.

## Which script makes which figure

| Figure | What it shows | Computed by | Drawn by |
|---|---|---|---|
| 1 | Delay-induced Hopf bifurcation, gentle set | `c1_gentle_hopf` | `fig1_gentle_hopf` |
| 2 | Hopf → torus; period doubling → folds, rich set | `c2_rich_hopf_torus`, `c3_rich_period_doubling` | `fig2_rich_bifurcations` |
| 3 | Floquet multipliers past the torus and period doubling | `c2`, `c3` | `fig3_floquet_multipliers` |
| 4 | Three periodic orbits and their symmetry | `c2`, `c3` | `fig4_orbit_profiles` |
| 5 | Two-parameter diagram in (b, τ) | `c4_two_parameter` (needs `c2`, `c3`) | `fig5_two_parameter` |
| 6 | Orbit diagram and Lyapunov exponent | `c5_route_to_chaos` | `fig6_route_to_chaos` |
| 7 | Periodic, quasi-periodic, mode-locked and chaotic solutions | `c5_route_to_chaos` | `fig7_chaos_examples` |

Run times on a laptop: `c1` about a minute, `c2` a few minutes, `c3` about 10 minutes, `c4` about 15 minutes, `c5` about 15 minutes.

## The two parameter sets

Both come from Coombes and Laing (2009) and are defined in `model/wc_parameters.m` (Table 1 of the paper).

| Set | a | b | c | d | P | Q | β | First Hopf delay |
|---|---|---|---|---|---|---|---|---|
| `'gentle'` | 10 | −10 | 10 | 2 | −2 | −4 | 1 | τ ≈ 0.108 |
| `'rich'` | −6 | 2.5 | 2.5 | −6 | 0.2 | 0.2 | 60 | τ ≈ 0.0486 |

The gentle set has a single, clean Hopf bifurcation and is the place to start. The rich set has a steep sigmoid and two identical populations, and its oscillations go on to undergo torus bifurcations, period doubling, folds and chaos.

## Conventions in the figures

- Solid lines are stable and dashed lines unstable.
- Colours mark the type of object: equilibria black, periodic orbits blue, period-doubled orbits orange, Hopf vermillion, torus purple, period doubling green, folds sky blue.
- Markers for bifurcation points: ○ Hopf, ◇ torus, □ period doubling, △ fold.

## Things to try

- In `hopfdelays`, change `b` and watch the critical delay move. Which mode, in-phase or anti-phase, loses stability first?
- In `eulermap`, increase `m` (smaller steps) and compare the oscillation amplitude with the `dde23` result from `wc_dde23`.
- In `c2_rich_hopf_torus`, change `'radius', 1e-4` back to the default `0.01` and see why the small starting orbit matters.
- After any continuation, look at the printed **trivial multiplier error**. One Floquet multiplier of a periodic orbit is exactly 1 in theory. If the computed value is not close to 1, the other multipliers cannot be trusted either, so increase `'intervals'` and rerun.
- In `c5_route_to_chaos`, try a history with u = v, for example `[0.3; 0.3]`. The model maps u = v to itself, so the solution never leaves that line: at τ = 0.1, where the dynamics are chaotic, it settles on the equilibrium instead, because only the anti-phase mode is unstable there. Why does `[0.3; 0.1]` behave so differently?

## Python versions

The `python/` folder has Python versions of everything that does not need DDE-BIFTOOL. They need Python 3.8 or later with NumPy, SciPy and Matplotlib. Run them from inside `python/`, e.g. `python hopfdelays.py`.

| File | What it does | Run time |
|---|---|---|
| `hopfdelays.py`, `eulermap.py`, `lyapdde.py` | Listings 1, 2 and 6 | seconds |
| `wcmodel.py` | The model: parameters, equilibrium, exact Hopf delays, Euler map, Lyapunov exponent | — |
| `hopf_curves.py` | The exact Hopf curves and double Hopf points of Figure 5 | seconds |
| `route_to_chaos.py` | Orbit diagram, Lyapunov exponents and example solutions (as `c5`) | about 15 minutes |
| `fig6_route_to_chaos.py`, `fig7_chaos_examples.py` | Figures 6 and 7 from the results of `route_to_chaos.py` | seconds |

We know of no Python package that does the continuation of DDEs that DDE-BIFTOOL does (branches of periodic orbits, Floquet multipliers, torus and period-doubling curves), so Listings 4–5 and Figures 1–5 are MATLAB only. DDE-BIFTOOL itself is free; in Julia, `DDEBifurcationKit.jl` offers bifurcation analysis of DDEs.

## References

- S. Coombes and C. Laing, *Delays in activity-based neural networks*, Phil. Trans. R. Soc. A 367 (2009) 1117–1129.
- K. Engelborghs, T. Luzyanina and D. Roose, *Numerical bifurcation analysis of delay differential equations using DDE-BIFTOOL*, ACM Trans. Math. Softw. 28 (2002) 1–21.
- J. Sieber, K. Engelborghs, T. Luzyanina, G. Samaey and D. Roose, *DDE-BIFTOOL manual*, arXiv:1406.7144.

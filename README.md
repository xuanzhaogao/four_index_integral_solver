# Solver for four-index integral

Research code and results for the boundary-integral four-index (ERI) solver paper. The solver
itself is `BoundaryIntegral.jl` (`~/codes/BoundaryIntegral.jl`); this repo holds the studies built
on it.

- `codes/article/` — everything behind the paper's numerical results (§5.1–§5.3); start with its
  `README.md`, which maps each figure and table to its script.
- `codes/fig_gen/` — figures for the method sections (§2–4) and the shared `fig_style.jl`.
- `codes/*` (others) — earlier method studies (box2d/box3d, near_eval, graphene, multi_rhs, …).
- `article/` — an old partial LaTeX draft; the current manuscript lives outside this repo.
- `notes/`, `results/` — derivation notes and early results write-ups.

# Numerical results of the article (§5)

Everything that produces a figure or table in §5 of the BIE four-index paper. The manuscript
itself lives on the Mac (`~/Articles/four_indices_bie`), not in this repo. Figures for the method
sections (§2–4) are made by `../fig_gen/`; its `fig_style.jl` is shared by every plot here.

```
article/
├── Project.toml          one environment for everything below (BoundaryIntegral, TKM3D = dev paths)
├── common/               Harness.jl / Lite.jl (§5.1, §5.2 drivers), _provenance.sh (sbatch helper)
├── sec51_convergence/    §5.1  convergence and dielectric contrast
├── sec52_performance/    §5.2  batched (multi-RHS) interface solve: thread + K scaling
├── sec53_lattice/        §5.3  lattice-scale four-index tensor, Si|SiO₂ heterojunction
└── docs/                 historical planning notes of the former numerical_results/ tree
```

## Environment

```sh
cd codes/article
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run every script with `--project=codes/article` (or `--project=.` from here). The §5.3 sbatch
scripts `cd` into `sec53_lattice/` and use a bare `--project`, which finds this `Project.toml`.

## Which script makes what

Data trees and figure folders carry a *tag* suffix (`data<TAG>/`, `figs<TAG>/`), selected with
`RERUN_TAG`. The empty tag is the data behind the **submitted** manuscript, kept for comparison;
the tags below are the ones the current text uses.

| § | output | script (in its section dir) | tag / variant |
|---|---|---|---|
| 5.1 | convergence figure: geometry, error vs N, GMRES iterations vs N (`figs_v2/fig61_fig1_combined.pdf`) | `plot_fig1_combined.jl` | `RERUN_TAG=_v2` |
| 5.1 | dielectric-contrast figure and Table 1 data (`figs_v2/fig63_contrast.pdf`) | `plot_contrast.jl`, `contrast/scripts/analyze_contrast.jl` | `RERUN_TAG=_v2` |
| 5.1 | the `_v2` data itself | `jobscripts/01_fig8_sweep.sbatch`, `jobscripts/02_tab1_contrast.sbatch` | writes `data_v2/` |
| 5.2 | Tables 2 and 3 (`figs_sec53seq/table_*.tex`) | `scripts/article_tables.jl` | `RERUN_TAG=_sec53seq EXTRA_TAG=_sec53ext THREADS_TAG=_sec53` (the K > 46 rows come from `data_sec53ext/`, whose `.jls` records are untracked) |
| 5.2 | thread + K scaling figure (`figs_sec53/fig66_multicube_scaling.pdf`) | `plot_scaling.jl` | `RERUN_TAG=_sec53` |
| 5.2 | the `_sec53*` data | `slurm/run_multicube_*_sec53.sbatch` | see `PAPER_CHANGES_sec53.md` |
| 5.3 | \|V\| heatmap + l_ec convergence (`figs/fig_63_eps2.4_mo64.pdf`) | `scripts/plot_fig_63.jl` | `FIG63_VARIANT=eps2.4_mo64` |
| 5.3 | junction figure: U − U^bare/ε_slab along x and vs r (`figs/fig_junction_lattice_conv_l3_eps2.4_k46.pdf`) | `scripts/plot_junction_analysis.jl` (needs `scripts/bare_reference.jl` once) | campaign `lattice_conv_l3_eps2.4_k46` |
| 5.3 | the V tensor (on ceph, `/mnt/ceph/users/xgao1/four_index/<campaign>`) | `driver.jl` + `jobscripts/run_all.sbatch` | see `sec53_lattice/README.md` |

Notes per section:

- **§5.1** — `convergence/` and `contrast/` each hold `scripts/` and their `data*/` trees;
  `convergence/NOTES.md` records the sweep design.
- **§5.2** — `NOTES.md` is the experiment description; `PAPER_CHANGES_sec53.md` records how every
  number in Tables 2–3 was measured, including withdrawn claims. Untracked run logs of the old
  `numerical_results/` tree are in `logs_numerical_results/`.
- **§5.3** — `sec53_lattice/README.md` describes the multi-node campaign pipeline
  (prepare → solve → consolidate → eval → assemble). Figure PDFs/PNGs under `sec53_lattice/figs/`
  are gitignored and regenerated from the scripts and the tracked `.tsv` data.

Slurm job scripts are templates: submit them yourself. See
https://wiki.flatironinstitute.org/SCC/Software/Slurm.

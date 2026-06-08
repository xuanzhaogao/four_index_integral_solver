# multi_rhs — many-RHS / four-index BIE experiments

Scripts, benchmarks, and results for the **multi-RHS** extension of the dielectric BIE
four-index-integral solver (`BoundaryIntegral.jl`, `multi_rhs` branch). A center *i* couples to
many neighbor pair densities ρ_ij = φ_i·φ_j; these share the same operator A, so they are solved
together as a block (block GMRES + `nd=K` batched FMM) and contracted into Vijkl.

The Julia env `dev`s `BoundaryIntegral.jl` (`~/codes/BoundaryIntegral.jl`). Run scripts from this
directory with `julia --project`. `Manifest.toml` is gitignored — run `julia --project -e 'using
Pkg; Pkg.instantiate()'` once to build it (re-`dev`s BoundaryIntegral by path).

## Layout

| path | contents |
|---|---|
| `*.bie` | system inputs: graphene monolayer, 8-neighbor, full-cell, 2-shell (K=36) |
| `run_four_index_bie.jl` | full `.bie → Vijkl` pipeline for a center group |
| `multi_rhs_graphene_5x5x1.jl`, `compare_block_vs_single.jl` | block-vs-separate validation |
| `bench_fmm.jl`, `bench_fmm_uniform.jl` | FMM3D kernel sweep (OMP × nd × N): slab vs 3D-uniform |
| `bench_per_rhs.jl`, `bench_solve_breakdown.jl` | block-solve per-RHS sweep + instrumented FMM/matmul/glue/GMRES breakdown |
| `plot_*.jl` | figures (read `data/*.csv`, write `figs/*.png`) |
| `jobscripts/` | Slurm batch (`ccm`/`genoa`, build-once + array) and shell drivers |
| `data/` | parsed results (`*.csv`), incl. `fmm_phase_breakdown.csv`; `data/raw/` = raw Slurm `.out` logs (gitignored) |
| `figs/` | generated figures |
| `docs/` | [`multi_rhs_report.md`](docs/multi_rhs_report.md), [`performance_findings.md`](docs/performance_findings.md) |
| `tools/` | `nd_timing.f90` — FMM3D phase-timing harness (links an FMM3D checkout; see header) |

## Key results (see `docs/`)

- **Block multi-RHS**: per-RHS solve amortization **~4×** (K=1→36), tracking the FMM (94–97 % of the
  solve); block GMRES does *not* cut iterations (`iter=nmv=7` ∀K) — the batched `nd=K` FMM is the
  only lever. (An earlier "18×" was an unpinned-threads artifact; corrected.)
- **Why `nd` batching helps** (FMM3D source + `tools/nd_timing.f90` phase breakdown): each kernel
  computes the per-interaction geometry/operators once (`1/r`, Legendre recurrences, M2L operators)
  and loops `do idim=1,nd` for cheap multiply-adds — ~84 % of an `lfmm3d` call is this shared work
  (M2L + P2P dominate), only ~16 % scales with `nd`. Tree construction is ~2 %.
- **FMM thread scaling**: saturates ~32–64 cores (96 ≈ 64); 1→96 ≈ 18–20×.
- **Geometry**: the thin slab interface costs up to **~7×** vs a 3D-uniform cloud (≈2.75× more
  octree work × ≈2.5× worse parallel scaling).

> Benchmarks: pin **both** `OMP_NUM_THREADS` and `OPENBLAS_NUM_THREADS`. Use the julia **1.12**
> binary (the manifest is resolved with 1.12; 1.10 fails to load it). Submit to `ccm` / `-C genoa`.

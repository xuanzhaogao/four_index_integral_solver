# Per-stage runtime — production pipeline (PrecomputedVolumeField + cache_fft)

Date: 2026-06-15. worker7001 (idle, 96 threads), warm. Real graphene monolayer
pz orbital (k_323201), calibrated parameters: L = 90, Lz = 3.35, eps_in = 3.5,
source_tol/rhs_tol 1e-3, lhs_tol 1e-5, gmres 1e-5, n_quad 6, edge_level 4,
correct_edges = true. Script: `scripts/bench_full_field.jl`, raw
`data/raw/bench_full_field.jls`.

Problem: 960,768 interface points, 355,862 source points, GMRES 8 iters,
residual 7.1e-6, u_onsite = 9.88 eV. Peak RSS 29.4 GB.

## Per-stage runtime

| phase | stage | runtime (s) |
|---|---|---:|
| load | XSF read + truncation + screening | 3.25 |
| precompute | field construction (cache_fft) | 8.82 |
| precompute | interface build (RHS-adaptive) | 6.76 |
| precompute | LHS assembly (near corrections) | 13.56 |
| solve | RHS assembly (field) | 2.18 |
| solve | GMRES solve (8 iters) | 22.21 |
| eval | volume potential u_int (field) | 0.12 |
| eval | scattered potential (FMM + hcubature) | 5.38 |

## Phase totals

| phase | runtime (s) | share |
|---|---:|---:|
| load | 3.25 | 5% |
| precompute (construction + build + LHS) | 29.14 | 47% |
| solve (RHS + GMRES) | 24.39 | 39% |
| eval (u_int + scatter) | 5.50 | 9% |
| **total** | **62.28** | 100% |

(end-to-end wall including startup + warm-up: 67.0 s)

## Notes

- Top costs: GMRES (22.2 s), LHS assembly (13.6 s), field construction
  (8.8 s), interface build (6.8 s). The pipeline is solve/assembly-bound.
- The RHS-independent precompute (field construction + interface build + LHS
  assembly = 29.1 s) is paid once; each additional RHS on the same
  geometry/source family costs RHS + GMRES + eval ~= 30 s (GMRES-dominated).
- Peak RSS 29.4 GB is set by the cache_fft fine-grid storage. The standard
  field (cache_fft = false) peaks ~14 GB with u_int ~0.4 s instead of 0.12 s.

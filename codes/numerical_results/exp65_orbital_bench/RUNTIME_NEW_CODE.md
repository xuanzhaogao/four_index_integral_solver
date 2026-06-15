# Per-stage runtime — new code (PrecomputedVolumeField + cache_fft)

Date: 2026-06-15. worker7001 (idle, exclusive use), 96 threads, warm.
Real graphene monolayer pz orbital (k_323201), calibrated parameters
(L = 90, Lz = 3.35, eps_in = 3.5, source_tol/rhs_tol 1e-3, lhs_tol 1e-5,
gmres 1e-5, n_quad 6, edge_level 4, correct_edges = true). Identical
configuration and stages to the original `run_single_rhs.jl`, so the columns
are a clean old-vs-new comparison. Script: `scripts/bench_full_field.jl`,
raw `data/raw/bench_full_field.jls`.

Result: interface points 960,768 (unchanged), GMRES 8 iters, residual
7.1e-6, **u_onsite = 9.88 eV** (matches the original pipeline). Peak RSS
29.4 GB (cache_fft fine-grid storage; standard-field variant peaks ~14 GB).

## Per-stage table

| phase | stage | new code (s) | original (s) | note |
|---|---|---:|---:|---|
| load | XSF read + truncation + screening | 3.25 | 3.2 | unchanged (I/O) |
| precompute | field construction (cache_fft) | 8.82 | — | new (replaces nothing; amortizes below) |
| precompute | **interface build (field overload)** | **6.76** | **201.4** | volume field reused; old hit the 96-thread FFT pathology per depth |
| precompute | LHS assembly (near corrections) | 13.56 | 16.9 | unchanged (FMM3D surface op) |
| solve | RHS assembly (field) | 2.18 | 2.2 | field gradient at interface points |
| solve | GMRES solve (8 iters) | 22.21 | 23.1 | unchanged (FMM3D matvec) |
| eval | **volume potential u_int (field)** | **0.12** | **44.2** | cache_fft interp-only; old was a full ltkm3dc on the pathological grid |
| eval | scattered potential (FMM + hcubature) | 5.38 | 4.5 | unchanged (FMM3D surface op) |

## Phase totals

| phase | new (s) | original (s) | speedup |
|---|---:|---:|---:|
| load | 3.25 | 3.2 | 1.0x |
| precompute (construction + build + LHS) | 29.1 | 218.7 | 7.5x |
| solve (RHS + GMRES) | 24.4 | 25.3 | 1.0x |
| eval (u_int + scatter) | 5.5 | 48.6 | 8.8x |
| **TOTAL** | **62.3** | **295.8** | **4.7x** |

(end-to-end wall incl. startup/warm-up overhead: 67.0 s)

## Reading the result

- **The two stages that collapsed are exactly the two that used TKM/FINUFFT
  volume transforms** and therefore hit the 96-thread FFT pathology in the old
  code: the RHS-adaptive interface build (201 -> 6.8 s) and the u_int volume
  potential (44 -> 0.12 s). The new code avoids them by (a) reusing one
  precomputed spectrum across all mesh-refinement depths and the evaluation
  (field), and (b) cache_fft's interp-only read-out, and (c) the field's box
  dimensions landing on a thread-friendly FFT grid.
- **The unchanged stages are unchanged** (to within run noise): LHS assembly
  (13.6 vs 16.9), GMRES (22.2 vs 23.1), scattered potential (5.4 vs 4.5).
  These are FMM3D surface operations on the (identical) interface — never
  affected by the FINUFFT thread bug or the field work.
- **The pipeline is now solve/assembly-bound, not volume-bound.** The top
  costs are GMRES (22 s), LHS assembly (14 s), field construction (9 s),
  interface build (7 s). The bare volume work (build + u_int) that dominated
  the old run is now ~7 s of 62 s.

## Caveat on the 4.7x headline

This 4.7x end-to-end is against the original `run_single_rhs.jl` baseline,
whose interface build and u_int were inflated by the FINUFFT 96-thread
pathology (see THREAD_PATHOLOGY.md). Most of the win is dodging that bug
(the field box happens to land on a thread-friendly grid), not algorithmic
superiority. The genuinely structural, bug-independent gains are: the field's
reuse of one spectrum across mesh depths and evaluation, and cache_fft's
interp-only read-out. Against a hypothetical thread-fixed old pipeline the
end-to-end gain would be smaller (the build and u_int would already be a few
seconds each), but the new code reaches that good regime automatically and
robustly.

## Amortization note (multi-RHS)

field construction (8.8 s) + interface build (6.8 s) + LHS assembly (13.6 s)
= the RHS-independent precompute (~29 s) is paid once. Each additional RHS
(same geometry/source family) costs RHS + GMRES + eval ~= 30 s, dominated by
GMRES. The volume field and the LHS operator both amortize across RHSs.

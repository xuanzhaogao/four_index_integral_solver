# Report: the "slow ltkm3dc" was a FINUFFT 96-thread FFT pathology

Date: 2026-06-15. Author: exp65 investigation. Hardware: Flatiron Rusty,
AMD EPYC Genoa 96-core nodes (worker7001/7018/7020/7049/7072/7093/7095, and
ccm Slurm nodes). Julia 1.12.6.

## TL;DR

The graphene-monolayer four-index pipeline spent ~50 s per
volume-potential (`u_int`) evaluation and ~200 s per adaptive mesh build. The
cause was **not** the algorithm, the tolerance, the grid size, or node
contention — it was a **FFTW threaded-FFT pathology**: on certain 3D mode-grid
dimensions, FINUFFT's FFT is ~**80× slower at 96 threads than at 64**, while a
neighbouring grid is unaffected. The production call happened to land on a bad
grid and ran at the full 96 threads. Reproduced standalone with only FINUFFT
(`repro/`, `RESULT.md`).

## How it was found, and the false leads (recorded so we don't relitigate)

The single-RHS production benchmark (`run_single_rhs.jl`) showed the adaptive
mesh build at ~200 s (68% of the solve) and `u_int` at ~44 s. Decomposing it:

| hypothesis | test | verdict |
|---|---|---|
| Grid too big (TKM box covers sources∪targets) | print both grids | **false** — old grid (413,375,325)=5.0e7 is *smaller* than the field grid (441,401,353)=6.2e7 |
| FINUFFT auto-picking upsampfac 2.0 | force 1.25/2.0 | **false** — auto picks 1.25 (the fast one) at both tols |
| Tolerance (1e-3 vs 1e-4) | fixed-grid tol sweep (`bench_tol_isolation.jl`) | **false** — 2.63 s vs 3.33 s type-1, 0.39 s vs 0.39 s type-2; tol irrelevant (tighter is marginally *slower*) |
| Node contention | dedicated exclusive Slurm node | **false** — pathology reproduced clean; the earlier "sum of pieces > full call" was a type-unstable top-level-global timing artifact, not contention |
| Per-call type-1 + type-2 recomputation (reuse) | — | **partly true** — real but small (~3×); does not explain the ~50 s |

The decisive experiment (`diag_ltkm_internal.jl`, worker7001) instrumented
`_ltkm3dc_eval` internally and swept the FINUFFT thread count:

```
REAL ltkm3dc(pgt=1) full call, 96 threads:  54.3 s
box        nthreads   type-1    scale    type-2    sum
field        96        3.61      0.44     0.43      4.47
field         1        2.77      0.44     2.26      5.48
ltkm3dc      96       27.51      0.35    27.23     55.09   = the 54 s
ltkm3dc      16        0.48      0.36     0.42      1.26    43x faster
ltkm3dc       1        1.96      0.35     1.90      4.22
```

Sum ≈ full call (decomposition faithful). The full 54 s = type-1 + type-2,
both at 27 s, only at 96 threads, only on the `ltkm3dc` box grid.

## The isolated reproduction (`repro/`, depends only on FINUFFT)

FINUFFT 3.5.2 / finufft_jll 2.5.1, synthetic uniform-random points, N=356k,
tol 1e-3. FFT execution time (`finufft_exec`; plan time ~0):

| grid (modes) | 1 | 8 | 16 | 32 | 64 | **96** |
|---|--:|--:|--:|--:|--:|--:|
| (413,375,325) | 2.02 | 0.51 | 0.40 | 0.35 | 0.34 | **27.33** |
| (441,401,353) | 2.38 | 0.63 | 0.49 | 0.40 | 1.98 | **0.39** |

- The cliff is **specific to the full core count (96)** and to the grid
  dimensions: (413,375,325) is fine at ≤64 threads then jumps ~80×; the
  neighbouring (441,401,353) is fine everywhere including 96.
- It is the **threaded FFT execution**, not planning (FFTW ESTIMATE, plan ~0).
- Both grids are 2,3,5-smooth after FINUFFT's upsampling (nf ≈ next-smooth of
  ~1.25× the mode count), so it is a bad threaded-FFT *decomposition* of a
  specific nf, not an obviously "ugly" size.
- A related symptom: at 96 threads with many plan create/destroy cycles, FFTW
  can **segfault** in `spawnloop` (seen in `bench_eps1e4.jl`) — the same
  threading fragility manifesting as a crash rather than a slowdown.

## Why the field path looked 100×+ faster (honest re-baselining)

`PrecomputedVolumeField` builds its box from the source bbox + a margin, which
lands on (441,401,353) — a thread-friendly grid — and reuses one type-1. So
its apparent 106–137× `u_int` speedup was **mostly this thread bug** in the old
path, not the field design. Against a *correctly threaded* old path:

- `u_int`: old (thread-fixed) ~1.3 s vs field query ~0.4 s → **~3×** (reuse).
- build / multi-RHS: the field's amortization (type-1 once, reused across mesh
  depths and RHSs) remains a genuine structural win, but the headline 9.5–22×
  build ratios must be re-measured against a thread-fixed `ltkm3dc` baseline.

The field path is still worth keeping (reuse + it dodges the pathology by
construction), but the dramatic factors in the earlier report sections are
inflated by the FFTW bug and are flagged as such there.

## Does changing eps to 1e-4 speed up the old code? No.

The mode grid is fixed by `kmax = π/h` (source spacing) and the box extent —
both independent of the solver tolerance. eps 1e-4 lands on the same
(413,375,325) grid, hits the same 96-thread cliff, same ~50 s. The fixed-grid
tol sweep confirms 1e-3 ≈ 1e-4 (tighter slightly slower). Tolerance is not a
lever for this cost; grid dimensions × thread count is.

## Recommendations (NOT a global nthreads cap)

1. **Report upstream to FINUFFT/FFTW.** A 2,3,5-smooth 3D transform of ~5e7
   modes should not be ~80× slower at 96 vs 64 threads. `repro/` is a
   self-contained reproducer (FINUFFT-only, pinned Project/Manifest).
2. **Pad mode-grid dimensions to thread-friendly sizes at the call site.** This
   is what the field box does for free; it keeps full 96-thread throughput on
   the well-behaved sizes. Preferable to capping threads globally (which would
   penalize the many grids that thread fine).
3. Optionally, a per-transform `nthreads` heuristic that lowers threads only
   for flagged-bad dimensions.

## Artifacts

- `repro/repro_finufft_threads.jl`, `repro/Project.toml`, `repro/Manifest.toml`,
  `repro/RESULT.md` — standalone FINUFFT reproducer + result.
- `scripts/diag_ltkm_internal.jl` — in-situ `_ltkm3dc_eval` decomposition +
  thread sweep (the decisive experiment).
- `scripts/bench_tol_isolation.jl` — fixed-grid tol × upsampfac sweep (ruled
  out tol).
- `REPORT.md` — full exp65 benchmark report; its profiling section carries the
  corrected, thread-aware conclusions.

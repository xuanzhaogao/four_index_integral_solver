# Multi-RHS / FMM Performance Findings

Notes from benchmarking the dielectric BIE multi-RHS solver and its FMM3D kernel on the
monolayer graphene model (cRPA data `k_323201_nb_144_c_15`; slab L=90, Lz=3.35, ε_in=3.5).
All kernel/solve timings below are from **isolated, exclusive genoa nodes** (96 cores, 1.5 TB),
one node per thread config, so they are free of core-sharing contention:

- slab FMM sweep — Slurm array job `6486631` → `../data/fmm_bench.csv`
- 3D-uniform baseline — Slurm array job `6483176` → `../data/fmm_bench_uniform.csv`
- solve breakdown — pinned-thread run → `../data/solve_breakdown.csv`

Scripts under `../` (`bench_fmm*.jl`, `bench_per_rhs.jl`, `bench_solve_breakdown.jl`,
`plot_*`); figures under `DATAFIG/`.

## TL;DR

- The solve cost is **almost entirely `lfmm3d`** (94–97 % at 1 thread, 70–85 % at 16). The
  near-correction sparse matmul is **negligible** (<0.5 %); GMRES is small and fixed (~1–3 s).
- **Batching RHS (`nd`) gives a flat ~3.7× per-RHS FMM win** at any thread count (the tree is
  built once; `nd=36` costs ~9.7× `nd=1`, not 36×). The full block **solve inherits the same
  ~4×** — and *only* that: block GMRES does **not** save iterations (`iter = nmv = 7` for every K),
  so the batched FMM is its sole lever.
- **FMM thread scaling saturates ~32–64 cores and is flat (no real gain) at 96.** 1→96 speedup is
  N-dependent: ~18× at N≈34 k, ~20× at N≈460 k.
- **The thin-slab surface geometry costs up to ~7× vs a 3D-uniform cloud** at the same N — split
  as **~2.75× more work** (deeper octree) **× ~2.5× worse thread scaling**.
- **An earlier "~18× per-RHS block solve" figure was a thread-contamination artifact** (unpinned
  threads); the correct number is ~4×. See §5.

## Setup

| benchmark | system | N (unknowns) | nd | OMP | fmm_tol |
|---|---|---|---|---|---|
| FMM kernel (`bench_fmm.jl`) | slab interface, far point source, edge-refine lvl 0–3 | 33,696 / 93,312 / 215,136 / 461,376 | 1,4,16,36 | 1…96 | 1e-4 |
| uniform baseline (`bench_fmm_uniform.jl`) | 3D-uniform random points | same N | same | same | 1e-4 |
| solve breakdown (`bench_solve_breakdown.jl`) | 2-shell graphene K=36 | 219,456 (6096 panels) | block 1 / 36 | 1 / 16 (pinned) | 1e-4 |
| per-RHS solve (`bench_per_rhs.jl`) | 2-shell graphene, K=1..36 | 219,456 | block sweep | pinned | 1e-5 |

## 1. `nd` amortization (batched RHS) — FMM kernel

Per-RHS FMM time `t_fmm/nd` falls as the block grows and plateaus. The factor is **~3.7× and
essentially independent of thread count and N** (clean slab data, N=461k, OMP=1):

| nd | total t_fmm (s) | per-RHS (s) | amortization vs nd=1 |
|---:|---:|---:|---:|
| 1 | 19.10 | 19.10 | 1.00× |
| 4 | 31.76 | 7.94 | 2.41× |
| 16 | 87.46 | 5.47 | 3.49× |
| 36 | 185.03 | 5.14 | **3.72×** |

Mechanism: one FMM tree build (~75 % of a single `nd=1` call) is shared across the `nd` densities;
only the far-field/gradient evaluation (~25 %) scales with `nd`. So `t_fmm(nd=36) ≈ 9.7·t_fmm(1)`,
giving 36/9.7 ≈ 3.7× per RHS, and the ceiling as nd→∞ is ~4× (only the shared 75 % can ever be
amortized).

## 2. Thread scaling — FMM kernel (nd=36)

Per-RHS FMM time vs OMP_NUM_THREADS (s), `fmm_tol=1e-4` (see `../figs/fmm_bench.png`):

| thr | N=33,696 | N=93,312 | N=215,136 | N=461,376 |
|---:|---:|---:|---:|---:|
| 1 | 0.350 | 0.843 | 2.013 | 5.140 |
| 16 | 0.037 | 0.081 | 0.288 | 0.578 |
| 64 | 0.019 | 0.040 | 0.125 | 0.261 |
| 96 | 0.020 | 0.039 | 0.129 | 0.252 |

- Near-ideal (∝ 1/threads) to ~16 threads, then **saturates**; **96 threads ≈ 64 threads** (flat,
  no meaningful gain — and marginally worse for small N).
- 1→96 speedup grows with N: **~18× (N=34k), ~20× (N=461k)** — bigger problems parallelize better.
- **Operating point ≈ 32–64 threads; 96 buys nothing.**

## 3. Combined nd + threads (per-RHS, vs nd=1 / 1 thread)

`36·t(nd=1,1thr) / t(nd=36, T)` at N=461k — folds both levers together:

| threads | combined per-RHS speedup |
|---|---|
| 1 | 3.7× (nd only) |
| 16 | 33× |
| 32 | 52× |
| 64 | **73×** |
| 96 | 76× |

Decomposes as **(~3.7× nd) × (thread factor)**; the nd part is constant, the thread part saturates.
64 threads (73×) is the robust operating point; 96 (76×) is within noise of it.

## 4. Slab geometry penalty (slab vs 3D-uniform)

Throughput (M source·RHS evals/s, nd=36; `../figs/fmm_heatmap.png`):

| | peak throughput | best config |
|---|---|---|
| slab interface | **2.37 M/s** | N=93k, 96 thr |
| 3D uniform | **12.87 M/s** | N=461k, 96 thr |

At N=461k the slab is **2.75× slower at 1 thread and 7.0× slower at 96 threads**. The penalty is
*not* a flat factor — it factors cleanly into two compounding causes:

$$7.0\times\;(\text{omp=96})\;=\;\underbrace{2.75\times}_{\text{more work}}\;\times\;\underbrace{2.54\times}_{\text{worse scaling}}$$

1. **More work (~2.75×, present even single-threaded).** FMM3D's octree is built on a *cubic* root
   box, but the interface lives on a 90×90×3.35 slab — points fill only ~3.7 % of the cube's
   z-extent. The first ~`log₂(90/3.35)≈5` levels split empty z-space, so the tree refines like a
   2D quadtree (branch ~4, not 8) and goes ~2 levels **deeper** than a volume-filling cloud to
   reach the same leaf occupancy. More levels × more boxes per useful point ⇒ a larger FMM
   constant. Edge-refinement clustering on the faces adds multiscale non-uniformity on top.
2. **Worse parallel scaling (~2.5×).** The non-empty boxes are crammed in a thin sheet, so OpenMP
   load-balances poorly and interaction lists are ragged: thread scaling 1→96 at N=461k is **20×
   (slab) vs 52× (uniform)** — ~21 % vs ~54 % parallel efficiency.

This — geometry meeting the octree — is the primary cause of the low slab throughput, not the
kernel, the tolerance, or the near correction.

> **Anomaly:** the **uniform N=93,312** point is anomalously slow (throughput non-monotonic in N,
> lower than both N=34k and N=215k across all 8 thread counts), making its slab/uniform ratio < 1.
> This is an unlucky seeded random-point tree; that single point should be re-run with a different
> `Random.seed!` before it goes in a figure.

## 5. Solve cost breakdown (block GMRES, K=1/36, N=219k, pinned threads, 1e-4)

`../bench_solve_breakdown.jl`, instrumented matvec, `OMP_NUM_THREADS` pinned
(`../data/solve_breakdown.csv`):

| OMP | K | iter | nmv | total (s) | FMM | matmul | glue | GMRES | shares (FMM/mm/glue/gmres) | per-RHS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| 1 | 1 | 7 | 7 | 61.4 | 57.9 | 0.015 | 0.70 | 2.83 | 94 / 0 / 1 / 5 % | 61.4 |
| 1 | 36 | 7 | 7 | 552.6 | 537.4 | 0.232 | 12.7 | 2.25 | 97 / 0 / 2 / 0 % | **15.3** |
| 16 | 1 | 7 | 7 | 11.6 | 8.1 | 0.010 | 0.65 | 2.85 | 70 / 0 / 6 / 25 % | 11.6 |
| 16 | 36 | 7 | 7 | 91.6 | 77.6 | 0.230 | 12.7 | 1.06 | 85 / 0 / 14 / 1 % | **2.54** |

(glue = the per-matvec dense loops: build the `(K,n)` charge array + contract `grad·normal`.)

- **`lfmm3d` dominates (94–97 % at 1 thread; 70–85 % at 16)**; **matmul negligible** (≤0.23 s).
- **`nmv = 7` for K=1 and K=36 alike.** Block GMRES applies the operator the same number of times
  regardless of K — it does **not** save iterations. The whole per-RHS amortization comes from the
  batched FMM.
- **Per-RHS solve amortization (K=1→36): 4.0× (omp=1), 4.5× (omp=16)** — tracking the FMM (which
  scales 57.9→537.4 = 9.3× for nd=36, i.e. ~3.9× per RHS). Same ~4× ceiling as the bare kernel
  (§1).
- **GMRES overhead ≈ fixed (~1–3 s)** → big *share* (25 %) at K=1/16 threads (Amdahl), amortized to
  ~1 % at K=36.
- **The dense "glue" is single-threaded**; its share grows as the FMM parallelizes away: 2 % @
  omp=1 → 14 % @ omp=16 (→ ~26 % at 96 in an earlier run). It is the next thing to thread if
  pushing past ~32 cores on the full solve.

> **Correction (was §6).** An earlier *unpinned* `bench_per_rhs.jl` run reported a per-RHS block
> speedup of **~18×** that "tracked the ideal 1/K line," attributed to fewer iterations. **Both
> claims were wrong.** It was a **thread artifact**: with `OMP_NUM_THREADS` unset, small-K solves
> ran ~single-threaded while large-K solves picked up OpenBLAS/FMM threads — folding a ~4× thread
> speedup onto the K axis. Decisive check: at `omp=1` the K=36 FMM *alone* is 537 s, so that run's
> reported 127 s total for K=36 is impossible single-threaded. With threads pinned the honest
> number is **~4×**, and `iter = nmv = 7` for every K (no iteration savings).

## 6. Reconciling the speedup numbers (avoid conflating them)

| number | what it is |
|---|---|
| ~3.7× | `nd=36` vs `nd=1` FMM, **same threads** (per-RHS) |
| **~4×** | full block **solve** per-RHS, K=1→36 at **fixed** threads (tracks the FMM; `nmv=7` ∀K) |
| ~18–20× | FMM **1→96 threads** (small N → N=461k) |
| ~33× | `nd=36` @ 16 threads vs `nd=1` @ 1 thread (combined) |
| ~73× | `nd=36` @ 64 threads vs `nd=1` @ 1 thread (combined) |
| ~7× | slab-vs-uniform **geometry** penalty @ 96 threads, N=461k (= 2.75× work × 2.5× scaling) |
| ~~18.5×~~ | **RETRACTED** — earlier unpinned per_rhs run; thread artifact, not a real speedup |

## 7. Recommendations / next steps

- **Threads:** run FMM at **~32–64 cores**, not 96 (saturates; 96 ≈ 64).
- **Batch RHS:** always solve co-located densities as a block — a free ~3.7–4× per-RHS (FMM tree
  built once). Note this is the *only* block-GMRES benefit here; it does not cut iterations.
- **Always pin `OMP_NUM_THREADS` *and* `OPENBLAS_NUM_THREADS`** in benchmarks — the spurious 18×
  came from leaving them unset, so the FMM/OpenBLAS thread count drifted with problem size.
- **Thread the glue:** the single-threaded charge-assembly + `grad·normal` contraction grows from
  2 % to ~14–26 % of the solve as threads increase; threading/`@simd` it (or folding `grad·n` into
  the FMM dipole evaluation) recovers that at high core counts.
- **Geometry:** the thin-slab surface costs up to ~7× vs uniform; an octree better matched to the
  near-2D, large-aspect distribution (or a different near/far split) is where the remaining factor
  lives.
- **Tooling:** build the FMM operators **once** (`slurm_bench_fmm_build.slurm`) and fan out a
  job-array, one **exclusive** node per thread config (`slurm_bench_fmm_array.slurm`), partition
  `ccm` / `-C genoa`. Use the 1.12 julia binary directly (the manifest is resolved with 1.12; the
  `~/.julia/juliaup/julia-1.10*` binary fails to load it).

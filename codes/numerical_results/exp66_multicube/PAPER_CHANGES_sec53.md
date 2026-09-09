# Changes §5.2 needs after moving it onto the §5.3 system (2026-09-07)

§5.2 (performance of the batched interface solve — Tables 2, 3 and the scaling
figure) was benchmarked on a substrate the physics section does not use: two
90×90×90 cubes with a 90×90×9 slab at ε = 10, a generic demo value. It now runs
the **same system §5.3 reports** — the ×3 converged Si|SiO₂ heterojunction
(270³, ε 11.9 | 3.9) with the graphene slab at its cRPA ε = 2.4.

Every number below is measured. Jobs `6998896` (K sweep, K = 1..31), `6998897`
(thread sweep, 8 tasks) and `6999183` (K sweep extension, K = 37, 40, 46), all on
rocky8 genoa nodes, 2026-09-07. Data in `data_sec53/`; `data/` (submitted) and
`data_v2/` (July rerun) are untouched.

Only the substrate changed. Discretization (`n_quad = 6`, `edge_level = 4`,
`l_ec = 0.568`, `rhs_tol = 1e-3`, `lhs_tol = gmres_rtol = 1e-5`) and the pipeline
are identical to `_v2`, and the K = 1..31 points use the same cutoffs, so the
deltas below are attributable to the geometry and permittivity alone. Three
points (K = 37, 40, 46) are new — see Table 3.

---

## The system, and the number in the abstract

| | `_v2` | new |
|---|---|---|
| cubes | 90³ | **270³** |
| slab | 90×90×9, ε = 10 | **42.95×44.06×9, ε = 2.4** |
| surface unknowns *N* | 1,617,264 | **2,751,264** |
| screening ratio `V_vac/V` | 8.3948 | **2.4437** |

**The abstract's "1.6 million surface unknowns" must become 2.75 million.** That
figure is this benchmark's *N* (it matches no §5.3 campaign — `conv_l3` runs at
1,349,640 and `het3x` at 651,420), as already noted in
`rerun_2026_07/PAPER_CHANGES.md`.

*N* = 2,751,264 at **every** K from 1 to 46 — the slab is sized from the largest
cluster and held fixed, so the multi-RHS study is not confounded by a changing
interface. This matters for the memory finding below: the growth in resident
memory is entirely the K right-hand sides, not a growing discretization.

---

## Table 2 — strong scaling (K = 1)

Reconstructed as `t_precompute + t_RHS + t_1iter · niter_ref` with
`niter_ref = 12` (was 14), per `plot_scripts/plot_scaling.jl`.

| threads | `_v2` total (s) | new total (s) | `_v2` speedup | new speedup | new efficiency |
|---|---|---|---|---|---|
| 1 | 1322.8 | 2474.5 | 1.00× | 1.00× | 100% |
| 2 | 1039.9 | 1883.8 | 1.27× | 1.31× | 66% |
| 4 | 684.2 | 1162.6 | 1.93× | 2.13× | 53% |
| 8 | 378.0 | 735.8 | 3.50× | 3.36× | 42% |
| 16 | 220.8 | 412.5 | 5.99× | 6.00× | 37% |
| 32 | 126.0 | 251.1 | 10.50× | 9.85× | 31% |
| 64 | 81.4 | 119.5 | 16.26× | 20.70× | 32% |
| 96 | 66.0 | 97.7 | **20.04×** | **25.34×** | **26%** |

Absolute times roughly double (tracking *N*), but **the speedup improves,
20.0× → 25.3×**, and parallel efficiency with it (21% → 26%). The mechanism is
in the stage times: the DOF increase costs a consistent ~2.1× at 1–32 threads
but only ~1.6× at 64 and 96. At N = 1.62M the high-thread end was
overhead-bound; the larger problem has enough parallel work to amortize it. This
is a claim the paper can make in its own favour, and it is the opposite of what
a naive "bigger problem, worse scaling" reading would predict.

**Unresolved: the 32/64 point.** Efficiency is non-monotone there (31% at 32
threads on `worker7038`, 32% at 64 on `worker7105`), with a superlinear 2.16×
step in per-iteration time; both the precompute and GMRES stages show it, and
the kink is visible in panel (a). `rerun_2026_07/PAPER_CHANGES.md` documents this
exact failure mode once already (the published 32-thread row, 31% slow on
`worker7176`). Re-run those two points on different nodes before quoting them:

```
sbatch --reservation=rocky8 --array=5-6 slurm/run_multicube_threads_pinned_sec53.sbatch
```

---

## Table 3 — K sweep (96 threads), extended to K = 46

| K | `_v2` niter | new niter | `_v2` total (s) | new total (s) | `_v2` amort | new amort | new RSS (GB) | % of node |
|---|---|---|---|---|---|---|---|---|
| 1 | 14 | 12 | 100.9 | 125.2 | 1.00× | 1.00× | 42.8 | 3% |
| 4 | 13 | 11 | 168.3 | 232.8 | 2.67× | 2.32× | 118.7 | 8% |
| 10 | 13 | 10 | 344.1 | 436.4 | 3.35× | 3.28× | 277.0 | 18% |
| 13 | 12 | 10 | 415.4 | 536.5 | 3.69× | 3.50× | 354.7 | 24% |
| 19 | 12 | 10 | 582.4 | 776.2 | 3.90× | 3.53× | 512.7 | 34% |
| 25 | 12 | 10 | 743.7 | 991.4 | 4.00× | **3.66×** | 673.2 | 45% |
| 31 | 12 | 10 | 920.0 | 1236.1 | 4.07× | **3.66×** | 829.5 | 55% |
| 37 | — | 10 | — | 1491.5 | — | 3.65× | 987.0 | 66% |
| 40 | — | 10 | — | 1742.5 | — | 3.33× | 1063.4 | 71% |
| 46 | — | 10 | — | 2140.8 | — | **3.04×** | 1223.1 | **81%** |

Iteration counts fall on the lower contrast (14→12, 12→10) and are then **constant
at 10 for every K ≥ 10**, so nothing below is a conditioning effect.

### The non-monotonicity is real, and it is a memory effect

**Amortization peaks at 3.66× (K = 25–31), holds through K = 37, then declines to
3.04× at K = 46.** The paper's original claim of a decline beyond the optimum was
therefore correct in kind; the July fixes did not remove the turnover, they *moved
it* from K ≈ 25 to K ≈ 37. Both earlier statements need revising: the submitted
text put the peak too early, and the July note ("the crossover is pushed beyond
the range measured here") was an artifact of the sweep stopping at K = 31.

**The stated mechanism is wrong and must be replaced.** The paper attributes the
decline to block-Arnoldi orthogonalization being O(K²N) per iteration. That term
is negligible at these sizes: ~2 N K² m² is 0.5 TFLOP at K = 31 and 1.2 TFLOP at
K = 46, i.e. seconds of BLAS-3 against a ~300 s excess. The cost that actually
grows is the per-iteration block solve per RHS:

| K | 13 | 19 | 25 | 31 | 37 | 40 | 46 |
|---|---|---|---|---|---|---|---|
| block/K/niter (s) | 3.034 | 3.044 | 2.960 | 2.970 | 2.990 | **3.286** | **3.617** |
| peak RSS (% of node) | 24% | 34% | 45% | 55% | 66% | 71% | 81% |

Flat to within 3% while resident memory stays below ~66% of the node, then +22% at
71% and 81% occupancy. The onset tracks memory occupancy, not K.

Julia GC is *not* the mechanism either, and the records prove it: at K = 46 the
live heap is **27.3 GB against 1223 GB resident**. Almost the entire footprint is
native allocation inside FMM3D/TKM3D plus the sparse near-correction structures,
so the GC has nothing to thrash on. What remains is memory-system pressure —
bandwidth saturation and NUMA/first-touch locality degrading as the per-iteration
working set grows with K on a nearly-full node. Bandwidth and NUMA have **not**
been separated; claim it at that level and no finer.

Recommended framing: the batch size is bounded by memory, not by a compute
crossover — a statement that is both true here and more useful to a reader sizing
their own runs.

### K = 46 is the last reachable point

Peak RSS is linear at 26.2 GB per RHS (RSS ≈ 16.6 + 26.2 K GB; predicts 828.8 at
K = 31 vs 829.5 measured, 1221 at K = 46 vs 1223 measured). The next graphene
neighbour shell is K = 58 — the shells jump 46 → 58 with nothing between — needing
~1536 GB against a 1538 GB node. 1538 GB is the largest memory tier in `ccm`, so
K = 46 cannot be exceeded on hardware comparable to the rest of the curve.

Evaluation per RHS creeps only mildly over the whole sweep (6.5 → 7.2 s) and is
19% of the K = 46 wall time (333.1 s of 2140.8 s), continuing the trend §5.3's
Table 4 already showed — the cost is no longer solve-dominated.

Figure regenerated at `figs_sec53/fig66_multicube_scaling.pdf` with all ten K
points. Both annotations are data-derived and update themselves (25.3×, 3.7×).

**Two figure issues.** `xlims!(ax2, 0, 33)` was hard-coded and silently clipped
every point past K = 31, including the turnover — **fixed**, now derived from the
data. `ylims!(ax2, 0, 250)` is still hard-coded while the data spans 30.3–111.1;
the turnover is a 6.3 s rise, which in a 250 s window is 2.5% of the panel height,
so the headline result is nearly invisible. Deriving it (or clipping to ~0–120)
would fix that; left alone pending a decision, since some panels use a
deliberately fixed window.

---

## Evaluation cost: now tabulated, and it does not amortize

§5.2 previously timed Algorithm~(solve) only, excluding evaluation, with §5.3 timing the full
pipeline. That understated the cost a reader actually pays: evaluation is 16% of the K = 46
wall time and its per-source cost is essentially flat in K (6.1 s at K = 1 to 7.3 s at K = 46).
Folding it in:

| K | 1 | 4 | 10 | 19 | 25 | 31 | 37 | 40 | 46 |
|---|---|---|---|---|---|---|---|---|---|
| solve-only speedup | 1.00× | 2.32× | 3.28× | 3.53× | **3.66×** | 3.66× | 3.65× | 3.33× | 3.04× |
| end-to-end speedup | 1.00× | 2.16× | 2.89× | 3.08× | **3.18×** | 3.16× | 3.12× | 2.89× | 2.67× |

**The end-to-end batching gain is 3.2×, not 3.7×**, and the turnover is deeper. `tab:multirhs`
now carries an `evaluation (s)` column and the end-to-end speedup; `EVAL_IN_TABLE=0` restores
the solve-only columns.

### Why evaluation cannot be batched here (measured, not argued)

The obvious fix — reuse the multi-RHS machinery for the evaluation, i.e. call
`BI.evaluate_batch_potential`, which shares one corrected `pottrg` map across all K columns and
collapses every far target into a single `nd = K` FMM — **does not work for this problem, and
the reason is structural.**

`run_multicube.jl` grew a `MULTICUBE_EVAL` knob (`percolumn` / `batched` / `both`; `both` runs
each on the same interface and σ and checks they agree, which they do). Measured:

* the batched path is **4% slower** (10.76 s vs 0.16 + 10.13 s at K = 4, smoke tolerances);
* over the whole K = 46 K×K block, **0 of 703 million (target, source) point pairs are
  far-field**, in all 2116 (target set, source) combinations.

Every pair density ρ₁ⱼ lies within ~6.5 Å of the others with ~10 Å support, so at c_pad = 5
every orbital's support is inside every other's near region. `far_idx` is empty, the `nd = K`
FMM never executes, and the batched path degenerates to K per-column `PrecomputedVolumeField`
builds — paying extra bookkeeping and losing `ltkm3dc`'s own `kmax` tuning.

So the per-column evaluation loop is already near-optimal here, and the amortization claim is
**solve-only by structure, not by omission**. This is worth one or two sentences in §5.2: it
explains why §5.3's stage split inverted toward evaluation (Table 4: 50% solve / 40% eval), and
it tells a reader that batching buys throughput in the solve and nothing in the evaluation.
Batched evaluation would pay only for well-separated pairs, which a localized Wannier basis
inside one neighbour cutoff does not contain — §5.3's wider 10×10 lattice does, so the two
sections are not in conflict.

The default is `percolumn`, so the recorded timings and the tables are unchanged by this work;
`both` is what produced the comparison above.

---

## New: §5.2 and §5.3 now cross-validate

Because both sections solve the same system, their on-site Coulomb energy is
directly comparable for the first time — through two independent code paths
(this experiment's multi-RHS pipeline vs. `lattice_scale`'s campaign):

| | on-site *U* (eV) | *N* |
|---|---|---|
| §5.2, `exp66` pipeline | **7.1765** | 2,751,264 |
| §5.3, `lattice_scale` level 4 | **7.1889** | 2,752,560 |

Agreement to **0.17%** at matched DOF. The bare value is 17.537 eV, and
17.537/7.176 = 2.4437 reproduces the reported screening ratio exactly.

The residual 0.17% is not noise. §5.2 truncates its *target* set at
`support_rtol` (keeping `|ρ| ≥ 1e-4·max|ρ|`) while integrating the full source,
dropping the density tail asymmetrically and biasing low; §5.3's `onsite_U` uses
every point on both sides. §5.2's looser `rhs_tol` (1e-3 vs 1e-5) contributes as
well. The split between the two causes has **not** been measured — quoting this
as a cross-validation is fair; quoting 0.17% as a converged discretization
difference is not.

Worth noting for §5.3 independently: at ε = 2.4 the single-orbital *U* is flat to
4e-4 eV across all seven refinement levels (7.18873–7.18917), far tighter
convergence than the ε = 10 sweep showed.

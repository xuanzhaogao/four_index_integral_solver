# Changes §5.2 needs after moving it onto the §5.3 system (2026-09-07)

§5.2 (performance of the batched interface solve — Tables 2, 3 and the scaling
figure) was benchmarked on a substrate the physics section does not use: two
90×90×90 cubes with a 90×90×9 slab at ε = 10, a generic demo value. It now runs
the **same system §5.3 reports** — the ×3 converged Si|SiO₂ heterojunction
(270³, ε 11.9 | 3.9) with the graphene slab at its cRPA ε = 2.4.

> **Provenance.** This note was revised as the measurements improved, and earlier revisions
> stated things later runs refuted. Current as of 2026-09-09: the accuracy set is the lattice
> campaign's own (edge_refine_level 3, all tolerances 1e-5, **max_order 64**), giving
> N ~ 1.35e6, and the timing tables come from job 7006943 (all K sequentially on ONE node).
> Sections marked WITHDRAWN or superseded describe claims that did not survive; they are kept
> rather than deleted so the same conclusions are not re-derived.

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
| surface unknowns *N* | 1,617,264 | **1,350,288** |
| screening ratio `V_vac/V` | 8.3948 | **2.4437** |

**The abstract's "1.6 million surface unknowns" must become 1.35 million.**

That figure moved twice. At the submitted accuracy (`rhs_tol` 1e-3, edge_level 4) the sec53
geometry gives N = 2,751,264, and an earlier revision of this note said 2.75 million. But §5.2
was then running at an accuracy that matched neither the submitted §6.3 nor §5.3: a *finer*
edge refinement (l_ec 0.568 vs 1.136 Å) with *looser* source tolerances (1e-3 vs 1e-5). Setting
it to the lattice campaign's own accuracy gives **N = 1,350,288**, which matches `lec_single`
level 3 (1,350,288) exactly — independent confirmation that the two sections now discretize the
same interface.

*N* varies by 0.05% across the sweep (1,349,856 at K = 1 to 1,350,504 at K = 46): at
`rhs_tol` 1e-5 the interface refinement follows the source envelope, where at 1e-3 it was
pinned. Immaterial for timings, but the tables should not claim a single N for all K.

---

## Table 2 — strong scaling (K = 1)

`prepare` measured directly; `block BIE solve` reconstructed as
`t_RHS + t_1iter · niter_ref` (niter_ref = 11); `evaluation` measured in the same run at the
fixed lattice target set (N_p = 7,168,390).

| threads | prepare (s) | block solve (s) | evaluation (s) | total (s) | speedup | efficiency |
|---|---|---|---|---|---|---|
| 1 | 133.6 | 759 | 127 | 1020 | 1.0× | 100% |
| 2 | 74.9 | 591 | 75 | 741 | 1.4× | 69% |
| 4 | 40.2 | 381 | 47 | 469 | 2.2× | 54% |
| 8 | 23.8 | 254 | 33 | 311 | 3.3× | 41% |
| 16 | 14.4 | 146 | 22 | 183 | 5.6× | 35% |
| 32 | 10.5 | 84 | 18 | 112 | 9.1× | 28% |
| 64 | 8.7 | 50 | 15 | 74 | 13.9× | 22% |
| 96 | 8.9 | 39 | 14 | **62** | **16.4×** | 17% |

End to end, 1020 s on one thread to 62 s on 96, **16.4×** at 17% parallel efficiency.

The evaluation column scales monotonically to 96 threads (127 → 14 s) only because of the
`pottrg` threading fix (BoundaryIntegral 471ce6e) and dropping `OPENBLAS_NUM_THREADS` from N to
1. Before those, the evaluation at 96 threads took **99–155 s against 21.6 s at 64** — 96
threads was the worst configuration on the node. See the threading section below.

---

## Table 3 — K sweep, and two claims withdrawn

Measured sequentially on ONE node (job 7006943) at the lattice campaign's own accuracy
(edge_refine_level 3, all tolerances 1e-5, max_order 64, N ~ 1.35e6), evaluated at the fixed
target set of the lattice model (N_p = 7,168,390):

| K | 1 | 4 | 10 | 19 | 31 | 37 | 40 | 46 |
|---|---|---|---|---|---|---|---|---|
| prepare (s) | 10.0 | 11.0 | 12.3 | 10.8 | 12.1 | 12.0 | 11.1 | 11.1 |
| block solve (s) | 49 | 88 | 175 | 317 | 508 | 552 | 588 | 675 |
| evaluation (s) | 16 | 40 | 98 | 180 | 321 | 371 | 395 | 472 |
| total (s) | 76 | 140 | 285 | 507 | 841 | 934 | 994 | 1158 |
| peak RAM (GB) | 24 | 62 | 143 | 263 | 428 | 512 | 550 | 634 |
| end-to-end speedup | -- | 2.2x | 2.6x | 2.8x | 2.8x | 3.0x | 3.0x | **3.0x** |
| solve-only speedup | 1.00x | 2.39x | 3.17x | 3.44x | 3.53x | 3.89x | 3.96x | **3.97x** |

### WITHDRAWN: the amortization turnover

Earlier revisions of this note reported a peak near K = 25--37 followed by a decline (4.07x
falling to 3.04x), and attributed it to memory pressure. **Both claims are withdrawn.**

The per-K points had each been measured on a separate node, one sample each. Comparing two such
runs at identical settings showed node-to-node differences of up to 26% in the block solve
(K = 4: 87 s vs 110 s; K = 46: 785 s vs 662 s) -- larger than the 16% "decline" being read as
signal, and in one run K = 46 came out CHEAPER than K = 40, which is impossible for strictly
more work. Re-measuring all ten points sequentially on one node removes that variation and the
curve is monotone: solve-only amortization rises to 3.97x at K = 46, end-to-end plateaus at
3.0x. There is no turnover in the measured range.

The memory explanation had already failed independently: the same apparent turnover appeared at
42% node occupancy as at 81%, so it tracked K rather than memory. And the O(K^2 N)
orthogonalization the submitted text cites is ~1.2 TFLOP at K = 46, seconds against a solve of
hundreds -- it cannot produce a turnover either.

For the paper: the submitted claim of a decline beyond an optimum is not supported by data at
this accuracy, and neither is the July note's recast ("the crossover is pushed beyond the range
measured here"). The defensible statement is that batching amortizes the fixed per-batch cost
monotonically over the measured range, saturating near 3x end-to-end, and that the practical
ceiling is the neighbour-shell structure (K = 46 is the last shell below K = 58), not memory --
634 GB at K = 46 is 42% of the node.

**Any single-node timing quoted from a per-node array run carries ~20% uncertainty.** Structure
smaller than that is not resolvable without repeated or sequential measurement.

### Evaluation amortizes too

The noisy data hid this. Per-source evaluation cost falls 16.35 s at K = 1 to 9.5--10.3 s and
then stays flat -- a ~1.65x gain saturating by K = 10. That is why the solve-only speedup
(3.97x) exceeds the end-to-end figure (3.0x): the evaluation stops improving earlier, not
because its cost is independent of K.

---

### Extension to K = 79, and a possible decline

`nrange` had to be fixed first: it was pinned at 3, which silently truncates the neighbour
search past cutoff ~7.0 A (66 sites found where 67 exist at 7.4, 74 of 79 at 8.0), so extending
the sweep would have run a wrong cluster with no error raised (BoundaryIntegral-side commit
2045f00 in this repo's tree).

Four new shells, one node each (job 7008840), against the sequential K = 46 point:

| K | 46 (seq) | 58 | 67 | 73 | 79 |
|---|---|---|---|---|---|
| block solve (s) | 675.0 | 846.1 | 1022.2 | 1130.1 | 1258.4 |
| eval (s) | 472.3 | 593.3 | 683.4 | 741.6 | 829.5 |
| total (s) | 1435.7 | 1849.2 | 2187.1 | 2343.8 | 2631.9 |
| peak RSS (GB) | 634 | 783 | 904 | 982 | 1064 |
| % of node | 42% | 52% | 60% | 65% | 71% |
| solve/RHS (s) | 14.9 | 14.8 | 15.5 | 15.7 | 16.1 |
| solve amortization | 3.97x | **3.99x** | 3.82x | 3.78x | **3.67x** |
| e2e amortization | 3.00x | **3.02x** | 2.94x | 2.93x | **2.84x** |

Solve-only amortization peaks at K = 58 and declines monotonically to 3.67x at K = 79, 8% off
peak. It is entirely in the solve: per-RHS-per-iteration cost rises 1.621 -> 1.770 s (+9%) at
constant N_iter = 9, while per-source evaluation stays flat at 10.2--10.5 s.

**NOT ESTABLISHED.** The four points are on four different nodes and the 8% decline is smaller
than the 20--26% node-to-node spread measured earlier -- the same confound that manufactured a
phantom turnover at K ~ 37. Four monotone points weigh more than one outlier, and the onset
tracking occupancy (52% -> 71%) is at least consistent with memory pressure, unlike the phantom
which appeared at 42% and 81% indifferently. Settling it needs the five points measured
sequentially on one node (~2.5 h).

**The memory model is confirmed to 0.7%** out to K = 79: RSS ~ 10.6 + 13.4 K GB predicts
788/908/989/1069 against 783/904/982/1064 measured. K = 79 is the practical ceiling -- the next
shell, K = 103, projects to 1391 GB (93% of the node). What occupies 13.4 GB per RHS
(~1240 doubles per unknown) is still unaccounted for.

## Evaluation cost: now tabulated (superseded in part -- see above)

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

### Why evaluation cannot be batched OVER A SINGLE ORBITAL'S SUPPORT (superseded)

**Superseded.** The conclusion below is correct for the target set it used -- one
orbital's support, where there is no far field and no shared cost worth dividing --
but not for the evaluation the pipeline actually performs. Measured at the lattice
model's own target set (N_p = 7,168,390), per-source evaluation cost falls 16.35 s to
~9.8 s, a 1.65x amortization saturating by K = 10. The original text follows.


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

## §5.3: max_order = 8 silently under-resolved the near field

`max_order` caps the near-pair upsampling order: p_up = ceil(-log(up_tol) / (2 log rho_min)),
clamped to [n_quad, max_order]. When the clamp binds the correction is computed at a lower
order than its own tolerance demands, there is **no adaptive fallback for non-touching pairs**
(numerical_results/PLAN.md), and **no warning is emitted**.

Both lattice campaign generators hardcoded `max_order = 8` with no comment. Measured on a real
solved conv_l3_eps2.4 interface (37,520 panels, N = 1,350,720, up_tol 1e-5, n_quad 6):

| max_order | upsample pairs | p_up median | at the cap |
|---|---|---|---|
| **8** | 12,208 | 8 | **11,296 (92.5%)** |
| 16 | 12,208 | 16 | 9,136 (74.8%) |
| 32 | 12,208 | 32 | 8,380 (68.6%) |
| 64 | 12,208 | 64 | 6,700 (54.9%) |

So **every published §5.3 campaign ran with 92.5% of its non-touching near corrections clamped
below the order its stated 1e-5 tolerance required.** Corrections apply out to
rho* = 2.61; p_up <= 8 holds only for rho_min >= 2.05, so the whole band 1 < rho_min < 2.05 --
the CLOSEST pairs, which dominate the near error -- was under-resolved.

**The l_ec self-convergence study cannot detect this.** It varies edge refinement; a systematic
near-field under-resolution is present at every level, so the sweep converges cleanly to a
slightly wrong answer. U flat to 4e-4 eV across seven levels is consistent with both
"converged and correct" and "converged and biased".

All campaigns and both generators now use 64. Effect on the tensor, conv_l3_eps2.4 at
k_target = 31 (job 7006478) against the max_order = 8 reference: max relative difference
4.7e-3 (474x gmres_rtol), rms 1.3e-5, on-site U mean 7.16798 -> 7.16281 eV, max |dU| 0.0348 eV.
Concentrated in few entries, as a near-field quadrature fix should be.

**Si/SiO2 contrast survives:** 0.15612 -> 0.15556 eV (0.36%). The position-dependent screening
claim is unaffected.

A residual gap remains: **55% of upsample pairs are still at the cap at max_order = 64.** Those
have rho_min < 1.094 -- nearly coincident but not edge-adjacent, so the touching test routes
them to upsampling where they would need p_up in the hundreds. Widening the touching criterion
to send them to the adaptive quadtree is the real fix; raising max_order further is not.

---

## §5.3 Table 4 — new phase timings (85 batches, 10 nodes, max_order 64)

| phase | time | share |
|---|---|---|
| prepare | 10 m 43 s | 5.6% |
| solve | 1 h 06 m 01 s | 34.6% |
| **eval** | **1 h 50 m 25 s** | **57.9%** |
| consolidate | 3 m 05 s | 1.6% |
| assemble | 37 s | 0.3% |
| **total** | **3 h 10 m 51 s** | **31.8 node-hours** |

Evaluation now dominates at 58% (was ~40%), and the campaign costs 31.8 node-hours against the
July rerun's 14. The batching and triangle work did reduce the solve side (85 batches instead
of 198, ~2x less evaluation work from symmetry), but `max_order = 64` more than absorbed it:
the near corrections now actually resolve. This is accuracy bought with time, and §5.2 should
say so rather than present it as a speedup.

It also relocates where optimisation is worth doing: the evaluation, not the solve.

Symmetry: 9.34e-6 over independently computed pairs, 3,399,588 of 6,880,129 entries mirrored.

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

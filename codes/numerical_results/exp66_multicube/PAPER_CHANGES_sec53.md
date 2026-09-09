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

---

## §5.3 rerun: proper K, and the accuracy findings applied (2026-09-09)

Three of the findings above were only half-applied to §5.3. The production tensor was rerun
(job 7006478, `lattice_conv_l3_eps2.4_kb31`) with `max_order = 64`, K-balanced batching, the
symmetry-restricted evaluation and `OPENBLAS_NUM_THREADS = 1`. Not yet applied:

1. **K was 31, not the measured optimum.** `k_target = 31` was chosen when the K sweep still
   showed a turnover near 37 -- the phantom that was later withdrawn. The sequential sweep puts
   the cost floor flat over K = 37-58.
2. **The l_ec self-convergence sweep still ran at `max_order = 8`.** `scripts/lec_conv_single.jl`
   hardcoded it; only the two full-campaign generators had been fixed. So §5.3's convergence
   figure (`figs/lec_single_conv_eps2.4.tsv`, 7 levels) rests on the under-resolved near field,
   and, as noted above, the sweep is structurally unable to see that.
3. **Only `run_all.sbatch` had the OpenBLAS fix.** The other ten jobscripts still exported
   `OPENBLAS_NUM_THREADS = $SLURM_CPUS_PER_TASK`, including the two that run these campaigns.

All three are now fixed in the repo. `k_target` is also emitted by both generators (it had been
hand-added to the TOMLs, so any regeneration silently dropped it back to per-anchor batching).

### Why K = 46

Re-partitioning the campaign's actual 2623-pair list at each candidate:

| k_target | batches | K spread | envelope median | rounds on 10 nodes |
|---|---|---|---|---|
| 31 | 85 | 30-31 | 7.9 A | 9 |
| 37 | 71 | 36-37 | 8.1 A | 8 |
| **46** | **57** | **46-47** | **8.9 A** | **6** |
| 58 | 45 | 58-59 | 9.3 A | 5 |

46 sits on the cost floor (25.2 s per pair against 27.1 s at K = 31), quantises best against 10
nodes, inflates the source envelope only 12%, and stays clear of the possible K > 58 decline.
Peak RSS is 634 GB, 42% of the node, against 428 GB at K = 31.

**Expect only ~6% off the wall clock, not the 7.5% the per-pair cost implies.** Scaling the
measured kb31 phases: solve 66 -> ~59 min (K amortizes), evaluation 110 -> ~107 min (per-source
cost is flat past K = 10, so its total work is fixed by the pair count, not by K), prepare and
consolidate unchanged. ~3 h 00 m against 3 h 11 m. The tensor should be unchanged -- the two
campaigns differ in `name`, `root` and `k_target` and nothing else -- which makes the run a
K-independence check on V as well as a timing measurement.

### The resume trap

`solve_batch` skips any batch whose result file already exists, so **a rerun that changes the
solver must change the campaign name**, or it silently re-emits the old numbers. Roots holding
`max_order = 8` results, none of which may be resumed:

    lattice_conv_l3_eps2.4          lattice_10x10_het3x_eps2.4      lec_single_l{1..7}_eps2.4

Hence `lattice_conv_l3_eps2.4_k46` and `RERUN_TAG=_eps2.4_mo64` below, both fresh.

### Result: the l_ec sweep at max_order = 64 (job 7009683, 34 min, 1 node)

| level | l_ec (A) | dof | U (mo8) | U (mo64) | dU (eV) | t8 | t64 |
|---|---|---|---|---|---|---|---|
| 1 | 4.5450 | 308,952 | 7.1890899 | 7.1894404 | +3.50e-04 | 25 | 31 |
| 2 | 2.2725 | 653,040 | 7.1891741 | 7.1893268 | +1.53e-04 | 26 | 38 |
| 3 | 1.1362 | 1,350,288 | 7.1890041 | 7.1890427 | +3.86e-05 | 61 | 64 |
| 4 | 0.5681 | 2,752,560 | 7.1888593 | 7.1888963 | +3.70e-05 | 104 | 122 |
| 5 | 0.2841 | 5,564,880 | 7.1887807 | 7.1888119 | +3.12e-05 | 209 | 232 |
| 6 | 0.1420 | 11,197,296 | 7.1887508 | 7.1887704 | +1.96e-05 | 418 | 461 |
| 7 | 0.0710 | 22,469,904 | 7.1887323 | 7.1887446 | +1.23e-05 | 875 | 976 |

**Convergence in l_ec survives**, and U is unchanged for reporting purposes: at level 4, the
value Sec. 5.3 quotes, the near-field fix moves U by +3.7e-05 eV (5e-06 relative). The bias
decays monotonically with refinement, 3.5e-04 at level 1 to 1.2e-05 at level 7 -- refinement
raises rho_min and unclamps the corrections, so edge refinement was already absorbing most of
the max_order = 8 error. Resolving the near field costs 11-24% more time per level.

So **the l_ec figure's numbers stand.** That is worth stating plainly: the max_order finding
invalidated the near-field accuracy those campaigns *claimed*, and for this observable it turns
out not to have moved the answer.

### But the single-orbital sweep is not a proxy for the campaign

For the SAME orbital (idx 104, the one the sweep picks) at the SAME level 3:

| | U at mo8 | U at mo64 | shift |
|---|---|---|---|
| campaign `conv_l3` (orbital in a 31-pair batch) | 7.1547889 | 7.1501961 | **-4.59e-03 eV** |
| l_ec sweep (that orbital alone, neighbor_cutoff = 0) | 7.1890041 | 7.1890427 | **+3.86e-05 eV** |

The values differ by 0.034 eV (0.48%) and the near-field sensitivity by 124x, opposite in sign.
Orbital 104 is typical of the campaign, not an outlier -- 82nd of 198 by |shift|, against a mean
of -5.18e-03 and a median of -3.68e-03 eV.

**This is a second and independent reason the l_ec study cannot certify the campaign's
accuracy.** The first (above) is that a bias present at every refinement level is invisible to a
refinement sweep. This one is that the sweep does not even carry the same bias: it is a
different configuration whose near-field error is ~100x smaller. Whatever it establishes, it
establishes for one isolated orbital.

Mechanism NOT established. Two candidates were eliminated by measurement:

* *interface size* -- the campaign's per-batch interface is 1,350,504-1,352,556 unknowns against
  the sweep's 1,350,288, identical to 0.1%, so the campaign is not coarser overall;
* *support_rtol truncation* -- applied in `solve_batch_core` (tasks.jl:188) to the batch support
  at solve time, and the sweep's `onsite_U` contracts over that same truncated support, so both
  paths truncate identically.

Still open: refinement *placement* (equal panel count says nothing about where the panels go --
the campaign refines one shared interface toward the joint envelope of 31 pairs, the sweep
toward one orbital's support), and the far larger set of near panel-target corrections in the
campaign's `pottrg` map, which spans all 198 orbitals' supports rather than one.

### Queued

    # production tensor at K = 46, ~3 h on 10 nodes (~30 node-hours)
    sbatch --nodes=10 jobscripts/run_all.sbatch campaigns/lattice_conv_l3_eps2.4_k46.toml

    # l_ec sweep at max_order = 64, 7 levels, ~1 h on 1 node
    sbatch --export=ALL,RERUN_TAG=_eps2.4_mo64 jobscripts/run_lec_conv_eps2.4.sbatch 1 2 3 4 5 6 7

The l_ec rerun is the one with scientific content: it says whether U's flatness to 4e-4 eV
across seven levels survives a resolved near field, and it re-anchors the §5.2/§5.3
cross-validation, whose §5.3 side (7.1889 eV at level 4) is a `max_order = 8` number.

**Done** (job 7009683, results above): flatness survives, U moves 3.7e-05 eV at level 4, and the
cross-validation's §5.3 side is unchanged to five figures. The K = 46 campaign (7009682) is
still running.

---

## The K = 46 rerun (job 7009682): k_target is not a performance knob

Finished 2026-09-09, 3 h 05 m on 10 nodes against kb31's 3 h 11 m -- **2.8% faster, not the 6%
projected.** prepare rose 10.7 -> 14.6 min (RCB plus wider envelopes), eating half the gain that
solve (66.0 -> 64.5) and eval (110.4 -> 102.3) delivered.

### The tensor is NOT batching-invariant

kb31 and k46 differ in `name`, `root` and `k_target` alone, so V was expected to reproduce to
solver tolerance. It does not:

| | |
|---|---|
| max relative difference | **1.555e-03** (155x gmres_rtol) |
| rms relative difference | 9.058e-06 |
| on-site U mean | 7.16281 -> 7.15933 eV |
| max abs on-site shift | 1.138e-02 eV |

The difference is **entirely in the on-site entries** -- all 200 worst entries are onsite-onsite,
1% of the matrix; the median entry agrees to 2.2e-08 and 99% to 3.3e-05. That is precisely the
quantity Sec. 5.3 reports.

The cause is structural, not a bug: each batch refines ONE shared interface on the envelope of
its K pair sources. Same panel count (1,350,504-1,352,556, within 0.1% of the single-orbital
1,350,288) but a different placement -- and the on-site integral, the most peaked and
near-field-sensitive entry, is what notices. On-site U for orbital 104, all at level 3 and
max_order 64:

| batching | U (eV) |
|---|---|
| K -> 1 (isolated, `neighbor_cutoff = 0` -- the l_ec sweep) | 7.18904 |
| K = 31 | 7.15020 |
| K = 46 | 7.14974 |

Monotone and saturating: -3.9e-02 eV from K=1 to 31, then -4.6e-04 to K=46. **The campaign sits
0.55% below the K -> 1 limit**, and that is where the 0.48% discrepancy reported earlier between
the campaign and the single-orbital sweep comes from. The mechanism I first proposed (refinement
placement) and then wrongly dismissed on the grounds of equal panel count is the surviving one:
equal N says nothing about where the panels go.

### The error hierarchy is inverted from what the paper claims

For on-site U at level 3, largest first:

| source | size |
|---|---|
| batching, K = 1 -> 31 | **3.9e-02 eV** (0.55%) |
| max_order 8 -> 64 | 5.2e-03 eV (mean) |
| batching, K = 31 -> 46 | 3.5e-03 eV (mean) |
| **l_ec self-convergence at level 3** | **3.4e-04 eV** (4.7e-05 rel) |

Sec. 5.3 quotes the l_ec sweep as its accuracy evidence, and it is the SMALLEST term by two
orders of magnitude. This is the third and strongest reason that sweep cannot certify the
campaign: it runs at K = 1 by construction, so it is blind to the dominant error by design.

### The physics survives

| | Si (eps 11.9) | SiO2 (eps 3.9) | contrast |
|---|---|---|---|
| max_order 8, K ~ 13 | 7.08992 | 7.24604 | 0.15612 eV |
| max_order 64, K = 31 | 7.08502 | 7.24059 | 0.15556 |
| max_order 64, K = 46 | 7.08180 | 7.23686 | 0.15506 |

Both sides move down together, so the contrast -- a small difference of two large numbers -- takes
a ~0.3% relative hit per perturbation while the absolute values move 0.05%. Position-dependent
screening is unaffected qualitatively and stable to under 1% across every perturbation tried.
**What must change is the error bar, not the claim.**

---

## The eval phase is 8x slower in the campaign than standalone (job 7010961)

The K sweep on ONE batch (kb31 batch 1, full 7.33M target set, current BI, single process at
-t 96), against the Aug 30 production rate of 8.23 s per source:

| K | 4 | 8 | 13 | 17 | 24 | 31 |
|---|---|---|---|---|---|---|
| s per source | 9.90 | 7.92 | 8.06 | 7.33 | 7.58 | 7.57 |
| vs Aug 30 | 1.20x | 0.96x | 0.98x | 0.89x | 0.92x | 0.92x |

Flat. **K does not explain the eval slowdown**, and neither does OpenBLAS (7.73 s per source at
1 thread vs 7.57). The same batch, same sigma, same targets, same K, same BI:

| | s per source |
|---|---|
| kb31 production | **61.4** |
| this A/B, standalone | **7.57** |

**8.1x.** So the eval cost is set by how the campaign runs evaluation, not by the mathematics.
The one known environment difference: `run_phase` launches workers with `-t $JULIA_GLUE_THREADS`
and both campaign jobscripts hardcode **8**, so every worker evaluates with 8 Julia threads on a
96-core node, while this A/B and the Sec. 5.2 benchmark -- both fast -- use 96.

NOT yet established, and one observation resists it: Aug 30 also ran 8-thread workers and was
fast at K <= 17. Only a thread ceiling that binds once K exceeds it fits both.
`jobscripts/ab_eval_threads.sbatch` tests 8/16/32/96 at K = 17 and 31.

If it holds, ~1.6 h of the k46 run's 1.9 h evaluation was waste, the campaign drops from ~3 h to
~1.3 h, and **"evaluation dominates at 58%, so optimise there" is an artefact of the thread
setting rather than a property of the algorithm** -- that sentence must come out. It also means
k_target was chosen against a badly distorted cost structure and should be re-derived afterwards,
which now matters for accuracy and not only for speed.

---

## Resolved: the eval cost was the symmetry restriction's contraction (BI b41ee80)

Every hypothesis in the section above was wrong, and each was killed by measurement:

| candidate | verdict |
|---|---|
| K | K sweep flat: 9.90, 7.92, 8.06, 7.33, 7.58, 7.57 s per source at K = 4..31 vs Aug 30's 8.23 |
| OpenBLAS threads | 1 vs 96 identical (7.73 / 7.50); and 8x96 vs 8x1 identical (231.9 / 234.7 s) |
| worker Julia threads | -t 8 and -t 96 identical -- all four cells 232-240 s |
| max_order | reaches `solve_batch_core` only |
| I/O from ceph | `t_setup` is 1% of eval (7-10 s per batch) |
| load imbalance | both parallel phases run at 96% utilisation |
| BI version | 2aac4ac..e95f8af touches four files; the only hot-path one is gated behind an empty neighbour list |

The A/B had been timing `evaluate_batch_potential`, but production records `t_phi` for all of
`eval_batch_core`. For kb31 batch 1 (K = 31, full target set): 235 s standalone against a
production `t_phi` of 1896 s. The 1661 s between them is what the symmetry restriction added.

Owning a subset of rows means owning a subset of the shared target set, so the contraction had
to remap every global target index into that subset -- a `Dict` of up to 7.3M entries built per
batch, probed once per target per row, inside a serial loop. The removed function's docstring
had asserted "the saving is NOT in the contractions -- those are cheap dot products." It was
never measured.

The trade it was making: **-2.4 node-hours of potential evaluation, +14 node-hours of
contraction.** Across the phase the potential evaluation is 3.2 of kb31's 18.1 node-hours.

Every row is now evaluated (BI b41ee80). Consequences:

* **eval work becomes `n_pairs x n_targets`, independent of the batching** -- 2623 x 7.33M =
  19.2 G source-points, ~5.5 node-hours, plus ~1.0 for the contraction at Aug 30's direct-index
  rate. Against kb31's 18.1 that is **2.8x**, and a campaign drops ~3 h 05 m to **~2 h**.
* `max_rel_asym` regains its meaning. Nothing is mirrored, so every off-diagonal entry is
  computed twice by different batches from different interfaces -- an end-to-end check rather
  than the tautology it had become. The pipeline test's 1e-8 symmetry assertion was measuring
  the mirror, not the solver.
* One accuracy mechanism is removed: batches no longer refine their interfaces
  (`_refine_interface_for_targets`) against different target sets, so an entry's quadrature no
  longer depends on which batch owned it. The solve-side envelope dependence remains.

### What this does to the k_target choice

With eval batching-independent, only prepare and solve still respond to K, and they barely do:

| | batches | prepare | solve |
|---|---|---|---|
| Aug 30, K ~ 13 | 198 | 9.1 min | 68.9 min |
| kb31, K = 31 | 85 | 10.7 | 66.0 |
| k46, K = 46 | 57 | 14.6 | 64.5 |

Solve improves 6% from K = 13 to K = 46 and prepare gets 5.5 min worse, so the whole K-balanced
batching idea is worth about **4 minutes in 185** -- while costing 0.55% on the on-site U, the
quantity Sec. 5.3 reports. That is a bad trade, and it now points the same way accuracy does:
**k_target should go back down.** K = 31 already beats K = 46 on total time once eval is fixed.

---

## Confirmed: full-matrix evaluation, benchmarked (job 7012551)

Eval + assemble re-run on the finished k46 solve, BI b41ee80:

| | symmetry-restricted | full matrix | |
|---|---|---|---|
| eval work | 16.39 node-hours | **7.11** | 2.31x |
| eval wall, 10 nodes | 102.2 min | **46.2 min** | **2.21x** |
| utilisation | 96% | 92% | |
| per-batch `t_phi` | 1902 -> 143 s | 390 / 414 / 556 (min/median/max) | |
| **spread** | **13.3x** | **1.43x** | |

**The ramp is gone**, which is the diagnosis confirming itself: it was the shrinking row count,
not anything physical. The residual 1.43x matches the node-to-node variation measured earlier.

Campaign total: prepare 14.6 + solve 64.5 + consolidate 3.5 + eval 46.2 + assemble 0.7 =
**2 h 10 m, 21.6 node-hours** against 3 h 05 m and 30.9.

Prediction was 39 min / 22,279 s, so I under-predicted by 15-18%. The potential-evaluation term
was fitted from six measured points and is solid; the contraction term was back-derived as a
residual of a residual from the Aug 30 run, and it is the part that missed -- 6,016 s actual
against 2,705 s predicted, 2.2x. That was the flagged weak link, and the flag was right about
which term would fail while understating by how much.

### Validation

    entries filled by symmetry: 0 of 6,880,129
    max rel asymmetry:          3.032e-05        (over every off-diagonal entry)

Against the tensor the symmetry-restricted run produced from the same sigmas and interfaces:

| | |
|---|---|
| max relative difference | 3.029e-05 |
| rms relative difference | 6.839e-07 |
| on-site U, max abs change | 2.99e-10 eV |

Three things fall out, and each is a check rather than a coincidence:

* **The max difference equals the asymmetry** (3.029e-05 vs 3.032e-05). It must: the entire
  difference is that the previously mirrored half is now computed independently, so its size is
  exactly the V[i,j] / V[j,i] disagreement. Not zero (the new code really is not mirroring), not
  larger (nothing else changed).
* **On-site U is unchanged to 3e-10 eV.** Diagonal entries are their own transpose and were
  never mirrored, so they had better not move. They do not.
* **`max_rel_asym` is now a real end-to-end error indicator.** The 1.71e-05 the previous run
  reported covered only entries two batches happened to share; 3.03e-05 covers all 6.88M
  off-diagonal entries, each computed twice from different interfaces and sigmas. At ~3x the
  solver tolerance this is the number Sec. 5.3 should quote as its pipeline-level accuracy --
  it is the only figure in this whole investigation that is a genuine independent check of the
  full chain.

### Still open

`k_target` should come down: with eval now batching-independent, K buys ~4 minutes in 185 and
costs 0.55% on the on-site U (see the section above). The k46 tensor stands as the current
production result, but it is not the most accurate one obtainable.

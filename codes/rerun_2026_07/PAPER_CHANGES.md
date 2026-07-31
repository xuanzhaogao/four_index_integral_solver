# Changes the paper needs after the July 2026 rerun

Every number below is measured, not projected. All seven rerun jobs completed
2026-07-31; no results remain outstanding.
Line numbers refer to `~/Articles/four_indices_bie/main.tex` as of 2026-07-31.

Two classes of change are mixed here and should be kept separate when editing:

* **(F)** consequences of the two LHS fixes — these are what the rerun was for.
* **(D)** pre-existing defects found while auditing, unrelated to either fix.
  These would have needed correcting even if the rerun had changed nothing.

---

## §6.1 — Convergence and dielectric contrast

**(F) `:1205` — GMRES counts.** "GMRES counts of $16$--$24$ that stay flat as
$N$ grows" → **9--13**. The flatness claim is unaffected; only the numerals move.

**(F) `:1212` — contrast growth.** "The GMRES count roughly doubles over the same
range, from $16$--$18$ to $35$--$43$" → grows by about **70%**, from **10--11 to
17--19**. The mesh-independence clause in the same sentence still holds.

**(F) Table 1 — `N_iter` column** (`E_r` column is unchanged):

| $\epsilon_2$ | published | new |
|---|---|---|
| 6 | 16–18 | 10–11 |
| 20 | 24–28 | 12–13 |
| 60 | 29–36 | 15–17 |
| 200 | 35–43 | 17–19 |

**(F) Fig 8** — panel (b) is unchanged (errors reproduce across all 15 points);
panel (c) drops as above. Regenerated at
`codes/numerical_results/figs_v2/fig61_fig1_combined.pdf`.

**(F) `:418`** — "the GMRES iteration count is expected to be essentially
independent of the resolution" still holds and is now better supported.

---

## §6.3 — Performance of the batched interface solve

**(F) Table 2 — every absolute time falls ~2.4×** via $N_{\mathrm{iter,ref}}$
$36 \to 14$. Speedup and parallel efficiency are ratios and barely move
(20.3× → 20.0× at 96 threads; 21% efficiency either way).

**(D) Table 2, the 32-thread row was contaminated.** Per-iteration cost matches
the rerun to within 1% at every other thread count and is **31% slow** at 32
(`worker7176`). `seff` confirms the task got its full exclusive node, so this is
node-level slowness, not misallocation. It is what produces the efficiency dip
36% → 25% → 25% in the submitted table; the rerun declines monotonically
37% → 33% → 25%. **The published row is wrong independently of either fix.**

**(F) Table 3** — iteration counts $36$–$29 \to 14$–$12$; times fall ~2.4×;
peak memory essentially unchanged (537.6 → 523.9 GB at $K=31$).

**(F) `:1266` — the non-monotonicity claim must be withdrawn.** The paper says
the amortization advantage "is not monotonic ... and beyond $K\approx25$ it
declines, to $3.2\times$ at $K=31$". It no longer does:

| $K$ | 1 | 4 | 10 | 13 | 19 | 25 | 31 |
|---|---|---|---|---|---|---|---|
| published | — | 1.95 | 3.19 | 3.33 | 3.56 | **3.62** | 3.16 |
| rerun | — | 2.68 | 3.36 | 3.69 | 3.90 | 4.01 | **4.07** |

The stated *mechanism* survives — orthogonalization is $\mathcal{O}(K^2N)$ *per
iteration*, so 2.4× fewer iterations defers the crossover past $K=31$. Recast as
"the crossover is pushed beyond the range measured here" rather than deleting the
reasoning.

---

## §6.4 — Lattice-scale four-index tensor

**(F) Table 4 — the campaign is 2.8× cheaper and the stage split inverts:**

| stage | published | new |
|---|---|---|
| prepare | 19 s *(see D below)* | 5 m 03 s |
| solve | 3 h 14 m | **41 m 20 s** |
| consolidate | 3 m 32 s | 3 m 19 s |
| eval | 32 m 17 s | 32 m 52 s |
| assemble | 39 s | 38 s |
| total | 3 h 51 m | **1 h 23 m** |
| node-hours | 39 | **14** |

**(F) `:1320` — the stage-percentage sentence.** "the solve and evaluation
stages, which account for $84\%$ and $14\%$" → **50% and 40%**. The conclusion
drawn from it — that cost is FMM-dominated block solve — no longer holds;
evaluation is now nearly co-equal. "about $0.2$ node-hours per orbital batch" →
about **0.07**.

**(D) Table 4's `prepare` row was measured on a resumed run.** The published
campaign's `manifest.tsv` is stamped 12:04:47 while its own prepare banner reads
12:10:26 — prepare had already been done by an earlier session, so the 19 s was a
no-op verifying existing files. A genuine cold prepare is ~5 min: both fresh
rerun campaigns took ~5 m, as did both published `conv` campaigns. Note this
correction runs *against* the paper's interest.

**(D) `:1334` mislabels Fig 10(b).** The text says "the $h_{\min}=2.27$~\AA{}
tensor of panel (b)", but panel (b) is `lattice_conv_l3`, $h_{\min} = 1.14$ Å.
Confirmed two ways: `plot_fig_63.jl` loads `lattice_conv_l3/V_full_eV.jls`, and
only that campaign reproduces all four quoted $U$ values and the 2.5e-3 symmetry
figure. Table 4's 2.27 Å is correct — it is the level-2 `het3x` run. The paper
legitimately mixes timings from one campaign with the tensor from another, but
must say so.

**(D) The abstract's "1.6 million surface unknowns" matches no §6.4 campaign.**
`conv_l3` batches run at 1,349,640 dof and `het3x` at 651,420. 1.62M is §6.3's
multicube benchmark. Independent of the fixes.

**(F) Fig 10(c) — convergence is far better than published, and the panel needs
rewriting rather than re-plotting:**

| $l_{ec}$ (Å) | $N$ | $U$ published | $U$ new | err published | err new |
|---|---|---|---|---|---|
| 4.545 | 307k | 2.07411 | 2.02698 | 2.3e-2 | **3.8e-4** |
| 2.272 | 652k | 2.05347 | 2.02727 | 1.2e-2 | 2.3e-4 |
| 1.136 | 1.35M | 2.04256 | 2.02745 | 7.1e-3 | 1.5e-4 |
| 0.568 | 2.75M | 2.03646 | 2.02757 | 4.1e-3 | 8.5e-5 |
| 0.284 | 5.56M | 2.03294 | 2.02764 | 2.4e-3 | **5.0e-5** |

Richardson $U_\infty$: 2.02814 → 2.02775, so "$U_\infty \approx 2.03$ eV"
survives verbatim. But *U* is now converged to ~4e-4 at the **coarsest** mesh,
and the published curve was mostly measuring the missing correction rather than
discretization error. Power law still clean (ratio 0.586, slope −0.700 vs −0.779).

**Caption must state the tolerance floor.** The new errors span 3.8e-4 to 5.0e-5
against `gmres_rtol = 1e-5`, so the finest points sit only ~5× above the solver's
own floor. Without saying so, a referee will read the flattening as saturation.

**(F) `:1332` — the on-site $U$ values all move.** Measured on `conv_l3`, the
campaign Fig 10 actually plots:

| | min | max | Si mean | SiO₂ mean |
|---|---|---|---|---|
| published | 1.9876 | 2.1754 | 2.0044 | 2.0962 |
| **new** | **1.9753** | **2.1389** | **1.9906** | **2.0735** |

Every orbital shifts downward by 0.012–0.036 eV. RMS $|\Delta U|/U$ = **0.92%**,
max 1.68%; full-tensor RMS $\Delta V$ = 5.8e-4. The sentence "increases from
about $1.99$ ... to about $2.18$ ... with side-averaged values of $2.00$ and
$2.10$" becomes **1.98 → 2.14, side averages 1.99 and 2.07**.

**The Si/SiO₂ contrast falls 0.0918 → 0.0829 eV, a 9.7% reduction.** The physical
claim survives — the buried junction produces position-dependent screening that
no homogeneous background reproduces — but its magnitude does not. This is the
paper's headline result, so the change should be stated rather than absorbed
silently.

(For scale: the coarser `het3x` campaign shifts 1.6% RMS with a 16% contrast
reduction, so the error is refinement-dependent as expected. Quote the `conv_l3`
figures, which are what Fig 10 shows.)

**(F) `:1333` — the symmetry check improves 36×**, from $2.5\times10^{-3}$ to
**6.8e-5** (`het3x`: 4.33e-3 → 1.11e-4). The paper offers this as "an end-to-end
check of the distributed assembly"; that number was largely reporting the missing
edge correction rather than assembly error. The check is now far more convincing
— this correction strengthens the paper's own argument.

---

## Methods sections

The two fixes are algorithmic and belong in the method description, not only in
the results. §2 currently describes the near-field correction without stating
that touching (edge-adjacent) pairs are corrected by the adaptive quadtree — the
distinction matters, because at $p=6$ with $\tau=10^{-4}$ the Bernstein radius
admits *no* non-touching pair, so the touching correction is the entire near
correction. The no-edge ablation makes this concrete: 5% error and 101 iterations
at $p=6$, versus 4e-4 and 23 with it.

The diagonal preconditioner is not described in the paper at all and now accounts
for roughly half the iteration reduction. It needs a paragraph in §5 (the
algorithm) covering: $A = G + D^{\mathsf T} + C^{\mathsf T}$ with $G$ constant per
interface, right-preconditioning by $P = G^{-1}$, the $|1/t| < 2$ bound, and that
it is a no-op at single contrast.

---

## Reproducibility artifact

`../../BIE_ERI_results/` currently freezes the **submitted** state and is
correct as such. Once §6.4 is settled, it needs: the rerun data swapped in, the
BoundaryIntegral pin moved from `fdb9c62` to `3cc1283`, and the per-directory
READMEs updated to record `correct_edges`/`precondition` per experiment.

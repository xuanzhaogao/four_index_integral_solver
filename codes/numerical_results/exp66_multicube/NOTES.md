# Experiment 6.6 — multi-cube heterojunction, multi-RHS with the graphene p_z orbital

The exp65 **multi-RHS** pipeline run on a multi-dielectric **heterojunction substrate**:
two L×L×L substrate cubes with different permittivities sharing the internal face
x = c_x, with a material slab on top, and the **real graphene Wannier p_z orbital**
(sublattice A, k_323201, from exp65) as the source. The central orbital plus its lattice
neighbours within a cutoff form K pair densities ρ = φ_center·φ_neighbour, block-solved on
ONE shared interface; the central row of the Coulomb matrix V[ρ_11, ρ_1j] is evaluated.
K grows with the cutoff, one K per Slurm-array task (exactly the exp65 array structure).

## System (`run_multicube.jl` + `system_multicube` in `common/Harness.jl`)

Substrate boxes are anchored on the central orbital's density centroid
`c ≈ (0.25, 0, 7.5)` (slab center = orbital centroid; junction directly beneath):

Two substrate variants, selected with `MULTICUBE_GEOM` (**default `sec53`**):

| | `sec53` (default) | `published` |
|---|---|---|
| Ω₁ (cube 1) | 270×270×270, ε = 11.9 | 90×90×90, ε = 11.9 |
| Ω₂ (cube 2) | 270×270×270, ε = 3.9 | 90×90×90, ε = 3.9 |
| slab | 42.95×44.06×9, **ε = 2.4** | 90×90×9, **ε = 10** |
| output tree | `data_sec53/` | `data/` (published), `data_v2/` (July rerun) |

`sec53` is the **same system the article's lattice section (Sec. 5.3) reports**
(`lattice_scale/campaigns/lattice_conv_l3_eps2.4.toml`): the ×3 converged Si|SiO₂
heterojunction and the graphene slab at its cRPA ε = 2.4. `published` reproduces the
geometry behind the submitted Tables 2/3, where ε_slab = 10 was a generic demo value and
L = 90 leaves a finite-size boundary artifact in the on-site U (lattice_scale measured a
bowl at 90 Å; ×3 is converged). Common to both:

| | value | notes |
|---|---|---|
| slab thickness | 9 | the graphene z-support, **not** L/10 — those coincided only at L = 90 |
| anchoring | slab center = orbital centroid ≈ (0, 0, 7.5); junction directly beneath | Sec. 5.3's system translated (cube tops at z = 3, slab z ∈ [3, 12]) |
| outside | vacuum, ε = 1 | |
| sources | graphene p_z orbital + neighbours → K pair densities ρ = φ₁φⱼ | all in the slab |
| eval | central row V[ρ_11, ρ_1j]; V[11] = onsite self-energy | target = ρ_11 = φ₁² core |

The `sec53` slab reuses Sec. 5.3's exact lateral dimensions rather than re-deriving them: the
largest cluster here (cutoff 5 Å) has a 5 Å position half-extent, so with the ~10 Å in-plane
φ² support it sits well inside the 21.5 Å half-width. The lateral size is **fixed across K**
on purpose — an interface that grew with K would confound the multi-RHS scaling study.

Cutoffs `0, 1.5, 2.5, 2.9, 3.9` bohr → **K = 1, 4, 10, 13, 19** (same as exp65). The φ²
support lies fully inside the slab; neighbours spread laterally over both cubes (the
ρ_1j thus feel the asymmetric Si/SiO₂ environment). Special features: buried cube–cube
face, two slab–cube contact faces, triple-junction lines.

**Multi-region screening.** The interface spans 3 ε regions, so exp65's interface-based
`screened_volume_source` (uniform `eps_in` only) does not apply. Each source is screened
**box-based** (`screened_volume_source(boxes, epses, eps_out, …)`, by which box it sits in
— all orbitals → the slab) in both the RHS assembly and the eval. This is the only change
from exp65's pipeline; the batched LHS operator + block GMRES + corrected `pottrg` eval are
identical. RHS is assembled per-source (the batched path loops anyway for distinct-position
sources). Runs at 96 threads (TKM3D's in-function FINUFFT nthreads cap = 16 avoids the FFTW
96-thread plan pathology; no DUCC override needed).

**Quantity.** `V[11]` = screened onsite Coulomb self-energy `⟨ρ_11 | u_inc[ρ_11] + S[σ]⟩`;
`V_vac[11]` = bare interface-free self-energy `⟨ρ_11 | TKM[ρ_11]⟩`; `V_vac/V` ≈ effective
screening of the heterojunction environment. `Vrow` = the full central row over K neighbours.

## Running

```
cd codes/numerical_results

# smoke (fast, validates the full multi-RHS pipeline at one K):
MULTICUBE_SMOKE=1 MULTICUBE_CUTOFF=1.5 \
  JULIA_NUM_THREADS=8 OMP_NUM_THREADS=8 \
  julia --project=. exp66_multicube/scripts/run_multicube.jl

# production K-sweep as parallel Slurm-array tasks (one K per node):
sbatch exp66_multicube/slurm/run_multicube_array.sbatch

# single K on one node (e.g. K=1 onsite):  sbatch exp66_multicube/slurm/run_multicube.sbatch
```

ENV: `MULTICUBE_CUTOFF` (single cutoff per task, the array sets it), `MULTICUBE_SMOKE`,
`MULTICUBE_CUTOFFS` (comma list override), `MULTICUBE_GEOM` (`sec53` default / `published`),
`RERUN_TAG` (output-tree suffix), `CORRECT_EDGES` (default 1), `PRECONDITION` (default 1),
`MULTICUBE_GEOM_ONLY` (print K per cutoff, no solve). Geometry/ε constants are at the top of
`run_multicube.jl`.

Paper-facing consequences of the `sec53` switch (measured Tables 2/3 deltas, the abstract's
*N*, and the new §5.2/§5.3 cross-validation) are in **`PAPER_CHANGES_sec53.md`**.

The Sec. 5.3-system campaigns (write `data_sec53/`, leaving `data/` and `data_v2/` untouched):

```
sbatch exp66_multicube/slurm/run_multicube_array_sec53.sbatch            # K sweep, 7 tasks
sbatch exp66_multicube/slurm/run_multicube_threads_pinned_sec53.sbatch   # thread sweep, 8 tasks
RERUN_TAG=_sec53 julia --project=. plot_scripts/plot_scaling.jl          # -> figs_sec53/
```

Run the K sweep **before** the thread sweep: `plot_scaling.jl` reconstructs the low-thread
block-solve time as `t_RHS + t_1iter * niter_ref`, and `niter_ref` comes from the converged
K=1 record. GMRES needs fewer iterations at ε_slab = 2.4 than at 10, so the `_v2` reference
is stale for this tree.

Outputs: per cutoff `data/raw/multicube_K<K>_cut<cutoff>.jls` (K, n_points, n_src, niter,
block_resid, v11_raw, v11_vac, screen_ratio, Vrow, per-stage timings, peak RSS). In a
single-process sweep (no `MULTICUBE_CUTOFF`) it also writes `data/multicube.csv`. After the
array finishes, gather the `.jls` (cat the per-task `.out` for the summary lines).

## Results

**Single-RHS production** (job 6522501, 96 threads, eps=1e-4, p=4, r=3; the earlier
Harness self-energy path, archived as `data/raw/multicube_orb_L90_eps1e-04_p4_r3.jls`):
N=347,840; GMRES 29 iters; **V=70.30, V_vac=590.3, V_vac/V=8.40**; 81 s.

**Multi-RHS smoke** (`ccmlin078`, 8 threads, K=4 / cutoff 1.5): interface 367,488 pts,
84,872 source pts, GMRES 20 iters; **V[11]=69.2, V_vac/V=8.42** (consistent with the
single-RHS value); block GMRES 116 s (8 threads, smoke tol).

_(production K-array pending)_

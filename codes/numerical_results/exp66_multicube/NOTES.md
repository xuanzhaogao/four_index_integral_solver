# Experiment 6.6 — multi-cube heterojunction, multi-RHS with the graphene p_z orbital

The exp65 **multi-RHS** pipeline run on a multi-dielectric **heterojunction substrate**:
two 90×90×90 substrate cubes with different permittivities sharing the internal face
x = c_x, with a material slab on top, and the **real graphene Wannier p_z orbital**
(sublattice A, k_323201, from exp65) as the source. The central orbital plus its lattice
neighbours within a cutoff form K pair densities ρ = φ_center·φ_neighbour, block-solved on
ONE shared interface; the central row of the Coulomb matrix V[ρ_11, ρ_1j] is evaluated.
K grows with the cutoff, one K per Slurm-array task (exactly the exp65 array structure).

## System (`run_multicube.jl` + `system_multicube` in `common/Harness.jl`)

Substrate boxes are anchored on the central orbital's density centroid
`c ≈ (0.25, 0, 7.5)` (slab center = orbital centroid; junction directly beneath):

| | value | notes |
|---|---|---|
| Ω₁ (cube 1) | 90×90×90, **ε = 11.9** | x ∈ [c_x−90, c_x] (Si) |
| Ω₂ (cube 2) | 90×90×90, **ε = 3.9** | x ∈ [c_x, c_x+90] (SiO₂) |
| slab | 90×90×9, **ε = 10** | centered on the orbital, z ∈ [c_z−4.5, c_z+4.5] |
| outside | vacuum, ε = 1 | |
| sources | graphene p_z orbital + neighbours → K pair densities ρ = φ₁φⱼ | all in the slab |
| eval | central row V[ρ_11, ρ_1j]; V[11] = onsite self-energy | target = ρ_11 = φ₁² core |

Cutoffs `0, 1.5, 2.5, 2.9, 3.9` bohr → **K = 1, 4, 10, 13, 19** (same as exp65). The φ²
support lies fully inside the slab; neighbours spread laterally over both cubes (the
ρ_1j thus feel the asymmetric Si/SiO₂ environment). Special features: buried cube–cube
face, two slab–cube contact faces, triple-junction lines.

**Multi-region screening.** The interface spans 3 ε regions, so exp65's interface-based
`screened_volume_source` (uniform `eps_in` only) does not apply. Each source is screened
**box-based** (`screened_volume_source(boxes, epses, eps_out, …)`, by which box it sits in
— all orbitals → slab ε=10) in both the RHS assembly and the eval. This is the only change
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
`MULTICUBE_CUTOFFS` (comma list override), `CORRECT_EDGES` (default 1),
`MULTICUBE_GEOM_ONLY` (print K per cutoff, no solve). Geometry/ε constants (L=90,
eps 11.9/3.9/10) are at the top of `run_multicube.jl`.

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

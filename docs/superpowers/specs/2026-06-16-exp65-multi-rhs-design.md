# exp65 multi-RHS extension — design

Date: 2026-06-16
Author: Xuanzhao Gao (with Claude Code)
Status: approved for planning

## Context

`exp65_orbital_bench` benchmarks the production four-index dielectric solve on the
real graphene monolayer slab (`L = 90`, `Lz = 3.35`, `eps_in = 3.5` calibrated) with
the real Wannier pz orbital. The current driver `scripts/run_single_rhs.jl` measures
**one** RHS (source `|phi_1|^2`, onsite eval) stage-by-stage with per-stage walltime +
peak-RSS snapshots. Its REPORT explicitly lists as a not-yet-run follow-up:

> Multi-RHS sweep reusing interface + LHS: total time vs number of RHS.

The revised `BoundaryIntegral.jl` (PVF merged to main; multi-RHS + campaign helpers
added) now provides the primitives to do this. This spec extends exp65 to the
multiple-RHS case using **real orbital pairs from a central orbital and its lattice
neighbors** (per `codes/lattice_scale`), not synthetic replication.

## Goal

Quantify how the production solve scales from 1 to K right-hand sides when the
RHS-independent precompute (shared interface mesh + near-corrected layer operator) is
built once and reused, and demonstrate the batched (`nd = K`) FMM / block-GMRES win
over running K independent single-RHS solves. Produce a real (physical) row of the
four-index Coulomb matrix `V[a,b]` for a central orbital.

## Multi-RHS API used (BoundaryIntegral.jl, all `BI.`-prefixed in the script)

- `assemble_lattice_batch(templates, instances, pairs; support_rtol)` → `LatticeBatch`:
  pair densities `rho_ij = phi_i * phi_j` on the **shared union-support grid** (all K
  columns share grid points → exercises the `nd = K` batched FMM), envelope-truncated.
- `solve_dielectric_lattice_batch(boxes, epses, eps_out, b; n_quad, rhs_atol, l_ec,
  fmm_tol, up_tol, max_order, gmres_rtol, max_depth)` → `(; sigma, interface, sources,
  stats)`: refines **one shared interface** on the rss-envelope
  (`multi_dielectric_box3d_rhs_adaptive`), then `solve_dielectric_box3d_block`
  (`batched_lhs_dielectric_box3d_fmm3d_corrected` + `Krylov.block_gmres`).
- `evaluate_batch_potential(interface, Sigma, sources, targets; lhs_tol, volume_tol,
  far_pad)` → `Phi` (n_targets × K): scattered map built once and applied per column;
  incident split near = TKM per source / far = one `nd = K` FMM.
- `four_index_matrix(interface, sources, Sigma; lhs_tol, volume_tol)` → K×K `V[a,b]`
  (independent contraction reference).
- Geometry helpers: `lattice_grid_steps`, `OrbitalInstance`, `density_centroid`,
  `true_cell_vectors`, `batch_volume_sources`, `envelope_volume_source`.

The single-RHS operator `lhs_dielectric_box3d_fmm3d_corrected` (+ `Krylov.gmres`) is
reused for the sequential baseline.

## Geometry — center + neighbor shells (the K sweep)

Real two-orbital graphene monolayer, identical frame to `lattice_scale/demo_2x2.toml`:

- templates: `graphene_00001.xsf` (sublattice A = type 1), `graphene_00002.xsf`
  (sublattice B = type 2), from
  `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15`.
- lattice vectors `a1 = [2.465, 0, 0]`, `a2 = [-1.2325, 2.1347526, 0]` (bohr).
- One **A orbital is the center** (id 1). Neighbors are A/B orbitals at lattice
  translations `n1*a1 + n2*a2`; `OrbitalInstance` offsets via `lattice_grid_steps`.

The K right-hand sides are the pair densities `rho = phi_center * phi_neighbor`
(onsite included) for all neighbors within a cutoff. Sweeping the cutoff grows K
through the physical graphene shells:

| cutoff (bohr) | shell added | K (pairs) |
|---|---|---:|
| 0 (onsite only) | rho_11 = phi_1^2 | 1 |
| ~1.5 | + 3 nearest-neighbor B (d ≈ 1.42) | 4 |
| ~2.5 | + 6 next-nearest A (d ≈ 2.465) | 10 |
| ~2.9 | + 3 third-shell B (d ≈ 2.84) | 13 |
| ~3.9 | + 6 fourth-shell B (d ≈ 3.76) | 19 |

Exact K is whatever the cutoff enumerates; the table is the expected progression.
**K = 1 reproduces the single-RHS benchmark** (the center's `rho_11` is exactly the
single-RHS source) — a built-in validation anchor.

The default production sweep is `CUTOFFS = [0.0, 1.5, 2.5, 2.9, 3.9]` (override via env).
Note: far-shell pair products `phi_center * phi_neighbor` are nearly disjoint, so
`support_rtol` truncation keeps their support (and the envelope/interface) small even
as K grows — the cost growth with K is dominated by the per-column solve and the
K-times eval, not the mesh.

## Parameters

Same production parameters as `run_single_rhs.jl` for apples-to-apples comparability
with the existing REPORT:

| parameter | value |
|---|---|
| slab | `L x L x Lz = 90 x 90 x 3.35`, `eps_in = 3.5`, `eps_out = 1.0` |
| `n_quad` / `edge_level` | 6 / 4 (`l_ec = Lz / 2^edge_level * 1.01`) |
| tolerances | `rhs_tol = 1e-3`, `lhs_tol = 1e-5`, `gmres_rtol = 1e-5` |
| caps | `max_order = 64`, `max_depth = 12` |
| batch | `support_rtol = 1e-4`, `volume_tol = rhs_tol`, `far_pad = 2 grid steps` |

Smoke mode (`ORBBENCH_SMOKE=1`): coarse tolerances (`rhs_tol = 1e-2`, `lhs_tol = 1e-3`,
gmres `1e-3`, `edge_level = 2`), and a reduced cutoff list (e.g. `[0.0, 1.5]`).

## Phases measured (extends the single-RHS phase table)

For each cutoff/K, run the pipeline inline so every stage gets a walltime +
`Sys.maxrss()`/`Base.gc_live_bytes()` snapshot (mirrors `run_single_rhs.jl`):

- **load** — read the 2 xsf templates (one-time, shared across the whole sweep).
- **assemble** — `assemble_lattice_batch`: build the K pair densities. (RHS-set
  dependent; cheap.)
- **precompute** — envelope source + tkm kmax + shared interface build
  (`multi_dielectric_box3d_rhs_adaptive` on the rss-envelope) + **batched LHS operator**
  (`batched_lhs_dielectric_box3d_fmm3d_corrected`, near corrections). *RHS-independent
  across the K columns — the amortized cost.*
- **solve** — batched `nd = K` RHS assembly + **block-GMRES**. *Scales with K — the
  marginal cost.*
- **eval** — `evaluate_batch_potential` at the target densities + contraction into
  `V[a,b]`. TKM near-eval runs at **96 threads as-shipped**, documenting the FFTW
  thread pathology (per decision below).

## Measurements — the full story

1. **Phase split vs K** — the per-stage table at each K; how each phase scales.
2. **Amortization** — fixed precompute vs per-RHS marginal cost `(t_solve + t_eval)/K`;
   total/K vs the single-RHS 296 s anchor; the amortization curve.
3. **Block vs sequential** — on the *same* shared interface + near corrections: one
   `solve_dielectric_box3d_block` (block-GMRES, `nd = K` matvec) vs K independent
   single-RHS GMRES solves (`lhs_dielectric_box3d_fmm3d_corrected` + `Krylov.gmres`,
   `nd = 1`). Report both walltimes and the speedup; this isolates the batched-FMM win.
   The single-RHS operator is built once and reused across the K sequential solves
   (the operator is RHS-independent), so the comparison measures matvec batching, not
   operator-rebuild overhead.
4. **Peak RAM vs K** — the N×K arrays (Sigma, RHS) + FMM workspaces; report peak RSS
   and the per-K increment.
5. **Eval at 96 threads** — report the TKM near-eval cost as-shipped (it hits the
   documented ~50x FFTW 96-thread pathology); flag it explicitly in the report and
   cross-reference `THREAD_PATHOLOGY.md`.

## Method decision

Direct-API **instrumented driver** (`scripts/run_multi_rhs.jl`), not the file-backed
campaign pipeline. Rationale: each stage must be individually timed with RSS snapshots,
and the block-vs-sequential comparison must run on a single shared interface in one
process. The driver calls the same `*_lattice_batch` primitives the campaign uses, so
the numerics are identical to a campaign `solve_batch`/`eval_batch`.

Two deliberate deviations from calling `solve_dielectric_lattice_batch` as a black box,
to keep the benchmark apples-to-apples with `run_single_rhs.jl` (production):

- Build the batched operator **explicitly** via
  `batched_lhs_dielectric_box3d_fmm3d_corrected(interface, lhs_tol, lhs_tol, max_order;
  correct_edges = true)` and call `Krylov.block_gmres` directly — because
  `solve_dielectric_box3d_block` hardcodes `correct_edges = false` and `max_order = 8`,
  whereas single-RHS production uses `correct_edges = true`, `max_order = 64`.
- Use `four_index_matrix(interface, sources, Σ; lhs_tol, volume_tol)` for the
  eval+contraction (the K×K `V` on the shared grid). `evaluate_batch_potential` is the
  campaign's arbitrary-target variant and is not needed for the self-contraction here.

## Decisions (from brainstorming)

- RHS content: **real orbital pairs** (central orbital + lattice neighbors), cutoff sweep.
- Headline: **full story** (amortization + block-vs-sequential + RAM).
- Eval threads: **run at full 96 threads, document the hit** (do not cap FINUFFT
  threads in this benchmark).

## Files

- `exp65_orbital_bench/scripts/run_multi_rhs.jl` — the instrumented driver (mirrors
  `run_single_rhs.jl`: warm-up, inline stages, smoke mode, env knobs).
- `exp65_orbital_bench/slurm/pvf_multi_rhs.sbatch` — copy of `pvf_field.sbatch`:
  `-p gen -C genoa --nodes=1 --exclusive`, juliaup julia
  `/mnt/home/xgao1/.juliaup/bin/julia`, `JULIA_NUM_THREADS = OMP_NUM_THREADS = 96`,
  `--project=/mnt/home/xgao1/codes/BoundaryIntegral.jl`. **User submits `sbatch`.**
- `exp65_orbital_bench/data/raw/multi_rhs_K<k>.jls` — per-K raw record (NamedTuple).
- `exp65_orbital_bench/data/multi_rhs.csv` — one row per K (K, n_points, n_src, niter,
  t_assemble, t_precompute, t_solve_block, t_solve_seq, t_eval, rss_peak_gb, max_rel_asym).
- REPORT addendum in `exp65_orbital_bench/REPORT.md` (table + observations).

No `module load`: pure-Julia/juliaup project (BI.jl env vendors FINUFFT/FMM3D/TKM3D);
mixing Lmod modules with juliaup artifacts risks an ABI mismatch.

## Validation

- **K = 1 onsite** `u` is the same physical quantity as the single-RHS REPORT value
  (~9.88 eV at `eps_in = 3.5`) and should land in the same ~10 eV range — but **exact
  agreement is NOT expected**: this benchmark follows the `lattice_scale/demo_2x2`
  convention (box centered at `cz = 7.5`, orbital at the box **midplane**), whereas
  `run_single_rhs.jl` centers the box at the origin and shifts the orbital to the box
  **top face** (`z = Lz/2`); the REPORT caveat notes this z-placement shifts the
  `eps_in` calibration. Plus the lattice batch uses `support_rtol` envelope truncation +
  grid-snapped placement vs the single-RHS `VolumeSource(tol)` + centroid shift. So the
  K = 1 value is a same-order **sanity check** (report it with the convention caveat);
  the comparable *cost* anchor to the single-RHS REPORT (precompute/solve/eval scale) is
  the firm one.
- **`V[a,b]` symmetry**: `max|V - V'| / max|V|` small (the campaign's `assemble_v`
  diagnostic) — a correctness check on the four-index contraction.
- **Block vs sequential** Sigma agreement to the GMRES tolerance (both solve the same
  system) — confirms the block path is numerically equivalent, not just faster.
- Smoke run completes end-to-end before the production submit.

## Runtime / risk notes

- Production precompute is ~200 s (interface build) regardless of K; the multi-RHS
  point is that it amortizes. The sweep over 5 cutoffs each rebuilds the interface, so
  the full sweep pays ~5 x precompute + the K-dependent solve + eval.
- **The eval phase dominates the wall time at large K.** `four_index_matrix` loops the
  TKM `u_inc` over the K columns, and each `ltkm3dc` call at 96 threads hits the ~50 s
  FFTW pathology (THREAD_PATHOLOGY.md) — so eval ≈ K x ~50 s (K=19 eval ≈ 16 min). This
  is the intended measurement, not a bug to fix here. Summed across the 5 cutoffs the
  job is ~60–90 min, so the sbatch `--time` should be **02:00:00**.
- Peak RAM is expected to stay well under the node's 1.5 TB (single-RHS was 17 GB; the
  N×K arrays add ~0.1 GB per K at production N).
</content>
</invoke>

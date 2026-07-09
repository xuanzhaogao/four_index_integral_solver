# Point-charge screening probe + single-ε control — design

**Date:** 2026-07-09
**Goal:** Determine whether the strong x-ramp in the `onsite_grid` heterojunction U map
(2.17 → 3.41 eV, Si→SiO₂) is *real ε-contrast screening* or a *finite-size / box-geometry
artifact*, using a cheap point-charge probe instead of the full orbital-density pipeline.

## Background

`onsite_grid` puts a graphene-π orbital at each of 1792 sites (z=7.5) over a **dielectric
heterojunction**: two 90³ substrate cubes, ε=11.9 (Si, center x=−39.45) and ε=3.9 (SiO₂,
center x=+50.55), touching at x≈5.55, plus an ε=10 top slab (z∈[3,12]). The on-site U
(diagonal of V, `neighbor_cutoff=0`) ramps monotonically across x and, suspiciously, does
**not** plateau to the single-cube bulk values (SiO₂ side hits 3.41 eV vs 2.25 eV bulk) and
transitions over the whole ~60 Å domain rather than the ~5 Å charge-to-substrate scale.
A symmetric U bowl along **y** (ε constant along that cut) is a genuine, smaller finite-size
edge effect (~2–6 %).

## Observable

Reaction-field (image) self-energy of a unit point charge. At each grid site place q=1 at
z=7.5; solve the BEM for the interface density σ; evaluate **only the induced potential**
Φ_ind at the charge's own location. The direct self-Coulomb diverges for a point but is
~constant across sites, so the induced part isolates exactly the spatially-varying screening.
`U_proxy[a] = Φ_ind(r_a) · 4π·E2` (E2 = 14.3996). Only the **spatial shape** is compared, so
the overall constant (and the ½ self-energy factor, omitted) is irrelevant.

Key implementation fact: in `evaluate_batch_potential`, the term `pottrg * Σ` (line 227,
`laplace3d_pottrg_fmm3d_corrected_hcubature`) **is** Φ_ind; the "incident part" is the direct
term we exclude. So no new BI physics is needed.

## Geometry (two runs, geometry-identical to the het union, ε made uniform)

Substrate: **one box**, center (5.547, 10.318, −42.0), L=(180, 90, 90) — identical footprint
to the touching het union — with uniform ε. Slab: `[5.547, 10.318, 7.5, 90, 90, 9, 10]`
(unchanged). Grid: the same 1792 sites as `onsite_grid` (read x,y,z from its `onsite_U.tsv`).

- **Run A:** substrate ε = 3.9
- **Run B:** substrate ε = 11.9

Solve params copied verbatim from `onsite_grid.toml`: n_quad=6, edge_refine_level=2,
rhs_tol=1e-3, lhs_tol=1e-5, gmres_rtol=1e-5, volume_tol=1e-5, max_order=8, max_depth=128.

## Method — standalone in-process script (no BI edits, no campaign driver)

The interface depends only on the dielectric geometry, so **all charges share one interface
and one LHS operator** — the whole cost saving.

1. Build boxes (uniform substrate + slab) as `BoxGeom`, `epses`, `eps_out=1.0`.
2. Read grid positions `P` (3×N) from `onsite_grid/onsite_U.tsv`.
3. Refinement envelope = `VolumeSource(P, ones, ones)` (unit density at every site). Grid
   spacing ~1.23 ≪ charge-to-substrate gap ~4.5, so this resolves every charge's RHS.
4. `interface = multi_dielectric_box3d_rhs_adaptive(n_quad, l_ec, boxes, epses, env, rhs_tol;
   eps_out, max_depth)` — mirrors `solve_dielectric_lattice_batch`. `l_ec` via the campaign
   rule (min box Lz / 2^level · 1.01).
5. `pottrg = laplace3d_pottrg_fmm3d_corrected_hcubature(interface, P, lhs_tol, lhs_tol, 5.0)`
   once (nt = N targets).
6. **Chunk** the charges (~128/chunk) to bound memory; each chunk shares the interface:
   - sources = single-point `VolumeSource(P[:,a], [1.0], [1.0])` per charge in the chunk.
   - `Σ, _ = solve_dielectric_box3d_block(interface, sources; fmm_tol=lhs_tol, up_tol=lhs_tol,
     max_order, rtol=gmres_rtol, itmax=500, screen_boxes=boxes, screen_epses=epses,
     screen_eps_out=eps_out)` (source screening /ε_local, uniform ε=10 slab → constant).
   - `U_proxy[a] = dot(pottrg[globalidx_a, :], Σ[:, local_a]) · 4π · E2` (diagonal only).
7. Write `figs/pointcharge_ctrl_eps{3.9,11.9}.tsv` (orbital, x, y, z, U_proxy).

## Analysis / decision

Reuse the x/y-cut analysis from the discovery step:
- **Flat in x (only a symmetric edge bowl at both x-ends)** ⇒ a uniform substrate has no
  lateral artifact ⇒ the heterojunction x-ramp is real ε-contrast screening.
- **Still ramps across x** ⇒ finite-size / box-geometry artifact, and the het ramp is
  (partly) not physical.

## Execution

Pilot first: build one interface (ε=11.9) + solve ~5 probe charges (center, far ±x, far ±y),
check timing + sane values, then scale to full 1792. Run on **worker7011** (user-directed).

## Out of scope

Heterojunction point-charge validation run and the ε-sweep single-cube anchor (considered,
deferred — the single-ε control alone answers the artifact question).

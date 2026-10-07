# ×3 heterojunction benchmark redesign + campaign-geometry 3D viz — design

**Date:** 2026-07-10

Two related deliverables, using the finite-size result (a genuine cube needs ~×3 lateral
size to converge; see memory `lattice_scale_onsite_finite_size`).

## Part A — redesigned benchmark `lattice_10x10_het3x`

Replaces `lattice_10x10_multicube` (200 orbitals on a ×1 het, hexagonal parallelogram cloud).

**Substrate — ×3 heterojunction (converged):**
```
Si   cube (ε=11.9): center (-129.453, 10.318, -132),  270×270×270
SiO₂ cube (ε=3.9) : center ( 140.547, 10.318, -132),  270×270×270
slab      (ε=10)  : center (   5.547, 10.318,  7.5),  270×270×9   (z∈[3,12])
junction plane x = 5.547
```

**Orbitals — ~198, real graphene honeycomb in a square window:** the true hexagonal lattice
(a₁=(2.465,0), a₂=(−1.2325,2.1347526), both sublattices A/B) clipped to a SQUARE window of
half-width 11.5 bohr (side ~23) centered on the junction (5.547, 10.318), z=7.5 → 99 orbitals
over Si, 99 over SiO₂. (An earlier "square super-lattice + 2-atom basis" gave dimerized columns,
not a natural lattice; switched to real honeycomb per review.)

**Slab is auto-sized:** Lx/Ly = orbital-block extent + 2×11 bohr (measured φ² support) so the
density stays inside ε=10 (no leak); thickness Lz=9 unchanged. Substrate cubes stay 270³.

**Pairing/solve:** neighbor_cutoff = 5.0 (V_ijkl for pairs within 5 bohr); solve params identical
to existing campaigns (n_quad=6, edge_refine_level=2, rhs_tol=1e-3, lhs_tol=1e-5, gmres_rtol=1e-5,
support_rtol=1e-4, volume_tol=1e-5, max_order=8, max_depth=128); n_centers_per_batch as in the
current benchmark. Output = full four-index V_ijkl via the campaign pipeline.

**Deliverable:** `scripts/gen_benchmark_het3x.jl` → `campaigns/lattice_10x10_het3x.toml`
(mirrors gen_hetline_campaign.jl / gen_onsite_campaign.jl).

## Part B — `plot_campaign_geometry` in BoundaryIntegral MakieExt

A reusable 3D geometry plot over the campaign data model (boxes/epses/eps_out/orbitals) —
consolidates the ad-hoc `viz_onsite_system.jl` (2D) and `plot_lattice_system_3d.jl`.

- Stub `function plot_campaign_geometry end` in BI core; implementation in `ext/MakieExt.jl`.
- `plot_campaign_geometry(c::CampaignInput; full=true, kwargs...) -> Figure`, plus a
  `plot_campaign_geometry(toml_path::AbstractString)` convenience (loads the campaign).
- Draws: every dielectric box as a **wireframe cuboid** (12 edges) colored by ε (colormap +
  colorbar / ε labels); orbital positions as a 3D scatter **colored by sublattice `type`**.
  `Axis3`, equal aspect (`aspect=:data`), full extent to scale (cubes large, orbital block a
  small cluster — intended), title = campaign name.
- Returns the `Figure`; caller saves.

**Implementation constraint:** BI edits are made in a **git worktree** (never the live checkout),
tested there, then merged to `main` so the live checkout / consumer env picks it up cleanly.

**Follow-up (out of scope now):** retire `viz_onsite_system.jl` / `plot_lattice_system_3d.jl`
once the package function is in.

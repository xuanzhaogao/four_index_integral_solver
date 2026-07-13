# Design: onsite-U vs. orbital position on the heterojunction substrate

Date: 2026-06-24

## Goal

Map the screened onsite Hubbard interaction $U_i \equiv V[\rho_{ii},\rho_{ii}]$ (the
screened Coulomb self-energy of the orbital density $\rho_{ii}=\phi_i^2$) as a function
of the orbital's position on the exp66/lattice Si/SiO2 heterojunction substrate. The
$10\times10$ V_ijkl campaign sampled only the central ~third of the slab, so its onsite-U
map (Fig.~12 right of the article) shows a gradient but never the bulk-Si / bulk-SiO2
plateaus or the slab-edge behavior. This experiment fills the slab and answers:

1. How does $U$ vary with position across the buried Si/SiO2 junction?
2. Does $U$ reach flat bulk plateaus in the interior before the slab edges — i.e. is the
   90³-bohr substrate large enough to define a bulk onsite U?

## Approach

Reuse the campaign's `prepare` + `solve` with **`neighbor_cutoff = 0`** (one self-pair
$(i,i)$ per site → one σ-solve per site, distributed one worker per node). Then **replace
`consolidate/eval/assemble`** — which build the full $n\times n$ matrix by evaluating each
site's Φ at the *shared union of all supports* — with a custom **onsite eval**
(`scripts/onsite_eval.jl`): for each solved batch, evaluate Φ of $\rho_{ii}$ at the
orbital's **own** support and contract → $U_i = \langle\rho_{ii}\,|\,\Phi_i\rangle$ (eV).

This is the onsite-ONLY cost the experiment actually needs: **1792 σ-solves + 1792 LOCAL
evals** (each at the site's own ~few-k support), *not* a 1792×1792 full-matrix campaign.
`onsite_eval.jl` calls only public BI functions (`evaluate_batch_potential`,
`grid_positions`, `VolumeSource`, …) — no BI changes — and is resumable (appends to
`<root>/onsite_U.tsv`, skips finished sites). The eV normalization matches `assemble_v`'s
onsite formula, $U_i^{\rm eV} = \langle\rho_{ii}|\Phi_i\rangle \cdot 4\pi E_2/\|\phi_i\|^4$
with $E_2 = 14.3996$ eV·Å and $\|\phi_i\|^2 = \int\rho_{ii}$.

Orbitals are placed on **real graphene lattice points** (not arbitrary positions): since
we skip the off-diagonal V_ijkl, we can afford a large number of onsite solves and a
dense, full-resolution map.

## System (unchanged from the 10×10 campaign)

Current 90³-bohr multicube heterojunction substrate: a Si cube ($\epsilon=11.9$) and a
SiO2 cube ($\epsilon=3.9$) sharing the buried junction at $x=5.547$, capped by a
$90\times90\times9$ host slab ($\epsilon=10$) spanning $z\in[3,12]$,
$x\in[-39.453,50.547]$, $y\in[-34.682,55.318]$. The orbital probe is one graphene $p_z$
Wannier orbital (sublattice A) from the campaign templates, translated to each lattice
site at the slab mid-plane $z=7.5$. Solve parameters identical to `lattice_10x10`
(`n_quad=6`, `edge_refine_level=2`, `rhs_tol=1e-3`, `lhs_tol=gmres_rtol=1e-5`,
`support_rtol=1e-4`, `volume_tol=1e-5`, `max_order=8`, `max_depth=128`),
`n_centers_per_batch = 1`.

## "Avoid touching the surface" — placement constraint

Every orbital's support must stay strictly inside the slab so the box-based screening
(orbital → slab $\epsilon=10$) holds and no density crosses a dielectric face. The φ²
support of the Wannier template (`graphene_00001.xsf`) was **measured** at the campaign's
`support_rtol = 1e-4` (centroid $z=7.5$): in-plane radius **≈ 10.0 bohr**, vertical
half-extent **≈ 3.83 bohr** (3D max radius 10.7).

- **Vertical:** $z=7.5$ (slab mid-plane). The 3.83-bohr half-extent clears the slab
  top/bottom faces (4.5 bohr away) by **≈ 0.67 bohr** — tight but valid at `support_rtol
  = 1e-4`. (Note: at a tighter `1e-6` truncation the half-extent is 7.34 bohr and would
  not clear; we stay at `1e-4`, consistent with the campaign.)
- **Lateral:** clip the supercell so every site is ≥ **margin = 11 bohr** (in-plane
  support 10 + 1 buffer) inside the slab $x,y$ faces. The generator computes this margin
  from the measured support and rejects sites outside the clipped footprint.

This placement was validated visually before the run by `scripts/viz_onsite_system.jl`
(→ `figs/fig_onsite_system.pdf`): a top view (coverage/margin/junction) and an $(x,z)$
side view (the φ²-band clearance to the slab faces).

**Resolution caveat:** because $U(\mathbf r)$ is the orbital convolved with the dielectric
environment, the junction step is smeared over ~the orbital footprint (~10–20 bohr); the
map varies only on ~10-bohr scales. Full graphene density (2.465 bohr) therefore
oversamples it ~16×, but is used deliberately for a crisp map / exact lattice registry.

## Phases

**Phase 1 — 1D line (`onsite_line` campaign).** One lattice row crossing the junction:
sites at ~fixed $y$ (nearest lattice row to the slab-center $y=10.318$), spanning the
margin-clipped slab width in $x$. `prepare` + `solve` + `onsite_eval` give $U_i$ along the
row (`onsite_U.tsv`) → $U(x)$.
Deliverable: a $U(x)$ curve with the junction $x=5.547$ marked, used to read off the
Si/SiO2 plateaus, transition width, and edge upturn — and to judge substrate size.

**Phase 2 — 2D map (`onsite_grid` campaign).** Full graphene density ($a=2.465$ bohr)
filling the 11-bohr-margin-clipped slab footprint: **1792 sites** (896 over Si / 896 over
SiO₂, symmetric about the junction, both sublattices), confirmed by the pre-flight viz.
`prepare` + `solve` + `onsite_eval` give $U(x,y)$ (`onsite_U.tsv`). Deliverable: a
$U(x,y)$ heatmap/scatter over the lattice.

## New artifacts (in `lattice_scale`)

- `scripts/viz_onsite_system.jl` (**done**): pre-flight viz of the substrate + margin-clipped
  probe lattice (top + side views) → `figs/fig_onsite_system.pdf`; shares the lattice/margin
  logic that `gen_onsite_campaign.jl` will reuse.
- `scripts/gen_onsite_campaign.jl` (**done**): generate `campaigns/onsite_line.toml` (28
  sites) and `campaigns/onsite_grid.toml` (1792 sites) — the margin-clipped lattice, both
  sublattices, same `[dielectrics]`/`[solve]` blocks as `lattice_10x10_multicube.toml`,
  `[pairing] neighbor_cutoff = 0`.
- `scripts/onsite_eval.jl` (**done**): per-batch onsite eval — Φ of $\rho_{ii}$ at the
  site's own support → $U_i$ (eV) into `<root>/onsite_U.tsv`; resumable. Replaces
  `consolidate/eval/assemble`.
- `scripts/plot_onsite.jl` (**done**): read `onsite_U.tsv` and plot $U(x)$ (1D) /
  $U(x,y)$ (2D) in the shared `fig_gen` style → `lattice_scale/figs/`.
- `jobscripts/run_onsite.sbatch` (**done**): `prepare → solve → onsite_eval`, resumable.
- Campaign roots on ceph (as for `lattice_10x10`), one per phase.

## Execution

`jobscripts/run_onsite.sbatch` runs `prepare → solve → onsite_eval` (user submits;
resumable, one worker per node for `solve`). Cost: **1792 σ-solves + 1792 local evals** —
each onsite eval evaluates Φ at the site's own support (~few-k targets), so there is no
full-matrix work. The distributed `solve` is the bottleneck; `onsite_eval` is the cheap
post-pass (runs on the head node looping batches; can be split to a `--nodes=1` follow-up
job to free the solve nodes for the big grid). Phase 1 (28 sites, 1 node) is the pilot —
**run it first**, read the per-site solve+eval cost from the logs, then size Phase 2's node
count. `[batching] n_centers_per_batch` can be raised to amortize the solves (multi-RHS
block solve over nearby sites). A future optimization is to fold the onsite eval into the
`solve` worker (one small BI change) so it is fully distributed with no second pass.

## Deliverables & out of scope

Deliverables: `U(x)` and `U(x,y)` figures, and a verdict on whether the 90³ substrate is
large enough (plateaus reached or not). Out of scope: enlarging the substrate (only a
follow-up if Phase 1 shows no interior plateau) and any off-diagonal $V_{ijkl}$.
Whether to add these results to the article is a separate, later decision.

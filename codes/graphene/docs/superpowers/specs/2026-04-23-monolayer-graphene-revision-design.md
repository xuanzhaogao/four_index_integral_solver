# Monolayer Graphene Revision — Design

Date: 2026-04-23

## Background

Previous work in this repo targeted a bilayer-slab dielectric model (scalar `eps_in = 2.4`, `eps_out = 1.0`, `Lz = 6.7`) with shell labels `U_00..U_05` aligned to Rösner et al. for bilayer graphene. The 2026-03-26 report explicitly notes that this is a dielectric-slab approximation, not a reproduction of the paper's momentum-dependent `ε_eff^2D(q)`.

A collaborator (Malte Rösner) has now produced a proper **monolayer** graphene QE+CoQui calculation that provides:

- Two Wannier XSF orbital files on a 150×150×192 grid over a ~12.24 × 10.60 × 14.92 Å supercell (primitive cell `a = 2.465 Å`, `c = 15 Å`, two C atoms on sublattices A and B separated by 1.4232 Å in xy at `z = 7.5 Å`).
- Reference Coulomb matrix elements (orbital-average):
  - **Bare**: intra = 17.434191 eV, inter = 8.839615 eV, Hund's (spin-flip) = 0.130804 eV, Hund's (pair-hopping) = 0.130804 eV.
  - **Static cRPA screened**: intra = 9.872342 eV, inter = 5.236802 eV, Hund's = 0.093121 eV each.

Malte has indicated that the bare channel is already testable. The screened reference depends on an effective dielectric "thickness" + effective ε(q) that he is still writing a script to extract. The bilayer monolayer example from him is expected later.

## Scope of This Revision

This spec covers **only the monolayer bare-interaction validation**, organized as a new parallel workflow sibling to the archived bilayer-slab work. Screened monolayer and bilayer revisions are out of scope here and will get their own specs once Malte delivers the missing reference artifacts.

## Repo Reorganization

1. Move the existing bilayer-slab workflow under a `bilayer_slab/` subdirectory:
   - `graphene/{src,scripts,data,results}` → `graphene/bilayer_slab/{src,scripts,data,results}`
   - Move existing bilayer-related test files from `graphene/test/` into `graphene/bilayer_slab/test/`.
   - Keep `Project.toml`, `Manifest.toml`, `docs/` at `graphene/` top level (shared environment and design docs).
   - Keep `graphene/test/runtests.jl` at the top level as a **dispatcher**: it `include`s `bilayer_slab/test/runtests.jl` and `monolayer/test/runtests.jl` so `julia --project=. test/runtests.jl` still runs the full suite. (Each subdirectory gets its own `test/runtests.jl`.)
   - Update the 2026-03-26 report: add a "Historical — bilayer-slab approximation" header noting it has been superseded for the purpose of reference comparison by the monolayer workflow.
2. Create a new `graphene/monolayer/` tree:
   - `graphene/monolayer/src/`
   - `graphene/monolayer/scripts/`
   - `graphene/monolayer/data/`
   - `graphene/monolayer/results/`
   - `graphene/monolayer/test/` (with its own `runtests.jl`)
3. Script entry points under the existing `Project.toml` environment (no new Julia project) so both workflows share dependencies.

Any path references in the moved bilayer scripts are updated to the new nested layout.

## Bare-Interaction Targets

Four quantities, matched against the cRPA reference (orbital-average):

| Name | Expression | Reference (eV) |
|------|------------|---------------:|
| `U_onsite` | `∫∫ \|φ₁(r)\|² \|φ₁(r')\|² / \|r-r'\|`                         | 17.434191 |
| `V_nn`     | `∫∫ \|φ₁(r)\|² \|φ₂(r')\|² / \|r-r'\|`                         | 8.839615  |
| `J_sf`     | `∫∫ φ₁(r)φ₂(r) · φ₂(r')φ₁(r') / \|r-r'\|` (spin-flip exchange)   | 0.130804  |
| `J_ph`     | `∫∫ φ₁(r)φ₂(r) · φ₁(r')φ₂(r') / \|r-r'\|` (pair-hopping)        | 0.130804  |

Notes:

- `φ₁, φ₂` are the two Wannier orbitals read from `graphene_00001.xsf`, `graphene_00002.xsf`.
- The A↔B nearest-neighbor pair is taken directly as the intra-cell pair in the new XSF (atoms at `(0,0,7.5)` and `(0, 1.4232, 7.5)` — distance `1.4232 Å`). No translated orbital copies are constructed for the first validation run.
- `J_sf` and `J_ph` are mathematically equal for real orbitals — CoQui confirms this (both `0.130804 eV`). We still compute both independently as a consistency check rather than hard-coding the identity.
- "Orbital-average" in the CoQui output means the intra- and inter-orbital numbers are averaged over the two Wannier orbitals; for graphene's two sublattice orbitals related by symmetry the two channels coincide, so comparing our single-pair computation to the average is consistent.

## Numerical Method

**Primary method: direct real-space volume integration** (option a in brainstorming).

- Reuse the existing `volume_integral` + source-density machinery used for the bilayer-slab `u_int` path.
- Each integrand is a *pair of real-space scalar fields* on the XSF grid:
  - density-density channels (`U_onsite`, `V_nn`): source and target fields are `|φ|²`;
  - Hund's channels (`J_sf`, `J_ph`): source and target fields are the **signed products** `φ₁φ₂` (these integrate to approximately zero, unlike a density).
- No dielectric solver, no GMRES, no surface density — only the direct part.
- A single shared centering shift is applied (same convention as the bilayer workflow) so that both orbitals sit at the box center for numerical safety.

**Signed-orbital pipeline.** The bilayer-slab workflow squared the orbitals at load time. The monolayer Hund's channels need signed products, so the loader is generalized:

- `load_orbital_xsf(path; square::Bool)` — if `square=true` returns `|φ|²` (density channel), if `square=false` returns the signed `φ` grid.
- `build_pair_source(φ_a, φ_b, channel::Symbol)` — for `:density_density` returns `(|φ_a|², |φ_b|²)`, for `:hund_exchange`/`:hund_pairhop` returns `(φ_a·φ_b, φ_b·φ_a)` and `(φ_a·φ_b, φ_a·φ_b)` respectively.
- Downstream volume-integral code treats the two fields symmetrically as "source" and "target" scalar fields; it does not need to know whether they are densities or signed products.

**Convergence check: mirror-padding in xy.** Since the monolayer XSF supercell is ~5×5 primitive cells in plane rather than the 90×90 Å isolated-molecule box used in the bilayer-slab work, we verify that the orbital tail has decayed at the xy boundary:

- Recompute each of the four bare integrals with the orbital mirror-padded once (and once more, if time permits) in x and y.
- Report the sensitivity as a diagnostic alongside the primary values. Expected: sub-percent drift for well-localized p_z Wannier orbitals.

**Fallback.** If mirror-padding shifts any value by more than ~2%, switch to an FFT-based Poisson solve on the periodic grid with a Gygi-Baldereschi madelung correction matching CoQui's treatment (option b in brainstorming). That fallback is not implemented in this spec — it gets its own spec if triggered.

## Components

New code lives under `graphene/monolayer/src/`:

1. **`MonolayerOrbitalLoader.jl`** — reads the two monolayer XSF files, returns signed `φ` grids and (optionally) `|φ|²`. Applies the shared centering shift. Exposes grid + lattice metadata.
2. **`MonolayerBareIntegrals.jl`** — given a pair `(f_source, f_target)` of scalar fields on the XSF grid, computes the bare `1/|r-r'|` volume integral using the existing `volume_integral` machinery. Returns raw units and eV.
3. **`MonolayerBareWorkflow.jl`** — top-level glue: loads orbitals, constructs the four channel pairs, calls the integral routine, writes results.
4. **`scripts/bare_monolayer_graphene.jl`** — run entry point (no CLI options for v1; all parameters are read from the script).
5. **`test/`** — unit tests covering:
   - loader correctness (grid shape, centering, signed vs squared);
   - pair-source construction for each channel;
   - a synthetic validation (two Gaussian densities with known analytic Coulomb integral) for the bare-integral routine;
   - reproducibility of `J_sf == J_ph` for real orbitals to within numerical tolerance.

## Data Flow

```
graphene_00001.xsf  ─┐
                     ├─► MonolayerOrbitalLoader ─► (φ₁, φ₂, grid metadata)
graphene_00002.xsf  ─┘
                                                    │
                                                    ▼
                   channel ∈ {onsite, nn, hund_sf, hund_ph}
                                                    │
                                                    ▼
                          build_pair_source(...) ─► (f_src, f_tgt)
                                                    │
                                                    ▼
                    MonolayerBareIntegrals ─► U_channel [raw, eV]
                                                    │
                                                    ▼
                 data/bare_monolayer_graphene.csv  (+ results/report.md)
```

## Output

- `graphene/monolayer/data/bare_monolayer_graphene.csv` columns:
  `channel, source_label, target_label, u_raw, u_ev, u_ref_ev, rel_err_pct, n_source_points, n_target_points, source_norm, target_norm, source_tol, mirror_pad_level, shift_x, shift_y, shift_z`
- `graphene/monolayer/results/2026-04-23-bare-monolayer-report.md` summarizing:
  - reference block copied from `_coqui_thc_crpa.out`;
  - our four values + relative errors;
  - mirror-pad convergence table;
  - pass/fail judgment against a chosen tolerance (see below).

## Success Criterion

Primary: all four bare channel values within **≤ 2%** relative error of the CoQui reference at mirror-pad level 1. Stretch: ≤ 1%. Anything outside 5% triggers the fallback method rather than a pass.

For the Hund's channels (~0.13 eV), 2% is ~2.6 meV — tight but achievable given the integrals use products of well-localized orbitals.

## Verification Plan

- Unit tests pass locally (`julia --project=. test/runtests.jl`).
- Synthetic Gaussian-pair test reproduces analytic `erf`-based Coulomb integral to the solver tolerance.
- Density-channel values sit inside the 2% band.
- Mirror-pad drift is sub-percent.
- Report committed alongside the CSV.

## Explicitly Out of Scope

- Screened (cRPA) monolayer comparison — deferred until Malte provides the effective-thickness / ε(q) extraction script.
- Bilayer monolayer re-run — deferred until Malte delivers the bilayer example.
- FFT + madelung fallback — only implemented if the convergence check fails.
- Additional shells beyond onsite/NN — no reference data yet.
- Any changes to the bilayer-slab workflow beyond moving it under `bilayer_slab/` and adding a historical header.

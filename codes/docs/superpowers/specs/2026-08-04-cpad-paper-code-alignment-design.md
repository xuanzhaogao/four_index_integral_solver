# Aligning the near/far incident-potential geometry with Section 3 of the paper

Date: 2026-08-04

## Problem

Section 3 of `~/Articles/four_indices_bie/main.tex` specifies a single near/far
split for the incident potential. `BoundaryIntegral.jl` implements it three
different ways, none of which matches the paper exactly, and the Section 3.4
validation experiment (Fig. 5) sidesteps the whole question by setting
`h_n = 0` and evaluating only targets inside the source box.

### The paper's specification

With source samples on a lattice with basis matrix `A_ρ`:

| quantity | paper | equation |
|---|---|---|
| characteristic spacing | `h = ‖A_ρ‖₂` | (3.3) |
| source box | `B`, axis-aligned, side lengths `l_α` | (3.7) |
| padding width | `h_n = c_pad · h` | (3.10) |
| near region | `B_pad = B ⊕ [-h_n, h_n]³`; target is near iff `x ∈ B_pad` | (3.9) |
| truncation radius | `L = sqrt(Σ_α (l_α + h_n)²)` | (3.11) |
| real-space period | `L_α ≥ l_α + h_n + L`, i.e. `Δk_α = 2π/(l_α + h_n + L)` | (3.16) |

`L` is the largest source-to-target distance: a source is confined to `B` and a
target to `B_pad`, so the per-axis reach is `l_α/2 + (l_α/2 + h_n) = l_α + h_n`.

`k_max` is deliberately **not** pinned by the paper. §3.3 states `k_Nyq = π/h`
only "for an isotropic Cartesian grid", and says that for a general lattice
`k_Nyq` "is determined by that discretization and is regarded as known here".
`k_max` is therefore out of scope for this work.

### The four mismatches

1. **`h` is the wrong end of the spectrum.** `_estimate_source_spacing`
   (`src/shape/box3d_fmm_helpers.jl:57-74`) returns the *minimum* nearest-neighbour
   distance. Eq. (3.3) asks for `‖A_ρ‖₂`, the largest singular value, described as
   "a conservative characteristic spacing". These agree for a cubic grid and
   diverge for the skewed graphene/xsf lattices, where the code's `h` — and hence
   `h_n` — comes out too small. `VolumeSource` (`src/core/sources.jl:20-24`) stores
   only `positions`/`weights`/`density`; the basis passed to the skew-lattice
   constructor (`src/core/sources.jl:170`) is discarded, so `‖A_ρ‖₂` is not
   recoverable at the point of use.

2. **`c_pad` has two different values.** `margin_h = 5.0` in
   `PrecomputedVolumeField` and `h_factor = 5.0` in `_classify_near_far_targets`,
   but `far_pad_steps = 2.0` (`src/campaign/toml_input.jl:75`) in the campaign path
   that produced the Section 6.4 results.

3. **One near region is a ball, not a box.** `_classify_near_far_targets`
   (`src/shape/box3d_fmm_helpers.jl:105-113`) classifies via a KDTree `inrange`
   query of radius `h * h_factor` around the source points, not via Eq. (3.9)'s box.

4. **`L` and the period are over-estimated.** `PrecomputedVolumeField`
   (`src/shape/volume_field.jl:106-108`) takes `L` to be the diagonal of `B_pad`
   (`sqrt(Σ(l_α + 2h_n)²)`) and the period to be `(l_α + 2h_n) + L`, as if sources
   could sit anywhere in `B_pad` rather than only in `B`. Also, it pads the
   bounding box of the *sample points*, so its effective `h_n` is `(c_pad - ½)h`.
   At `n = 64` in the Fig. 5 setup this puts the shipped default at `η ≈ 1.06`
   rather than exactly 1. Every deviation is in the safe direction, so the Lemma
   still holds — but the numbers do not reproduce the paper.

## Decisions

- Unify all three paths on the paper's geometry with `c_pad = 5`.
- `h = ‖A_ρ‖₂` for `h_n`, obtained by storing the basis on `VolumeSource`.
  `k_max` keeps `_estimate_source_spacing`, which the paper leaves open.
- Define `B` as the **cell box** of the samples: the sample bounding box inflated
  by half a lattice cell, i.e. per-axis half-extent `½ Σ_j |A_ρ[α,j]|`. This is
  the box the quadrature cells tile, it contains the full support of the
  piecewise-constant density, and it makes the Fig. 5 experiment's `B = [-1,1]³`
  with `l = 2` literally correct for a cell-centered grid. It is also strictly
  larger than the sample bounding box, so it is the safe choice.
- Fig. 5 keeps its two panels showing the same two quantities. Only `h_n` and the
  target set change.
- Panel (b)'s `η` uses the paper's formula. No remark about implementation
  headroom is needed, because after this change the code sits at exactly `η = 1`.

## Design

### Component 1 — `VolumeSource` carries its lattice basis

`src/core/sources.jl`. Add one field:

```julia
struct VolumeSource{T, D} <: AbstractSource
    positions::Matrix{T}
    weights::Vector{T}
    density::Vector{T}
    lattice_basis::Union{Nothing, NTuple{3, NTuple{3, T}}}
end
```

The three inner tuples are the lattice vectors `a₁, a₂, a₃`, matching the
`basis` argument of `sources.jl:170` and the `p = origin + u*a + v*b + w*c`
convention of `_volume_source_positions` (`sources.jl:120`). They are the
*columns* of `A_ρ` in Eq. (3.2), so `‖A_ρ‖₂` is the largest singular value of the
3×3 matrix formed by stacking them as columns.

Populated by the skew-lattice constructor (`sources.jl:170`) and by both xsf
constructors (`src/utils/xsf_reader.jl:142`, `:226`) from `true_cell_vectors`.
`nothing` for sources built from bare point lists, where no lattice exists.

Two accessors, new file `src/core/source_geometry.jl`:

- `lattice_spacing(vs)` → `‖A_ρ‖₂` when `lattice_basis` is present, else falls
  back to `_estimate_source_spacing(vs)` with a `@debug` note. This is the `h` of
  Eq. (3.3).
- `source_box(vs)` → `(lo, hi, l)` for the cell box `B`, using the basis
  half-extents when present and `h/2` per axis otherwise.

Every existing `VolumeSource(...)` constructor keeps its current signature and
passes `nothing` or the derived basis, so no caller changes. `screened_volume_source`
and `batch_volume_sources` must propagate the field rather than rebuilding
without it — this is the main place a silent regression could hide.

### Component 2 — one shared geometry helper

New in `src/core/source_geometry.jl`, so all three paths compute the same numbers
from the same code:

```julia
"""
    near_field_geometry(vs; c_pad = 5.0) -> (lo, hi, l, h, hn, L, dks)

Section 3 near/far geometry: `B` = cell box of the samples (Eq. 3.7),
`hn = c_pad * ‖A_ρ‖₂` (Eqs. 3.3, 3.10), `B_pad = [lo, hi] = B ⊕ hn` (Eq. 3.9),
`L = sqrt(Σ(l_α + hn)²)` (Eq. 3.11), `dks[α] = 2π/(l_α + hn + L)` (Eq. 3.16).
"""
```

and `in_near_region(lo, hi, targets, i)` implementing Eq. (3.9)'s box test.
`dks` uses `prevfloat` as today, to stay strictly inside the aliasing-free set.

### Component 3 — the three call sites

**`PrecomputedVolumeField`** (`src/shape/volume_field.jl:83-140`). Rename
`margin_h` → `c_pad` (default 5.0); replace lines 102-108 with a
`near_field_geometry` call; `center` becomes the centre of `B`, which is also the
centre of `B_pad`, so the source and target phase shifts at lines 114-116 and in
`_field_scaled_targets` stay consistent. Pass the new `L` to
`TKM3D.truncated_laplace3d_hat`. `in_field_box` (`:304`) delegates to
`in_near_region`. FINUFFT's coordinate range stays safe: the largest scaled
target coordinate is `2π(l/2 + h_n)/(l + h_n + L) < π` because `h_n < L`.

**`evaluate_batch_potential`** (`src/solver/lattice_batch.jl:198-278`). Replace
the required absolute `far_pad::Float64` with `c_pad::Float64 = 5.0`, and lines
232-240 with `near_field_geometry` + `in_near_region`.

Its near branch currently calls `TKM3D.ltkm3dc`, which cannot be made to follow
the paper: `src/continuous.jl:133-146` derives `lengths, center` from
`combined_box_geometry_3xn(src, trg)` and sets `L = sqrt(l_x² + l_y² + l_z²)`
internally, exposing only `kmax` and `eps`. There is no hook for an externally
supplied `L` or `Δk`. This branch therefore routes through
`PrecomputedVolumeField` instead, leaving exactly one implementation of the
Fourier geometry.

Two consequences:

- **`L` stops depending on the target set.** `ltkm3dc` builds its box from the
  actual targets passed, so today's batch `L` varies with whichever targets land
  in the batch. Under Eq. (3.11) `L` depends only on `B_pad`, so results become
  reproducible independent of batch composition. This is a behaviour change in the
  production path, additional to the `c_pad` 2 → 5 change.
- **Memory must not scale with `K`.** `PrecomputedVolumeField` is per-density, and
  the batch has `K` densities on shared positions. Building `K` fields at once
  would hold `K` coefficient arrays, which is untenable at production scale
  (the struct docstring cites ~16 GB for a single field with gradients). So split
  the type: `near_field_geometry` returns a small geometry struct that
  `PrecomputedVolumeField` holds, and the batch loop constructs, evaluates, and
  discards one field per `a in 1:K` sequentially, reusing the shared geometry.
  Peak memory then matches today's single `ltkm3dc` call. Pass
  `compute_grad = false` — this path needs only the potential.

**`_classify_near_far_targets`** and `_classify_near_far_panels`
(`src/shape/box3d_fmm_helpers.jl:85-113`). Replace the KDTree ball with the
Eq. (3.9) box, `h_factor` → `c_pad`. Drop the now-unused `NearestNeighbors`
import if nothing else needs it.

**Campaign plumbing.** `far_pad_steps` → `c_pad` in `src/campaign/toml_input.jl:39,75`
(default 5.0) and `src/campaign/tasks.jl:311-320`, which no longer needs
`max_step`. Update the three fixtures under `test/fixtures/`.

### Component 4 — Fig. 5 data script

`fig_gen/fig5_data.jl`. Drop the direct `TKM3D.ltkm3dc` / `FMM3D.lfmm3d` calls
in favour of `PrecomputedVolumeField(vs; tol, c_pad = 5.0)` +
`volume_field_potential`, which now performs the paper's classification and the
FMM far branch internally.

Target set, fixed across every `n` and both methods:

- 200 **near**: uniform in `[-0.9, 0.9]³`, inside `B_pad` for all `n`.
- 200 **far**: `max|coord| ∈ [2.5, 4]`, outside `B_pad` for all `n` — the
  coarsest `n = 8` has `h_n = 1.25` and `B_pad ≈ [-2.25, 2.25]³`.

Holding the classification fixed across `n` is what keeps the curves clean; a
target at `|x| = 1.2` would otherwise switch branches between `n = 8` and
`n = 64`. Reported error is the combined relative `L²` over all 400 targets.

Panel (a): the hybrid at `τ ∈ {10⁻³, 10⁻⁶, 10⁻⁹, 10⁻¹²}` against the pure
particle sum, both on all 400 targets. The far branch also converges
spectrally — a Gaussian on a uniform grid with a smooth kernel has vanishing
Euler–Maclaurin boundary terms — so the existing narrative holds, and the
particle-everywhere curve still stalls near `10⁻⁴` on account of its near targets.

Panel (b): a standalone tunable-`Δk` hybrid, because the production path now
enforces `η = 1` exactly and cannot be driven below it. Its geometry mirrors
`near_field_geometry` with `Δk_α = 2π/(l_α + h_n + L)/η`, at `n = 64`, on the same
400 targets. Far targets are `η`-independent and set the post-threshold plateau.
`η = min_α L_α/(l_α + h_n + L)`; the `= min_α L_α/(l_α + L)` simplification in
Eq. (3.30) goes away.

`fig5_plot.jl` needs only the panel-(b) x-label and the hardcoded `[4:17]` slice
adjusted.

### Component 5 — paper edits

`main.tex`, all in §3.4 and the Fig. 5 caption:

- `\paragraph{TKM convergence}` (1277-1288): `h_n = 0` → `h_n = c_pad h`,
  `c_pad = 5`; describe the mixed target set.
- Lines 1231-1233 and 1296-1297: delete the two sentences disclaiming that the
  near–far switch and `c_pad` are untested. They now are.
- `\paragraph{Periodization threshold}` (1305-1318): `h_n = 0` → `h_n = 5h`;
  Eq. (3.30) loses its simplification.
- Caption (1332-1342): new target set, `h_n = c_pad h`.

No changes to §3.1-3.3 — the code is moving to them, not the reverse.

## Staging

The work is cohesive but too large to land atomically. Three stages, each with a
gate, so a failure in one does not strand the others:

- **Stage A — `BoundaryIntegral.jl`.** Components 1-3 plus their tests, in a git
  worktree. Gate: full test suite green, and the cross-path equivalence and
  cubic-lattice-unchanged tests passing. Nothing downstream starts until this
  merges.
- **Stage B — Fig. 5 and the paper.** Components 4-5. Gate: panel (b) reproduces
  the periodization threshold and panel (a) shows the plateau structure at all
  four tolerances.
**No Stage C.** The Section 6.4 campaign is not rerun. The quoted numbers were
produced at `far_pad_steps = 2`, and the new default is `c_pad = 5`, so in
principle they correspond to a superseded configuration. In practice the change
moves targets between two evaluators that already agree far below the quoted
precision: `test/solver/lattice_batch.jl:88` asserts `< 1e-5` maximum relative
difference between the near and far branches on a well-resolved source at
`far_pad = 2h`, against §6.4's `U ≈ 2.0–2.2 eV` and symmetry residual `0.0043`.
Revisit only if Stage A's tests reveal a discrepancy larger than that on a
realistic lattice.

## Verification

Order matters; each step gates the next.

1. `Pkg.test("BoundaryIntegral")` in the worktree. Expect real changes in
   `test/shape/volume_field.jl` (the `η ≈ 1.06` → `1.0` shift moves mode counts
   and therefore tolerances) and in `test/solver/lattice_batch.jl:88-151`, whose
   `far_pad = 2.0 * h` becomes `c_pad = 5`. Tolerances may need loosening or
   tightening; each adjustment must be justified, not fitted.
2. An equivalence test asserting the three paths agree on classification and on
   `(l, h_n, L, dks)` for one skewed and one cubic lattice. This is the guard
   against the three implementations drifting apart again.
3. A regression check that `‖A_ρ‖₂` and the old min-nearest-neighbour `h` agree
   on a cubic grid, so cubic-lattice results are provably unchanged.
4. Fig. 5 smoke run at `n ∈ {8, 16}` and `τ = 10⁻⁶` only, confirming the near and
   far target groups classify as intended at both ends of the `n` sweep, before
   the full sweep.
5. Full `fig5_data.jl`, then `fig5_plot.jl`, then check the transition in panel (b).

   Expect the observed transition slightly *below* `η = 1`, and do not tune
   anything to move it. Eq. (3.16) is a sufficient condition derived from
   `supp(ρ) ⊆ B`, but the test Gaussian has `s = 0.10`, so its numerical support is
   a ball of radius ≈ 0.6 rather than the full `B = [-1,1]³`. The periodic image of
   `supp(G_L * ρ)` therefore clears the target region a little early: with targets
   in `[-0.9, 0.9]³` the crossing sits near `η ≈ 0.91`. The claim the panel
   supports is that `η ≥ 1` is sufficient and that `O(1)` contamination appears
   below the threshold — not that the transition is exactly at 1. The caption
   wording must match what the data shows.

### Bounding the un-rerun Section 6.4 change

Since §6.4 is not being regenerated, Stage A must produce evidence that the
change is below its quoted precision rather than assuming it. Add a test that
evaluates the same realistic lattice batch at `c_pad = 2` and `c_pad = 5` and
records the maximum relative difference in `Φ`. If that number is not comfortably
below `0.0043`-level significance, the no-rerun decision has to be revisited.

## Risks

- **§6.4 stands on un-regenerated numbers.** `c_pad` 2 → 5 enlarges the near
  region, moving targets from the FMM branch to the TKM branch, and `L` also
  stops depending on batch composition. The argument that this is immaterial rests
  on the two branches agreeing to `< 1e-5`; the bounding test above is what turns
  that from an assumption into a measurement.
- **The `VolumeSource` field addition is the widest blast radius.** Any
  constructor or transformation that rebuilds a `VolumeSource` without
  propagating `lattice_basis` silently reverts `h` to the fallback. The
  propagation audit of `screened_volume_source` / `batch_volume_sources` /
  `envelope_volume_source` is mandatory, not optional.
- **Reducing `L` and the period reduces conservatism.** This is intended — it is
  what matching Eq. (3.16) means, and it cuts the mode array by roughly 16% — but
  it removes headroom that was previously absorbing any error in the `B_pad`
  containment argument. Test 2 is what protects this.
- **`volume_field.jl`'s `cache_fft` path depends on eight TKM3D private
  functions** (enumerated at `src/shape/volume_field.jl:151-163`). Mode-count
  changes flow into `nfdim`, so the `cache_fft` testset at
  `test/shape/volume_field.jl:140` must be run, not skipped.

## Process notes

All `BoundaryIntegral.jl` edits happen in a git worktree, never the live
checkout — an in-flight edit there breaks precompilation for concurrent runs.
Consumer projects path-dev `BoundaryIntegral`, so after the `VolumeSource` change
`fig_gen` will need `Pkg.resolve()` if its Manifest goes stale.

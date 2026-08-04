# c_pad / Paper-Code Alignment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make BoundaryIntegral.jl's near/far incident-potential geometry match Section 3 of `~/Articles/four_indices_bie/main.tex` exactly on all three code paths at `c_pad = 5`, then rework the Fig. 5 validation to use it with `h_n = 5h` and mixed near+far targets.

**Architecture:** One new file, `src/core/source_geometry.jl`, owns the paper's geometry: `lattice_spacing` (Eq. 3.3), `source_box` (Eq. 3.7), and `near_field_geometry` returning a `NearFieldGeometry` struct with `(lo, hi, center, l, h, hn, L, dks)` per Eqs. (3.9), (3.11), (3.16). All three existing near/far implementations delegate to it, so the geometry exists in exactly one place. `VolumeSource` gains a `lattice_basis` field so `‖A_ρ‖₂` is available at the point of use.

**Tech Stack:** Julia. `BoundaryIntegral.jl` (dev-pathed at `~/codes/BoundaryIntegral.jl`), `TKM3D.jl` (dev-pathed), `FMM3D`, `FINUFFT`. `CairoMakie` for figures. Tests via `Pkg.test("BoundaryIntegral")`.

## Global Constraints

- **`c_pad = 5.0`** is the single default everywhere. It replaces `margin_h = 5.0`, `h_factor = 5.0`, and `far_pad_steps = 2.0`.
- **`h = ‖A_ρ‖₂`** (largest singular value of the primitive sampling-cell basis) is the spacing used for `h_n = c_pad · h`. It is NOT used for `k_max`, which keeps `_estimate_source_spacing` (min nearest-neighbour distance) because §3.3 deliberately leaves `k_Nyq` open for general lattices.
- **`A_ρ` column `j` = `(axes[j] step) × basis[j]`.** The `basis` argument of `VolumeSource(axes, weights, density, origin, basis)` holds full cell vectors while `axes` hold fractional coordinates (see `src/utils/xsf_reader.jl:145-155`), so the primitive step vector must be reconstructed by multiplying by the axis step.
- **`B` is the cell box:** sample bounding box inflated by half a lattice cell, per-axis half-extent `½ Σ_j |A_ρ[α,j]|`.
- **`L = sqrt(Σ_α (l_α + h_n)²)`** and **`Δk_α = 2π/(l_α + h_n + L)`**, with `prevfloat` applied to `Δk_α` to stay strictly inside the aliasing-free set.
- **All BoundaryIntegral.jl edits happen in a git worktree**, never `~/codes/BoundaryIntegral.jl` directly — an in-flight edit there breaks precompilation for concurrent runs.
- **Do not loosen a test tolerance to make it pass.** Every tolerance change must be justified by an identified mechanism, recorded in the commit message.
- Spec: `docs/superpowers/specs/2026-08-04-cpad-paper-code-alignment-design.md`.

---

## File Structure

**Stage A — `~/codes/BoundaryIntegral.jl` (in a worktree)**

| file | responsibility |
|---|---|
| `src/core/source_geometry.jl` (create) | The paper's Section 3 geometry, and nothing else: `NearFieldGeometry`, `lattice_spacing`, `source_box`, `near_field_geometry`, `in_near_region` |
| `src/core/sources.jl` (modify) | `VolumeSource.lattice_basis` field, `_primitive_basis`, basis propagation through every constructor and `screened_volume_source` |
| `src/shape/volume_field.jl` (modify) | `PrecomputedVolumeField` delegates geometry; `margin_h` → `c_pad` |
| `src/shape/box3d_fmm_helpers.jl` (modify) | KDTree ball → Eq. (3.9) box; `h_factor` → `c_pad` |
| `src/solver/dielectric_box3d.jl` (modify) | caller rename `h_factor` → `c_pad` |
| `src/solver/lattice_batch.jl` (modify) | `LatticeBatch.lattice_basis`; `evaluate_batch_potential` `far_pad` → `c_pad`, near branch via `PrecomputedVolumeField` |
| `src/campaign/toml_input.jl`, `src/campaign/tasks.jl` (modify) | `far_pad_steps` → `c_pad` plumbing |
| `test/core/source_geometry.jl` (create) | geometry unit tests + cross-path equivalence + cubic-unchanged + `c_pad` 2-vs-5 bound |

**Stage B — `~/work/four_index_integral_solver/codes` and `~/Articles/four_indices_bie`**

| file | responsibility |
|---|---|
| `fig_gen/fig5_data.jl` (modify) | hybrid evaluator, `c_pad = 5`, mixed near+far targets |
| `fig_gen/fig5_plot.jl` (modify) | panel (b) x-label and `eta_list` slice |
| `main.tex` (modify) | §3.4 and the Fig. 5 caption |

---

## Stage A — BoundaryIntegral.jl

### Task 0: Create the worktree

**Files:** none (environment setup)

**Interfaces:**
- Produces: a worktree path, referred to below as `$WT`. All Stage A tasks run inside it.

- [ ] **Step 1: Create the worktree and branch**

```bash
cd ~/codes/BoundaryIntegral.jl
git worktree add ~/codes/bi-cpad-align -b cpad-paper-alignment
cd ~/codes/bi-cpad-align
git status
```

- [ ] **Step 2: Confirm the baseline test suite passes before changing anything**

```bash
cd ~/codes/bi-cpad-align
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()' 2>&1 | tail -40
```

Expected: all testsets pass. If anything fails on a clean checkout, STOP and report — you cannot attribute later failures without a green baseline.

Record the wall-clock time of this run; later tasks re-run the same command.

---

### Task 1: `VolumeSource` carries its lattice basis

**Files:**
- Modify: `$WT/src/core/sources.jl:20-24` (struct), `:106-118` (flat ctor), `:157-182` (grid ctors), `:205-217` (point-list ctor), `:343`, `:361`, `:470` (`screened_volume_source`)
- Test: `$WT/test/core/source_geometry.jl` (create)
- Modify: `$WT/test/runtests.jl` (register the new test file)

**Interfaces:**
- Produces:
  - `VolumeSource{T,3}` gains field `lattice_basis::Union{Nothing, NTuple{3, NTuple{3, T}}}`, holding the **primitive sampling-cell** basis `A_ρ` as three column vectors, or `nothing` when the source is not a lattice.
  - `VolumeSource{T,D}(positions, weights, density)` — 3-arg form still valid, sets `lattice_basis = nothing`.
  - `VolumeSource(positions, weights, density; tol = 0, lattice_basis = nothing)` — flat ctor gains the kwarg.
  - `BoundaryIntegral._primitive_basis(axes::NTuple{3,Vector{T}}, basis::NTuple{3,NTuple{3,T}}) -> Union{Nothing, NTuple{3,NTuple{3,T}}}`
  - `BoundaryIntegral.with_density(vs::VolumeSource{T,3}, density::Vector{T}) -> VolumeSource{T,3}` — copies positions/weights/basis, swaps density.

- [ ] **Step 1: Write the failing test**

Create `$WT/test/core/source_geometry.jl`:

```julia
using Test
using LinearAlgebra
using BoundaryIntegral
const BI = BoundaryIntegral

@testset "VolumeSource lattice_basis" begin
    # --- identity-basis grid ctor: axes hold ABSOLUTE coords with step h ---
    n = 8
    h = 2.0 / n
    xs = collect(-1.0 + h/2 .+ h .* (0:n-1))
    w  = fill(h^3, n, n, n)
    d  = fill(1.0, n, n, n)
    vs = VolumeSource((xs, xs, xs), w, d)
    @test vs.lattice_basis !== nothing
    # A_rho = h * I
    for j in 1:3, i in 1:3
        @test isapprox(vs.lattice_basis[j][i], i == j ? h : 0.0; atol = 1e-14)
    end

    # --- skew ctor: axes hold FRACTIONAL coords, basis holds full cell vectors ---
    nx = 4
    frac = collect((i - 1) / nx for i in 1:nx)
    At = (2.0, 0.0, 0.0); Bt = (1.0, 2.0, 0.0); Ct = (0.0, 0.0, 3.0)
    wg = fill(1.0, nx, nx, nx); dg = fill(1.0, nx, nx, nx)
    vsk = VolumeSource((frac, frac, frac), wg, dg, (0.0, 0.0, 0.0), (At, Bt, Ct))
    # primitive step vectors are the cell vectors divided by nx
    @test isapprox(collect(vsk.lattice_basis[1]), collect(At) ./ nx; atol = 1e-14)
    @test isapprox(collect(vsk.lattice_basis[2]), collect(Bt) ./ nx; atol = 1e-14)
    @test isapprox(collect(vsk.lattice_basis[3]), collect(Ct) ./ nx; atol = 1e-14)

    # --- flat ctor from bare points: no lattice ---
    pts = rand(3, 10)
    vsf = VolumeSource(pts, fill(1.0, 10), fill(1.0, 10))
    @test vsf.lattice_basis === nothing

    # --- flat ctor accepts an explicit basis ---
    vsb = VolumeSource(pts, fill(1.0, 10), fill(1.0, 10);
                       lattice_basis = ((h, 0.0, 0.0), (0.0, h, 0.0), (0.0, 0.0, h)))
    @test vsb.lattice_basis !== nothing

    # --- with_density preserves the basis ---
    vs2 = BI.with_density(vs, fill(2.0, length(vs.density)))
    @test vs2.lattice_basis == vs.lattice_basis
    @test all(vs2.density .== 2.0)
    @test vs2.positions == vs.positions

    # --- non-uniform axes are not a lattice ---
    bad = [0.0, 0.1, 0.5, 1.0]
    vsn = VolumeSource((bad, bad, bad), fill(1.0, 4, 4, 4), fill(1.0, 4, 4, 4))
    @test vsn.lattice_basis === nothing
end

@testset "screened_volume_source preserves lattice_basis" begin
    gsrc = BI.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 8, 1e-6)
    @test gsrc.lattice_basis !== nothing
    boxes = [(center = (0.0, 0.0, 0.0), Lx = 1.0, Ly = 1.0, Lz = 1.0)]
    sc = BI.screened_volume_source(boxes, [2.0], 1.0, gsrc, BI.SharpScreening())
    @test sc.lattice_basis == gsrc.lattice_basis
end
```

Register it in `$WT/test/runtests.jl` alongside the other includes (match the surrounding style, e.g. `@testset "source geometry" begin include("core/source_geometry.jl") end` — read the file and follow its existing convention exactly).

- [ ] **Step 2: Run the test to verify it fails**

```bash
cd ~/codes/bi-cpad-align
julia --project=. -e 'using Pkg; Pkg.test(; test_args = ["source_geometry"])' 2>&1 | tail -30
```

If the suite does not support filtering, run `julia --project=. test/runtests.jl` and read the `source geometry` testset output.

Expected: FAIL — `type VolumeSource has no field lattice_basis`, and `UndefVarError: with_density`.

- [ ] **Step 3: Add the field and the basis helper**

In `$WT/src/core/sources.jl`, replace the struct at lines 20-24:

```julia
struct VolumeSource{T, D} <: AbstractSource
    positions::Matrix{T}
    weights::Vector{T}
    density::Vector{T}
    # Primitive sampling-cell basis A_rho (three column vectors), or `nothing` when
    # the source is not a lattice. This is the paper's Eq. (3.2) basis: column j is
    # the step vector along lattice direction j, i.e. (axis step) * (cell vector).
    lattice_basis::Union{Nothing, NTuple{3, NTuple{3, T}}}
end

# 3-arg form: no lattice information.
VolumeSource{T, D}(positions::Matrix{T}, weights::Vector{T}, density::Vector{T}) where {T, D} =
    VolumeSource{T, D}(positions, weights, density, nothing)
```

Add after `_is_uniform_axis` (currently ending at line 203):

```julia
"""
    _primitive_basis(axes, basis) -> Union{Nothing, NTuple{3, NTuple{3, T}}}

Primitive sampling-cell basis `A_rho` (Eq. 3.2): column `j` is `(axis j step) * basis[j]`.
`axes` may hold fractional or absolute coordinates; `basis` holds the corresponding
cell vectors. Returns `nothing` unless all three axes are uniform with >= 2 points,
in which case the samples do not form a lattice.
"""
function _primitive_basis(
    axes::NTuple{3, Vector{T}},
    basis::NTuple{3, NTuple{3, T}},
) where {T}
    steps = ntuple(3) do d
        ax = axes[d]
        (length(ax) >= 2 && _is_uniform_axis(ax)) ? (ax[2] - ax[1]) : nothing
    end
    any(isnothing, steps) && return nothing
    return ntuple(d -> ntuple(i -> T(steps[d]) * basis[d][i], 3), 3)
end

"""
    with_density(vs, density) -> VolumeSource

Copy of `vs` with `density` replaced, preserving positions, weights, and
`lattice_basis`. Use this instead of rebuilding through the flat constructor,
which would silently drop the lattice basis.
"""
function with_density(vs::VolumeSource{T, 3}, density::Vector{T}) where {T}
    length(density) == length(vs.density) ||
        throw(ArgumentError("density length must match the source"))
    return VolumeSource{T, 3}(copy(vs.positions), copy(vs.weights), density, vs.lattice_basis)
end
```

- [ ] **Step 4: Thread the basis through every constructor**

`sources.jl:106-118`, flat ctor — add the kwarg:

```julia
function VolumeSource(
    positions::AbstractMatrix{T},
    weights::AbstractVector{T},
    density::AbstractVector{T};
    tol::Real = 0,
    lattice_basis::Union{Nothing, NTuple{3, NTuple{3, T}}} = nothing,
) where {T}
    _validate_volume_source_flat(positions, weights, density)
    pos = Matrix{T}(positions)
    w = Vector{T}(weights)
    rho = Vector{T}(density)
    pos_t, w_t, rho_t = _truncate_volume_source(pos, w, rho, T(tol))
    return VolumeSource{T, 3}(pos_t, w_t, rho_t, lattice_basis)
end
```

`sources.jl:157-164`, identity-basis grid ctor — pass the derived basis:

```julia
function VolumeSource{T, 3}(axes::NTuple{3, Vector{T}}, weights::Array{T, 3}, density::Array{T, 3}; tol::Real = 0) where {T}
    _validate_volume_source_grid(axes, weights, density)
    origin = (zero(T), zero(T), zero(T))
    basis = _identity_basis(T, Val(3))
    positions = _volume_source_positions(axes, origin, basis)
    w, rho = _flatten_volume_source_data(weights, density)
    return VolumeSource(positions, w, rho; tol = tol,
                        lattice_basis = _primitive_basis(axes, basis))
end
```

`sources.jl:170-182`, skew ctor — same, with the supplied basis:

```julia
function VolumeSource(
    axes::NTuple{3, Vector{T}},
    weights::Array{T, 3},
    density::Array{T, 3},
    origin::NTuple{3, T},
    basis::NTuple{3, NTuple{3, T}};
    tol::Real = 0,
) where {T}
    _validate_volume_source_grid(axes, weights, density)
    positions = _volume_source_positions(axes, origin, basis)
    w, rho = _flatten_volume_source_data(weights, density)
    return VolumeSource(positions, w, rho; tol = tol,
                        lattice_basis = _primitive_basis(axes, basis))
end
```

The point-list ctor at `sources.jl:205-217` stays as-is; bare points have no lattice.

- [ ] **Step 5: Fix the three `screened_volume_source` overloads**

Replace the trailing `return VolumeSource(copy(vs.positions), copy(vs.weights), rho)` at `sources.jl:343`, `:361`, and `:470` with:

```julia
    return with_density(vs, rho)
```

Truncation is not applied by these overloads, so positions are preserved and the basis remains valid.

- [ ] **Step 6: Run the new test to verify it passes**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | tail -40
```

Expected: the two new testsets PASS, and every pre-existing testset still passes. `GaussianVolumeSource` routes through the identity-basis grid ctor, so `gsrc.lattice_basis` is populated without touching it.

If any pre-existing test fails here, it is a basis-propagation regression, not a tolerance issue — find the constructor that dropped the field.

- [ ] **Step 7: Commit**

```bash
cd ~/codes/bi-cpad-align
git add src/core/sources.jl test/core/source_geometry.jl test/runtests.jl
git commit -m "$(cat <<'EOF'
sources: VolumeSource carries its primitive lattice basis

Eq. (3.3) needs h = ||A_rho||_2, but the basis passed to the skew-lattice
constructor was discarded. Store it as the primitive sampling-cell basis
(axis step times cell vector) so it is available at the point of use.

with_density replaces the rebuild-through-flat-constructor pattern in the
three screened_volume_source overloads, which silently dropped the basis.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: The shared Section 3 geometry helper

**Files:**
- Create: `$WT/src/core/source_geometry.jl`
- Modify: `$WT/src/BoundaryIntegral.jl` (include + exports)
- Test: `$WT/test/core/source_geometry.jl` (append)

**Interfaces:**
- Consumes: `VolumeSource.lattice_basis`, `_estimate_source_spacing` (from `src/shape/box3d_fmm_helpers.jl:57`)
- Produces:
  - `struct NearFieldGeometry{T}` with fields `lo, hi, center, l :: NTuple{3,T}` and `h, hn, L :: T` and `dks :: NTuple{3,T}`
  - `lattice_spacing(vs::VolumeSource{T,3}) -> T` — Eq. (3.3)
  - `source_box(vs::VolumeSource{T,3}) -> (lo, hi, l)` each `NTuple{3,T}` — Eq. (3.7), the cell box
  - `near_field_geometry(vs::VolumeSource{T,3}; c_pad::Real = 5.0) -> NearFieldGeometry{T}`
  - `in_near_region(g::NearFieldGeometry, targets::AbstractMatrix, i::Integer) -> Bool` — Eq. (3.9)

- [ ] **Step 1: Write the failing test**

Append to `$WT/test/core/source_geometry.jl`:

```julia
@testset "near_field_geometry matches Section 3 by hand" begin
    # Cell-centered cubic grid on B = [-1,1]^3, so l = 2 exactly.
    n = 64
    h = 2.0 / n
    xs = collect(-1.0 + h/2 .+ h .* (0:n-1))
    w  = fill(h^3, n, n, n)
    d  = fill(1.0, n, n, n)
    vs = VolumeSource((xs, xs, xs), w, d)

    @test isapprox(BI.lattice_spacing(vs), h; rtol = 1e-14)

    lo, hi, l = BI.source_box(vs)
    for a in 1:3
        @test isapprox(lo[a], -1.0; atol = 1e-14)
        @test isapprox(hi[a],  1.0; atol = 1e-14)
        @test isapprox(l[a],   2.0; atol = 1e-14)
    end

    c_pad = 5.0
    g = BI.near_field_geometry(vs; c_pad = c_pad)
    hn_exp = c_pad * h
    L_exp  = sqrt(3 * (2.0 + hn_exp)^2)
    @test isapprox(g.hn, hn_exp; rtol = 1e-14)
    @test isapprox(g.L,  L_exp;  rtol = 1e-14)
    for a in 1:3
        @test isapprox(g.lo[a], -1.0 - hn_exp; atol = 1e-14)
        @test isapprox(g.hi[a],  1.0 + hn_exp; atol = 1e-14)
        @test isapprox(g.center[a], 0.0; atol = 1e-14)
        # Eq. (3.16) at equality
        @test isapprox(g.dks[a], 2π / (2.0 + hn_exp + L_exp); rtol = 1e-12)
        # prevfloat keeps us strictly inside the aliasing-free set
        @test g.dks[a] <= 2π / (2.0 + hn_exp + L_exp)
    end

    # Eq. (3.9) box test
    tg = [0.0  1.0 + hn_exp/2   1.0 + 2*hn_exp;
          0.0  0.0              0.0;
          0.0  0.0              0.0]
    @test BI.in_near_region(g, tg, 1)
    @test BI.in_near_region(g, tg, 2)
    @test !BI.in_near_region(g, tg, 3)
end

@testset "lattice_spacing is the 2-norm, not the shortest step" begin
    # Skew lattice: ||A_rho||_2 must exceed the longest column norm.
    nx = 4
    frac = collect((i - 1) / nx for i in 1:nx)
    At = (1.0, 0.0, 0.0); Bt = (0.9, 0.4, 0.0); Ct = (0.0, 0.0, 1.0)
    vsk = VolumeSource((frac, frac, frac), fill(1.0, nx, nx, nx), fill(1.0, nx, nx, nx),
                       (0.0, 0.0, 0.0), (At, Bt, Ct))
    A = hcat(collect.(collect(vsk.lattice_basis))...)
    @test isapprox(BI.lattice_spacing(vsk), opnorm(A, 2); rtol = 1e-12)
    @test BI.lattice_spacing(vsk) > maximum(norm.(collect.(collect(vsk.lattice_basis))))
end

@testset "lattice_spacing falls back for non-lattice sources" begin
    # Cubic point cloud with no basis: fallback must equal the grid spacing,
    # which is what guarantees cubic-lattice results are unchanged.
    n = 6; h = 0.25
    pts = Matrix{Float64}(undef, 3, n^3); m = 0
    for k in 1:n, j in 1:n, i in 1:n
        m += 1
        pts[1, m] = i * h; pts[2, m] = j * h; pts[3, m] = k * h
    end
    vsf = VolumeSource(pts, fill(h^3, n^3), fill(1.0, n^3))
    @test vsf.lattice_basis === nothing
    @test isapprox(BI.lattice_spacing(vsf), h; rtol = 1e-12)
end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | tail -30
```

Expected: FAIL — `UndefVarError: lattice_spacing` / `source_box` / `near_field_geometry`.

- [ ] **Step 3: Write the geometry module**

Create `$WT/src/core/source_geometry.jl`:

```julia
# src/core/source_geometry.jl
#
# The near/far geometry of Section 3 of the four-index BIE paper, in one place.
# All three near/far call sites (PrecomputedVolumeField, evaluate_batch_potential,
# _classify_near_far_targets) delegate here so the geometry cannot drift apart.
#
# Equation references are to ~/Articles/four_indices_bie/main.tex:
#   (3.3)  h   = ||A_rho||_2, a conservative characteristic spacing
#   (3.7)  B   = cuboidal source box, side lengths l_alpha
#   (3.10) h_n = c_pad * h
#   (3.9)  B_pad = B + [-h_n, h_n]^3; a target is near iff x in B_pad
#   (3.11) L   = sqrt(sum_alpha (l_alpha + h_n)^2), the largest source-target distance
#   (3.16) L_alpha >= l_alpha + h_n + L, i.e. dk_alpha = 2pi / (l_alpha + h_n + L)

"""
    NearFieldGeometry{T}

Section 3 near/far geometry for one volume source. `lo`/`hi` are the corners of
the padded near region `B_pad` (Eq. 3.9), `center` is the centre of `B` (which is
also the centre of `B_pad`), `l` the side lengths of `B` (Eq. 3.7), `h` the
characteristic spacing (Eq. 3.3), `hn = c_pad * h` (Eq. 3.10), `L` the truncation
radius (Eq. 3.11), and `dks` the Fourier spacings (Eq. 3.16 at equality).
"""
struct NearFieldGeometry{T}
    lo::NTuple{3, T}
    hi::NTuple{3, T}
    center::NTuple{3, T}
    l::NTuple{3, T}
    h::T
    hn::T
    L::T
    dks::NTuple{3, T}
end

"""
    lattice_spacing(vs) -> T

The characteristic spacing `h = ||A_rho||_2` of Eq. (3.3): the largest singular
value of the primitive sampling-cell basis. Falls back to
`_estimate_source_spacing` (minimum nearest-neighbour distance) when the source
carries no lattice basis; the two agree for a cubic grid.

This is the spacing for `h_n = c_pad * h` only. `k_max` keeps
`_estimate_source_spacing`, because Section 3.3 leaves `k_Nyq` open for a
general sampling lattice.
"""
function lattice_spacing(vs::VolumeSource{T, 3}) where {T}
    b = vs.lattice_basis
    b === nothing && return _estimate_source_spacing(vs)
    A = T[b[1][1] b[2][1] b[3][1];
          b[1][2] b[2][2] b[3][2];
          b[1][3] b[2][3] b[3][3]]
    s = T(opnorm(A, 2))
    return s > zero(T) ? s : _estimate_source_spacing(vs)
end

"""
    source_box(vs) -> (lo, hi, l)

The source box `B` of Eq. (3.7): the bounding box of the samples inflated by half
a lattice cell, so that the quadrature cells tile `B` and `B` contains the full
support of the piecewise-constant density. The per-axis half-extent is
`(1/2) * sum_j |A_rho[alpha, j]|`, or `h/2` when no basis is available.
"""
function source_box(vs::VolumeSource{T, 3}) where {T}
    pos = vs.positions
    size(pos, 2) >= 1 || throw(ArgumentError("source_box requires at least one source point"))
    b = vs.lattice_basis
    half = if b === nothing
        hh = _estimate_source_spacing(vs) / 2
        ntuple(_ -> hh, 3)
    else
        ntuple(a -> (abs(b[1][a]) + abs(b[2][a]) + abs(b[3][a])) / 2, 3)
    end
    lo = ntuple(a -> minimum(view(pos, a, :)) - half[a], 3)
    hi = ntuple(a -> maximum(view(pos, a, :)) + half[a], 3)
    l  = ntuple(a -> hi[a] - lo[a], 3)
    return lo, hi, l
end

"""
    near_field_geometry(vs; c_pad = 5.0) -> NearFieldGeometry

Assemble the Section 3 geometry for `vs`. `dks` uses `prevfloat` so the real-space
period is strictly greater than the Eq. (3.16) bound, keeping the trapezoidal sum
inside the aliasing-free set rather than exactly on its boundary.
"""
function near_field_geometry(vs::VolumeSource{T, 3}; c_pad::Real = 5.0) where {T}
    c_pad >= 0 || throw(ArgumentError("c_pad must be >= 0"))
    loB, hiB, l = source_box(vs)
    h = lattice_spacing(vs)
    h > zero(T) || throw(ArgumentError("lattice spacing must be positive"))
    hn = T(c_pad) * h
    lo = ntuple(d -> loB[d] - hn, 3)
    hi = ntuple(d -> hiB[d] + hn, 3)
    center = ntuple(d -> (loB[d] + hiB[d]) / 2, 3)
    L = sqrt((l[1] + hn)^2 + (l[2] + hn)^2 + (l[3] + hn)^2)
    dks = ntuple(d -> prevfloat(T(2π) / (l[d] + hn + L)), 3)
    return NearFieldGeometry{T}(lo, hi, center, l, h, hn, T(L), dks)
end

"""
    in_near_region(g, targets, i) -> Bool

Eq. (3.9): is column `i` of the `3 × n` matrix `targets` inside `B_pad`?
"""
@inline in_near_region(g::NearFieldGeometry, targets::AbstractMatrix, i::Integer) =
    (g.lo[1] <= targets[1, i] <= g.hi[1]) &&
    (g.lo[2] <= targets[2, i] <= g.hi[2]) &&
    (g.lo[3] <= targets[3, i] <= g.hi[3])
```

- [ ] **Step 4: Wire it into the package**

In `$WT/src/BoundaryIntegral.jl`, add the include. `source_geometry.jl` calls `_estimate_source_spacing`, which lives in `src/shape/box3d_fmm_helpers.jl`; Julia resolves method calls at run time, not include time, so the include order does not matter as long as both are included. Place it immediately after the `src/core/sources.jl` include to keep `core/` together.

Add to the export list, next to the existing `PrecomputedVolumeField` export (`src/BoundaryIntegral.jl:23-24`):

```julia
export NearFieldGeometry, near_field_geometry, in_near_region, lattice_spacing, source_box
```

Confirm `LinearAlgebra` (for `opnorm`) is already imported at the top of `src/BoundaryIntegral.jl`; it is used by `lattice_batch.jl:138` via `det`, so it should be. If it is only imported inside a submodule, add `using LinearAlgebra: opnorm` to the main module.

- [ ] **Step 5: Run to verify it passes**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | tail -40
```

Expected: the three new testsets PASS. All pre-existing testsets still pass — nothing consumes the new code yet.

- [ ] **Step 6: Commit**

```bash
cd ~/codes/bi-cpad-align
git add src/core/source_geometry.jl src/BoundaryIntegral.jl test/core/source_geometry.jl
git commit -m "$(cat <<'EOF'
core: single source of truth for the Section 3 near/far geometry

near_field_geometry implements Eqs. (3.3), (3.7), (3.9), (3.11), (3.16)
directly: h = ||A_rho||_2, B = the cell box of the samples, B_pad = B + h_n,
L = the largest source-target distance, dk_alpha = 2pi/(l_alpha + h_n + L).

No call site uses it yet; the three existing near/far implementations are
converted in the following commits.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: `PrecomputedVolumeField` delegates its geometry

**Files:**
- Modify: `$WT/src/shape/volume_field.jl:1-140` (docstring, struct, constructor), `:304-314` (`in_field_box`, `_field_scaled_targets`), `:324-387` (evaluators)
- Test: `$WT/test/shape/volume_field.jl:44-93`

**Interfaces:**
- Consumes: `near_field_geometry`, `in_near_region`, `NearFieldGeometry` from Task 2
- Produces:
  - `PrecomputedVolumeField(vs; tol, kmax = nothing, c_pad = 5.0, compute_pot = true, compute_grad = true, cache_fft = false, cache_fft_pad = 1.25)` — `margin_h` is **gone**, replaced by `c_pad`
  - `PrecomputedVolumeField` field layout: `lo`, `hi`, `center`, `dks` are **replaced** by a single `geom::NearFieldGeometry{T}`. `sources`, `charges`, `nmodes`, `kmax`, `tol`, `prefactor`, `coeff`, `grad_coeff`, `nfdim`, `pot_grid`, `grad_grid`, `interp_nspread`, `interp_hc` are unchanged.
  - `in_field_box(f, targets, i)` keeps its signature and now delegates to `in_near_region(f.geom, targets, i)`

- [ ] **Step 1: Write the failing test**

In `$WT/test/shape/volume_field.jl`, append a new testset at the end of the file:

```julia
@testset "PrecomputedVolumeField geometry matches near_field_geometry" begin
    gsrc = BoundaryIntegral.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 12, 1e-6)
    g = BoundaryIntegral.near_field_geometry(gsrc; c_pad = 5.0)
    f = PrecomputedVolumeField(gsrc; tol = 1e-6, c_pad = 5.0)
    @test f.geom.lo == g.lo
    @test f.geom.hi == g.hi
    @test f.geom.dks == g.dks
    @test f.geom.L == g.L
    # classification agrees between the struct helper and the free function
    tg = [0.0 5.0; 0.0 5.0; 0.0 5.0]
    @test BoundaryIntegral.in_field_box(f, tg, 1) == BoundaryIntegral.in_near_region(g, tg, 1)
    @test BoundaryIntegral.in_field_box(f, tg, 2) == BoundaryIntegral.in_near_region(g, tg, 2)
    # c_pad widens the near region
    f2 = PrecomputedVolumeField(gsrc; tol = 1e-6, c_pad = 10.0)
    @test f2.geom.hi[1] > f.geom.hi[1]
    # margin_h is gone
    @test_throws MethodError PrecomputedVolumeField(gsrc; tol = 1e-6, margin_h = 5.0)
end
```

Note: an unrecognised keyword raises `MethodError` in Julia, not `ArgumentError`.

- [ ] **Step 2: Run to verify it fails**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | grep -A5 "geometry matches near_field_geometry"
```

Expected: FAIL — `type PrecomputedVolumeField has no field geom`.

- [ ] **Step 3: Replace the struct fields**

In `$WT/src/shape/volume_field.jl`, in the struct at lines 58-81, delete the four lines

```julia
    lo::NTuple{3, T}
    hi::NTuple{3, T}
    center::NTuple{3, T}
    dks::NTuple{3, T}
```

and insert in their place:

```julia
    geom::NearFieldGeometry{T}
```

Keep `nmodes`, `kmax`, `tol`, `prefactor` and everything below unchanged.

- [ ] **Step 4: Rewrite the constructor geometry block**

In the constructor, change the keyword `margin_h::Real = 5.0` (line 87) to:

```julia
    c_pad::Real = 5.0,
```

Then replace lines 97-112 (from `sources, charges = ...` through the `nmodes = ...` line) with:

```julia
    sources, charges = _volume_source_fmm_sources(vs)
    geom = near_field_geometry(vs; c_pad = c_pad)
    # k_max keeps the min-nearest-neighbour spacing: Section 3.3 leaves k_Nyq open
    # for a general lattice, so only h_n is pinned to Eq. (3.3)'s ||A_rho||_2.
    km = isnothing(kmax) ? T(_estimate_tkm3dc_kmax(_estimate_source_spacing(vs))) : T(kmax)
    km > zero(T) || throw(ArgumentError("kmax must be positive"))
    center = geom.center
    dks = geom.dks
    kx = TKM3D.centered_mode_axis(dks[1], km)
    ky = TKM3D.centered_mode_axis(dks[2], km)
    kz = TKM3D.centered_mode_axis(dks[3], km)
    nmodes = (length(kx), length(ky), length(kz))
```

In the mode-filter loop (lines 119-124), change `TKM3D.truncated_laplace3d_hat(k, Lbig)` to `TKM3D.truncated_laplace3d_hat(k, geom.L)`.

In the `return PrecomputedVolumeField{T}(...)` call (lines 133-139), replace the `lo, hi, (center[1], center[2], center[3]), dks,` arguments with a single `geom,`, keeping the argument order aligned with the new field order.

- [ ] **Step 5: Update the three consumers of the removed fields**

`volume_field.jl:304-307`:

```julia
@inline in_field_box(f::PrecomputedVolumeField, targets::AbstractMatrix, i::Integer) =
    in_near_region(f.geom, targets, i)
```

`volume_field.jl:309-314`:

```julia
function _field_scaled_targets(f::PrecomputedVolumeField{T}, targets, idxs) where {T}
    dks = f.geom.dks
    c = f.geom.center
    txn = T[dks[1] * (targets[1, i] - c[1]) for i in idxs]
    tyn = T[dks[2] * (targets[2, i] - c[2]) for i in idxs]
    tzn = T[dks[3] * (targets[3, i] - c[3]) for i in idxs]
    return txn, tyn, tzn
end
```

Then grep for any other use of the removed fields and fix each:

```bash
cd ~/codes/bi-cpad-align
timeout 60 grep -rn "\.lo\b\|\.hi\b\|\.dks\|\.center" src/shape/volume_field.jl
timeout 60 grep -rn "margin_h" src/ test/
```

Both greps must come back clean of `PrecomputedVolumeField` field accesses and of `margin_h` before moving on.

- [ ] **Step 6: Update the constructor docstring**

The docstring at `volume_field.jl:1-57` documents `margin_h` at lines 2 and 11. Replace those mentions with `c_pad`, and state the geometry source:

```
by `c_pad * h` where `h = ||A_rho||_2` (Eq. 3.3); the Fourier box follows
Eqs. (3.11) and (3.16) via `near_field_geometry`
```

- [ ] **Step 7: Run the full suite**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | tail -60
```

Expected: the new geometry testset PASSES.

**Expect real changes in `test/shape/volume_field.jl:44-93`.** The old geometry ran at `η ≈ 1.06` on the paper's scale with `h_n = 4.5h`; the new one runs at `η = 1` with `h_n = 5h`. The mode count therefore changes and the `B_pad` corners move outward slightly.

Two things to check specifically:

1. `volume_field.jl` test line 53, `count(...) == 4`. The near box grew by `0.5h` per side. The far targets in that fixture are at `3.0`-`4.0` while `GaussianVolumeSource((0,0,0), 0.3, 12, 1e-6)` has support radius `sqrt(2*0.09*log(1e6)) ≈ 1.58`, so `h ≈ 2*1.58/12 ≈ 0.26` and `B_pad` reaches about `1.58 + 0.13 + 1.32 ≈ 3.03`. **This count may now be 5, not 4.** If it changes, that is correct behaviour, not a bug: fix the fixture by moving the nearest far target further out (e.g. `3.0 → 5.0` in all three coordinates of column 5) so the intent — 4 near, 4 far — is preserved. Do not change the assertion to `== 5`; the surrounding loops at lines 63-79 assume indices 1:4 are near and 5:8 are far.
2. The `rtol` values at lines 69-70 and 77-78. Mode count went down ~16%, so the in-box spectral error may rise slightly. If `rtol = 5e-4` or `1e-2` now fails, report the actual observed error before changing anything.

- [ ] **Step 8: Run the `cache_fft` testset explicitly**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | grep -B2 -A10 "cache_fft"
```

`cache_fft` derives `nfdim` from `nmodes`, which just changed, and it reaches eight TKM3D private functions (enumerated at `src/shape/volume_field.jl:151-163`). This testset must pass, not be skipped.

- [ ] **Step 9: Commit**

```bash
cd ~/codes/bi-cpad-align
git add src/shape/volume_field.jl test/shape/volume_field.jl
git commit -m "$(cat <<'EOF'
volume_field: use the paper's near/far geometry; margin_h -> c_pad

PrecomputedVolumeField took L to be the diagonal of B_pad and the period to
be (l + 2h_n) + L, as if sources could sit anywhere in B_pad rather than only
in B, and padded the sample bounding box so its effective h_n was (c_pad-1/2)h.
It now delegates to near_field_geometry, which follows Eqs. (3.9)/(3.11)/(3.16)
exactly. The tighter period cuts the mode array by roughly 16%.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: The KDTree ball becomes the Eq. (3.9) box

**Files:**
- Modify: `$WT/src/shape/box3d_fmm_helpers.jl:85-113`
- Modify: `$WT/src/solver/dielectric_box3d.jl:162-182`, `:253-255`
- Test: `$WT/test/core/source_geometry.jl` (append)

**Interfaces:**
- Consumes: `near_field_geometry`, `in_near_region` from Task 2
- Produces:
  - `_classify_near_far_targets(targets::Matrix{T}, vs::VolumeSource{T,3}; c_pad::Real = 5.0) -> Vector{Bool}` — note the change from positional `(targets, vs, h, h_factor)` to keyword `c_pad`, and that `h` is no longer passed in
  - `_classify_near_far_panels(panels, vs::VolumeSource{T,3}; c_pad::Real = 5.0) -> Vector{Bool}`
  - `rhs_dielectric_box3d_hybrid(interface, vs, eps_src, fmm_tol; tkm_kmax = nothing, c_pad = 5.0)` — `h_factor` renamed

- [ ] **Step 1: Write the failing test**

Append to `$WT/test/core/source_geometry.jl`:

```julia
@testset "_classify_near_far_targets uses the Eq. (3.9) box" begin
    gsrc = BI.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 12, 1e-6)
    g = BI.near_field_geometry(gsrc; c_pad = 5.0)

    # A corner of B_pad is inside the box but far from every source point, so the
    # old KDTree ball of radius 5h classified it FAR while Eq. (3.9) calls it NEAR.
    corner = [g.hi[1] - 1e-9; g.hi[2] - 1e-9; g.hi[3] - 1e-9]
    inside = [0.0; 0.0; 0.0]
    outside = [g.hi[1] + g.hn; 0.0; 0.0]
    targets = hcat(corner, inside, outside)

    is_near = BI._classify_near_far_targets(targets, gsrc; c_pad = 5.0)
    @test is_near == [true, true, false]

    # agrees with in_near_region on every column, by construction
    for i in 1:size(targets, 2)
        @test is_near[i] == BI.in_near_region(g, targets, i)
    end
end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | grep -A8 "Eq. (3.9) box"
```

Expected: FAIL — a `MethodError` on the `c_pad` keyword, since the current signature is positional `(targets, vs, h, h_factor)`.

- [ ] **Step 3: Replace both classifiers**

In `$WT/src/shape/box3d_fmm_helpers.jl`, replace `_classify_near_far_panels` (lines 85-103) and `_classify_near_far_targets` (lines 105-113) with box tests. Read the existing bodies first to see exactly what each returns (a `Vector{Bool}` over targets, and over panels keyed on panel centroids) and preserve that contract:

```julia
# Section 3 classification (Eq. 3.9): a target is near iff it lies in B_pad.
# This replaces a KDTree ball of radius c_pad*h, which disagreed with the paper
# near the corners of B_pad.
function _classify_near_far_targets(
    targets::Matrix{T},
    vs::VolumeSource{T, 3};
    c_pad::Real = 5.0,
) where {T}
    g = near_field_geometry(vs; c_pad = c_pad)
    n = size(targets, 2)
    return [in_near_region(g, targets, i) for i in 1:n]
end
```

For `_classify_near_far_panels`, keep whatever per-panel representative point the current implementation uses (read lines 85-103; it builds a point set from the panels and calls `inrange`). Convert it to build the same representative points into a `3 × n_panels` matrix and then apply `in_near_region`:

```julia
function _classify_near_far_panels(
    panels::Vector{TempPanel3D{T}},
    vs::VolumeSource{T, 3};
    c_pad::Real = 5.0,
) where {T}
    g = near_field_geometry(vs; c_pad = c_pad)
    pts = _panel_representative_points(panels)   # 3 × n, same points the old body used
    return [in_near_region(g, pts, i) for i in 1:size(pts, 2)]
end
```

If the old body inlined the representative-point construction rather than calling a helper, extract it into `_panel_representative_points(panels)` in the same file so both classifiers share it.

- [ ] **Step 4: Update the callers**

`$WT/src/solver/dielectric_box3d.jl:168` — rename the keyword:

```julia
    c_pad::Float64 = 5.0,
```

Line 173 becomes `c_pad > 0 || throw(ArgumentError("c_pad must be positive"))`, and line 182:

```julia
    is_near = _classify_near_far_targets(targets, vs; c_pad = c_pad)
```

The local `h = ...` that fed the old call is now unused at that site; check whether `h` is still needed for `tkm_kmax` a few lines below before deleting it.

`dielectric_box3d.jl:253-255` — rename the forwarded keyword from `h_factor = h_factor` to `c_pad = c_pad` and its parameter at line 253.

Then find every remaining caller:

```bash
cd ~/codes/bi-cpad-align
timeout 60 grep -rn "h_factor\|_classify_near_far" src/ test/
```

Update each hit. `src/shape/box3d_rhs_adaptive.jl` and `src/shape/box3d_multi.jl` were flagged in the spec as carrying `tkm_kmax` kwargs; check whether they also forward `h_factor`.

- [ ] **Step 5: Drop the now-unused import if nothing else needs it**

```bash
cd ~/codes/bi-cpad-align
timeout 60 grep -rn "KDTree\|inrange\|NearestNeighbors\|knn" src/
```

`_estimate_source_spacing` still uses `KDTree`/`knn`, so the dependency stays. Only remove an `inrange` import if it is now unreferenced.

- [ ] **Step 6: Run the full suite**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | tail -60
```

Expected: the new classification testset PASSES.

**Expect changes in `test/shape/box3d.jl:223, 248, 260, 289, 309, 343`**, which exercise `_rhs_volume_targets_hybrid`. The near set grew (box corners that the ball excluded are now near), so more targets go through the TKM branch. Since TKM is the accurate branch near the source, errors should improve or hold. If an assertion tightens past its tolerance in the *good* direction, that is fine; if any degrades, report the observed number before touching the tolerance.

- [ ] **Step 7: Commit**

```bash
cd ~/codes/bi-cpad-align
git add src/shape/box3d_fmm_helpers.jl src/solver/dielectric_box3d.jl src/shape/box3d_rhs_adaptive.jl src/shape/box3d_multi.jl test/core/source_geometry.jl
git commit -m "$(cat <<'EOF'
box3d: near/far classification is Eq. (3.9)'s box, not a KDTree ball

_classify_near_far_targets used inrange with radius c_pad*h around the source
points, which disagrees with the paper near the corners of B_pad. Both
classifiers now delegate to in_near_region. h_factor -> c_pad throughout.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: `evaluate_batch_potential` on the paper's geometry

**Files:**
- Modify: `$WT/src/solver/lattice_batch.jl:53-59` (`LatticeBatch`), `:139` (`assemble_lattice_batch` return), `:147-158` (`envelope_volume_source`), `:165-170` (`batch_volume_sources`), `:174-278` (`evaluate_batch_potential`)
- Modify: `$WT/src/campaign/toml_input.jl:39`, `:75`; `$WT/src/campaign/tasks.jl:306-321`
- Modify: `$WT/test/fixtures/system_overrides.toml:40`, `system_small.toml:39`, `system_smooth_lat.toml:39`
- Test: `$WT/test/solver/lattice_batch.jl:63-152`

**Interfaces:**
- Consumes: `near_field_geometry`, `in_near_region` (Task 2); `PrecomputedVolumeField(...; c_pad, compute_grad)`, `volume_field_potential` (Task 3); `with_density` (Task 1)
- Produces:
  - `LatticeBatch` gains a final field `lattice_basis::Union{Nothing, NTuple{3, NTuple{3, Float64}}}`
  - `evaluate_batch_potential(interface, Σ, sources, targets; lhs_tol, volume_tol, c_pad = 5.0, range_factor = 5.0, screen_boxes = nothing, screen_epses = nothing, screen_eps_out = 1.0)` — `far_pad` is **gone**
  - `CampaignInput` field `far_pad_steps` renamed to `c_pad`, default `5.0`, TOML key `[eval] c_pad`

- [ ] **Step 1: Write the failing test**

Append to `$WT/test/solver/lattice_batch.jl`, inside the same outer testset that defines `res_e`, `b_e`, `dg_e`:

```julia
    @testset "evaluate_batch_potential geometry and the c_pad 2-vs-5 bound" begin
        # LatticeBatch propagates the lattice basis, so h = ||A_rho||_2 is available.
        @test b_e.lattice_basis !== nothing
        for v in batch_volume_sources(b_e)
            @test v.lattice_basis == b_e.lattice_basis
        end
        @test envelope_volume_source(b_e).lattice_basis == b_e.lattice_basis

        # The batch split agrees with the shared helper on the screened source.
        targets = b_e.positions
        sa = BoundaryIntegral.screened_volume_source(res_e.interface, res_e.sources[1],
            BoundaryIntegral.SharpScreening())
        g = BoundaryIntegral.near_field_geometry(sa; c_pad = 5.0)
        @test all(BoundaryIntegral.in_near_region(g, targets, i) for i in 1:size(targets, 2))

        # Bound the un-rerun Section 6.4 change: how much does c_pad 2 -> 5 move Phi?
        far = hcat(([6.0 * cos(t), 6.0 * sin(t), 0.0] for t in range(0, 2π; length = 9)[1:8])...)
        mixed = hcat(targets, far)
        Φ2 = evaluate_batch_potential(res_e.interface, res_e.sigma, res_e.sources, mixed;
            lhs_tol = 1e-6, volume_tol = 1e-8, c_pad = 2.0)
        Φ5 = evaluate_batch_potential(res_e.interface, res_e.sigma, res_e.sources, mixed;
            lhs_tol = 1e-6, volume_tol = 1e-8, c_pad = 5.0)
        scale = maximum(abs.(Φ5))
        bound = maximum(abs.(Φ2 .- Φ5)) / scale
        @info "c_pad 2 vs 5: max relative difference in Phi = $bound"
        # Section 6.4 quotes U to ~0.1 eV out of ~2 eV and a symmetry residual of
        # 4.3e-3; 1e-3 keeps the c_pad change an order of magnitude below both.
        @test bound < 1e-3

        # far_pad is gone
        @test_throws MethodError evaluate_batch_potential(res_e.interface, res_e.sigma,
            res_e.sources, targets; lhs_tol = 1e-6, volume_tol = 1e-8, far_pad = 0.1)
    end
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | grep -A12 "c_pad 2-vs-5 bound"
```

Expected: FAIL — `type LatticeBatch has no field lattice_basis`, and a `MethodError` on the `c_pad` keyword.

- [ ] **Step 3: Propagate the basis through `LatticeBatch`**

`$WT/src/solver/lattice_batch.jl:53-59`, add the field:

```julia
struct LatticeBatch
    pair_ids::Vector{Tuple{Int,Int}}     # (orbital id i, orbital id j) per column
    gidx::Vector{NTuple{3,Int}}          # n shared global grid indices (sorted)
    positions::Matrix{Float64}           # 3 × n
    weights::Vector{Float64}             # n (uniform: |det cell| / (nx ny nz))
    densities::Matrix{Float64}           # n × K raw (unscreened) pair densities
    # Primitive sampling-cell basis of the virtual global grid, for Eq. (3.3)'s h.
    lattice_basis::Union{Nothing, NTuple{3, NTuple{3, Float64}}}
end
```

`assemble_lattice_batch` already has the cell vectors in scope at line 137. Change lines 137-139 to:

```julia
    At, Bt, Ct = true_cell_vectors(t1)
    w = abs(det(hcat(collect(At), collect(Bt), collect(Ct)))) / (nx * ny * nz)
    # primitive step vectors: the cell vectors divided by the grid counts
    lb = ((At[1] / nx, At[2] / nx, At[3] / nx),
          (Bt[1] / ny, Bt[2] / ny, Bt[3] / ny),
          (Ct[1] / nz, Ct[2] / nz, Ct[3] / nz))
    return LatticeBatch(copy(pairs), gk, positions, fill(w, m), dk, lb)
```

`envelope_volume_source` (line 157):

```julia
    return VolumeSource(copy(b.positions), copy(b.weights), env; lattice_basis = b.lattice_basis)
```

`batch_volume_sources` (lines 166-169):

```julia
    return VolumeSource{Float64, 3}[
        VolumeSource(copy(b.positions), copy(b.weights), b.densities[:, k];
                     lattice_basis = b.lattice_basis)
        for k in 1:num_pairs(b)
    ]
```

Then grep for any other `LatticeBatch(` construction that now needs the extra argument:

```bash
cd ~/codes/bi-cpad-align
timeout 60 grep -rn "LatticeBatch(" src/ test/
```

- [ ] **Step 4: Rewrite the near/far split and the near branch**

In `$WT/src/solver/lattice_batch.jl`, change the signature at lines 198-204: replace `far_pad::Float64,` with `c_pad::Float64 = 5.0,`.

Replace lines 230-253 (the split plus the `ltkm3dc` near loop) with:

```julia
    # incident part: Section 3 near/far split (Eq. 3.9) on the screened source.
    # screened_volume_source preserves positions, so the geometry is identical for
    # every column and is built once.
    geom = near_field_geometry(screened[1]; c_pad = c_pad)
    near_idx = findall(i -> in_near_region(geom, targets, i), 1:nt)
    far_idx = setdiff(1:nt, near_idx)

    # near targets: TKM via PrecomputedVolumeField, one column at a time.
    # ltkm3dc cannot be used here: it derives its own truncation radius from the
    # combined source+target box (TKM3D src/continuous.jl:133-146) and exposes no
    # hook for Eq. (3.11)'s L. Building one field per column and discarding it
    # keeps peak memory at a single coefficient array, matching the old call.
    if !isempty(near_idx)
        near_targets = targets[:, near_idx]
        for a in 1:K
            fld = PrecomputedVolumeField(screened[a];
                tol = volume_tol, c_pad = c_pad, compute_grad = false)
            Φ[near_idx, a] .+= volume_field_potential(fld, near_targets)
        end
    end
```

The far branch (lines 256-275) is unchanged.

- [ ] **Step 5: Update the docstring**

The docstring at lines 174-196 describes `far_pad` as an absolute length and says "`far_pad` ≳ 2 grid steps". Replace those two passages:

```
- `u_inc[rho_a]`: batch-level near/far split on the Section 3 padded near region
  B_pad (Eq. 3.9), `h_n = c_pad * ||A_rho||_2`. Near targets: TKM via
  PrecomputedVolumeField per source, on the paper's Fourier geometry
  (Eqs. 3.11, 3.16). Far targets: ONE `nd = K` point-charge FMM (potential/4pi)
  over the screened quadrature points.

`c_pad` is DIMENSIONLESS, a multiple of the lattice spacing; default 5.0.
```

- [ ] **Step 6: Update the campaign plumbing**

`$WT/src/campaign/toml_input.jl:39`, rename the struct field `far_pad_steps::Float64` to `c_pad::Float64`. Line 75, change the default and the TOML key:

```julia
        Float64(get(get(d, "eval", Dict()), "c_pad", 5.0)), toml_path)
```

`$WT/src/campaign/tasks.jl:306-321` — `max_step` is no longer needed:

```julia
function eval_batch_core(br::BatchResult, targets, store, dg, c::CampaignInput)
    K = length(br.pair_ids)
    pos = grid_positions(dg, br.gidx)
    At, Bt, Ct = true_cell_vectors(dg)
    lb = ((At[1] / dg.nx, At[2] / dg.nx, At[3] / dg.nx),
          (Bt[1] / dg.ny, Bt[2] / dg.ny, Bt[3] / dg.ny),
          (Ct[1] / dg.nz, Ct[2] / dg.nz, Ct[3] / dg.nz))
    sources = [VolumeSource(copy(pos), copy(br.weights), br.densities[:, k];
                            lattice_basis = lb) for k in 1:K]

    Φ = evaluate_batch_potential(br.interface, br.sigma, sources, targets.positions;
        lhs_tol   = c.solve["lhs_tol"],
        volume_tol = c.solve["volume_tol"],
        c_pad      = c.c_pad,
        screen_boxes = c.boxes, screen_epses = c.epses, screen_eps_out = c.eps_out)
```

Keep the rest of the function unchanged. Note `norm` may now be unused in this file — check before removing the import.

In the three fixtures `$WT/test/fixtures/system_overrides.toml`, `system_small.toml`, `system_smooth_lat.toml`, replace `far_pad_steps = 2.0` with `c_pad = 5.0`.

Then sweep for stragglers:

```bash
cd ~/codes/bi-cpad-align
timeout 60 grep -rn "far_pad" src/ test/ docs/
```

Must come back empty.

- [ ] **Step 7: Update the existing `lattice_batch` tests**

In `$WT/test/solver/lattice_batch.jl`, replace every `far_pad = ...` local and every `far_pad = far_pad` argument with `c_pad = 5.0`. Affected lines: 67, 70, 105, 109, 128, 144, 151, 158, 160.

The reference in those tests is `ltkm3dc` called directly, which now uses a *different* Fourier geometry than the code under test. Both are aliasing-free representations of the same integral, so they should agree at the Fourier-tail plus NUFFT level. Run and read the `@info` lines:

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | grep -A3 "max_rel_diff"
```

If `max_rel_diff < 1e-5` still holds, leave the tolerances alone. If it does not, the coarse `b_e` fixture is under-resolved (its own comment at line 64-65 says so) and the two geometries need not agree tightly on it — in that case switch the reference from `ltkm3dc` to `PrecomputedVolumeField` on the same source and keep `1e-5`, rather than loosening the number. Record which you did and why in the commit message.

- [ ] **Step 8: Run the full suite**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | tail -60
```

Expected: all testsets pass, including the new `c_pad 2-vs-5 bound`. Note the reported bound — it is the evidence for skipping the Section 6.4 rerun.

If `bound >= 1e-3`, do NOT loosen the assertion. Stop and report: that is the signal that the no-rerun decision needs revisiting.

- [ ] **Step 9: Commit**

```bash
cd ~/codes/bi-cpad-align
git add src/solver/lattice_batch.jl src/campaign/toml_input.jl src/campaign/tasks.jl test/fixtures/ test/solver/lattice_batch.jl
git commit -m "$(cat <<'EOF'
lattice_batch: Section 3 near/far split at c_pad = 5

evaluate_batch_potential padded the sample bbox by an absolute far_pad
(campaign default 2 grid steps) and let ltkm3dc pick its own truncation radius
from the combined source+target box, so L depended on batch composition. It now
uses near_field_geometry and routes the near branch through
PrecomputedVolumeField, one column at a time so peak memory is unchanged.

far_pad_steps -> c_pad (dimensionless, default 5.0) through CampaignInput and
the fixtures. LatticeBatch carries the primitive lattice basis so Eq. (3.3)'s
h = ||A_rho||_2 survives batch_volume_sources and envelope_volume_source.

Adds a test bounding the max relative change in Phi between c_pad 2 and 5,
which is the evidence for not regenerating the Section 6.4 campaign.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Cross-path equivalence guard

**Files:**
- Test: `$WT/test/core/source_geometry.jl` (append)

**Interfaces:**
- Consumes: everything from Tasks 1-5.
- Produces: no new API. This task exists so the three implementations cannot silently drift apart again.

- [ ] **Step 1: Write the test**

Append to `$WT/test/core/source_geometry.jl`:

```julia
@testset "all near/far paths agree on classification" begin
    # One cubic and one skew lattice; every path must classify identically.
    function cubic_source(n)
        h = 2.0 / n
        xs = collect(-1.0 + h/2 .+ h .* (0:n-1))
        dens = Array{Float64,3}(undef, n, n, n)
        for k in 1:n, j in 1:n, i in 1:n
            r2 = xs[i]^2 + xs[j]^2 + xs[k]^2
            dens[i,j,k] = exp(-r2 / (2 * 0.25^2))
        end
        return VolumeSource((xs, xs, xs), fill(h^3, n, n, n), dens)
    end

    function skew_source(n)
        frac = collect((i - 1) / n for i in 1:n)
        At = (2.0, 0.0, 0.0); Bt = (0.7, 1.9, 0.0); Ct = (0.0, 0.0, 2.2)
        dens = fill(1.0, n, n, n)
        jac = abs(det([2.0 0.7 0.0; 0.0 1.9 0.0; 0.0 0.0 2.2]))
        return VolumeSource((frac, frac, frac), fill(jac / n^3, n, n, n), dens,
                            (0.0, 0.0, 0.0), (At, Bt, Ct))
    end

    for vs in (cubic_source(12), skew_source(12))
        g = BI.near_field_geometry(vs; c_pad = 5.0)
        # a spread of targets straddling the B_pad boundary in every direction
        pts = Float64[]
        for sx in (-1.0, 0.0, 1.0), sy in (-1.0, 0.0, 1.0), sz in (-1.0, 0.0, 1.0)
            for f in (0.5, 0.999, 1.001, 2.0)
                append!(pts, [g.center[1] + sx * f * (g.hi[1] - g.center[1]),
                              g.center[2] + sy * f * (g.hi[2] - g.center[2]),
                              g.center[3] + sz * f * (g.hi[3] - g.center[3])])
            end
        end
        targets = reshape(pts, 3, :)

        ref = [BI.in_near_region(g, targets, i) for i in 1:size(targets, 2)]

        # path: box3d classifier
        @test BI._classify_near_far_targets(targets, vs; c_pad = 5.0) == ref

        # path: PrecomputedVolumeField
        f = PrecomputedVolumeField(vs; tol = 1e-6, compute_grad = false, c_pad = 5.0)
        @test [BI.in_field_box(f, targets, i) for i in 1:size(targets, 2)] == ref
        # and its stored geometry is the same object's worth of numbers
        @test f.geom.lo == g.lo && f.geom.hi == g.hi
        @test f.geom.L == g.L && f.geom.dks == g.dks
    end
end

@testset "cubic lattice: ||A_rho||_2 equals the old spacing estimate" begin
    # Guarantees that results on cubic grids are unchanged by the h redefinition.
    for n in (8, 13, 24)
        h = 2.0 / n
        xs = collect(-1.0 + h/2 .+ h .* (0:n-1))
        vs = VolumeSource((xs, xs, xs), fill(h^3, n, n, n), fill(1.0, n, n, n))
        @test isapprox(BI.lattice_spacing(vs), BI._estimate_source_spacing(vs); rtol = 1e-12)
    end
end
```

- [ ] **Step 2: Run it**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | grep -A10 "all near/far paths agree"
```

Expected: PASS. If the box3d classifier and `in_field_box` disagree, one of Tasks 3-4 did not fully delegate — find it rather than adjusting the test.

The `f = 0.999` / `1.001` points sit just inside and just outside `∂B_pad` and are the ones that catch an off-by-half-a-cell error in `source_box`.

- [ ] **Step 3: Run the whole suite one final time and record the timing**

```bash
cd ~/codes/bi-cpad-align
time julia --project=. test/runtests.jl 2>&1 | tail -30
```

Compare against the Task 0 baseline. A large slowdown in the `lattice_batch` testsets would indicate the per-column `PrecomputedVolumeField` construction is more expensive than the old single `ltkm3dc` call; note the number either way.

- [ ] **Step 4: Commit**

```bash
cd ~/codes/bi-cpad-align
git add test/core/source_geometry.jl
git commit -m "$(cat <<'EOF'
test: guard against the three near/far paths drifting apart

Asserts that _classify_near_far_targets, PrecomputedVolumeField.in_field_box,
and in_near_region classify identically on cubic and skew lattices, including
targets just inside and just outside dB_pad. Also pins ||A_rho||_2 equal to the
old min-nearest-neighbour spacing on cubic grids, which is what makes cubic
results provably unchanged.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 5: Report Stage A results before starting Stage B**

Report: the `c_pad` 2-vs-5 bound, any tolerance that changed and why, the near-count fixture change in `test/shape/volume_field.jl`, and the test-suite timing delta. Stage B consumes this code, so it must be green first.

---

## Stage B — Figure 5 and the paper

### Task 7: Rewrite `fig5_data.jl`

**Files:**
- Modify: `~/work/four_index_integral_solver/codes/fig_gen/fig5_data.jl`

**Interfaces:**
- Consumes: `PrecomputedVolumeField(vs; tol, c_pad, compute_grad)`, `volume_field_potential`, `near_field_geometry`, `in_near_region`, `VolumeSource` (Tasks 1-3)
- Produces: `fig_gen/fig5_data.jls` with the same field names as before plus `c_pad`, `n_near`, `n_far`, `hn_by_n`. `eta_list` now spans `[0.5, 1.6]` with 20 points and is plotted in full (no slice).

- [ ] **Step 1: Point `fig_gen` at the worktree and resolve**

`fig_gen/Project.toml` already lists `BoundaryIntegral`. Its Manifest path-devs the live checkout, so switch it to the worktree for this work:

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
julia --project=. -e 'using Pkg; Pkg.develop(path = "/mnt/home/xgao1/codes/bi-cpad-align"); Pkg.resolve(); Pkg.precompile()' 2>&1 | tail -20
```

A "does not have X in its dependencies" precompile error here means a stale consumer Manifest; `Pkg.resolve()` in this project is the fix.

- [ ] **Step 2: Verify the geometry the figure will use, before touching the plot**

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
julia --project=. -e '
using BoundaryIntegral
const BI = BoundaryIntegral
for n in (8, 64)
    h = 2.0/n
    xs = collect(-1.0 + h/2 .+ h .* (0:n-1))
    vs = VolumeSource((xs, xs, xs), fill(h^3, n, n, n), fill(1.0, n, n, n))
    g = BI.near_field_geometry(vs; c_pad = 5.0)
    println("n=$n  h=$h  hn=$(g.hn)  l=$(g.l[1])  Bpad_hi=$(g.hi[1])  L=$(g.L)")
end'
```

Expected: `n=8 h=0.25 hn=1.25 l=2.0 Bpad_hi=2.25`, and `n=64 h=0.03125 hn=0.15625 l=2.0 Bpad_hi=1.15625`. These are the numbers the target set is designed around. If `l` is not exactly `2.0`, `source_box` is not producing the cell box and Task 2 needs revisiting.

- [ ] **Step 3: Rewrite the script**

Replace `fig_gen/fig5_data.jl` with:

```julia
#=
Data generation for Figure 5 (TKM/FMM hybrid incident-potential validation).

We evaluate
    u(x) = (1 / 4 pi) int rho(y) / |x - y| dy
for a normalized 3D Gaussian rho on B = [-1, 1]^3, sampled on a cell-centered
uniform Cartesian grid, and compare two evaluators at the SAME targets:
  - the Section 3 hybrid: TKM inside B_pad, FMM outside (PrecomputedVolumeField)
  - the pure particle sum of Eq. (3.13), at every target

The padding is h_n = c_pad * h with c_pad = 5, so both panels exercise the
near/far switch. Targets span both regions and keep their classification for
every n in the sweep:
  near: 200 points in [-0.9, 0.9]^3   -- inside B_pad even at n = 64 (B_pad ~ 1.156)
  far:  200 points with max|coord| in [2.5, 4] -- outside B_pad even at n = 8 (B_pad ~ 2.25)

The analytic reference is u_ref(x) = (1 / 4 pi) erf(|x| / (sqrt(2) s)) / |x|,
with limit (1 / 4 pi) sqrt(2 / pi) / s at x = 0, valid at every target.

Panels:
  (a) relative L2 error vs grid resolution n, for several tolerances tau
  (b) relative L2 error vs eta = min_alpha L_alpha / (l_alpha + h_n + L), n fixed.
      Uses a standalone TKM+FMM evaluator with a tunable dk, because the
      production path now sits exactly at eta = 1 and cannot be driven below it.
=#

using LinearAlgebra
using SpecialFunctions
using Serialization
using Printf
using Random
using FMM3D
using TKM3D
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const s_gauss = 0.10
const y0      = (0.0, 0.0, 0.0)

const box_half = 1.0
const l_box    = 2 * box_half
const C_PAD    = 5.0

const eps_list = (1e-3, 1e-6, 1e-9, 1e-12)
const n_list   = collect(8:8:64)

const n_for_eta = 64
const eta_list  = collect(10 .^ range(log10(0.5), log10(1.6); length = 20))

# ---------------------------------------------------------------------------
# Fixed target set: near + far, classification independent of n
# ---------------------------------------------------------------------------
const N_near = 200
const N_far  = 200

function build_targets()
    Random.seed!(42)
    t = Matrix{Float64}(undef, 3, N_near + N_far)
    # near: uniform in [-0.9, 0.9]^3, inside B_pad for every n (B_pad >= 1.156)
    for k in 1:N_near
        t[1, k] = 1.8 * (rand() - 0.5)
        t[2, k] = 1.8 * (rand() - 0.5)
        t[3, k] = 1.8 * (rand() - 0.5)
    end
    # far: max|coord| in [2.5, 4], outside B_pad for every n (B_pad <= 2.25)
    k = N_near
    while k < N_near + N_far
        p = 8.0 .* (rand(3) .- 0.5)          # uniform in [-4, 4]^3
        m = maximum(abs, p)
        (m >= 2.5 && m <= 4.0) || continue
        k += 1
        t[1, k] = p[1]; t[2, k] = p[2]; t[3, k] = p[3]
    end
    return t
end

const targets  = build_targets()
const N_targets = size(targets, 2)
const near_rng = 1:N_near
const far_rng  = (N_near + 1):(N_near + N_far)

# ---------------------------------------------------------------------------
# Density, source construction, analytic reference
# ---------------------------------------------------------------------------
@inline function rho(y::NTuple{3,Float64})
    r2 = (y[1]-y0[1])^2 + (y[2]-y0[2])^2 + (y[3]-y0[3])^2
    return exp(-r2 / (2 * s_gauss^2)) / (2π * s_gauss^2)^(3/2)
end

@inline function u_ref(x::NTuple{3,Float64})
    r = sqrt((x[1]-y0[1])^2 + (x[2]-y0[2])^2 + (x[3]-y0[3])^2)
    r < 1e-14 && return sqrt(2 / π) / s_gauss / (4π)
    return erf(r / (sqrt(2) * s_gauss)) / r / (4π)
end

tail_mass(R, s) = erfc(R / (sqrt(2)*s)) + sqrt(2/π) * (R/s) * exp(-R^2 / (2*s^2))

const u_r_t    = [u_ref((targets[1,k], targets[2,k], targets[3,k])) for k in 1:N_targets]
const u_r_norm = norm(u_r_t)

relerr(u) = norm(u .- u_r_t) / u_r_norm
relerr(u, rng) = norm(u[rng] .- u_r_t[rng]) / norm(u_r_t[rng])

# Cell-centered grid on B = [-1,1]^3. The identity-basis grid constructor stores
# A_rho = h * I, so BI.source_box returns exactly [-1,1]^3 and BI.lattice_spacing
# returns exactly h.
function grid_source(n::Int)
    h  = l_box / n
    xs = collect(-box_half + h/2 .+ h .* (0:n-1))
    weights = fill(h^3, n, n, n)
    density = Array{Float64,3}(undef, n, n, n)
    for k in 1:n, j in 1:n, i in 1:n
        density[i,j,k] = rho((xs[i], xs[j], xs[k]))
    end
    return VolumeSource((xs, xs, xs), weights, density)
end

# ---------------------------------------------------------------------------
# Standalone hybrid with a tunable Fourier spacing (panel b).
# Mirrors BI.near_field_geometry, then scales dk by 1/eta so eta < 1 is reachable.
# ---------------------------------------------------------------------------
function hybrid_eval_eta(vs::VolumeSource{Float64,3}, trg::Matrix{Float64},
                        eta::Float64; eps::Float64)
    g   = BI.near_field_geometry(vs; c_pad = C_PAD)
    src, q = BI._volume_source_fmm_sources(vs)
    km  = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(vs))
    dk  = ntuple(d -> 2π / (g.l[d] + g.hn + g.L) / eta, 3)

    kx = TKM3D.centered_mode_axis(dk[1], km)
    ky = TKM3D.centered_mode_axis(dk[2], km)
    kz = TKM3D.centered_mode_axis(dk[3], km)

    out  = Vector{Float64}(undef, size(trg, 2))
    inb  = [BI.in_near_region(g, trg, i) for i in 1:size(trg, 2)]
    bidx = findall(inb); oidx = findall(!, inb)

    if !isempty(bidx)
        srcx = dk[1] .* (vec(view(src,1,:)) .- g.center[1])
        srcy = dk[2] .* (vec(view(src,2,:)) .- g.center[2])
        srcz = dk[3] .* (vec(view(src,3,:)) .- g.center[3])
        coeff0 = TKM3D.FINUFFT.nufft3d1(srcx, srcy, srcz, complex.(q), -1, eps,
                                        length(kx), length(ky), length(kz))
        coeff = ndims(coeff0) == 4 ? dropdims(coeff0; dims = 4) : coeff0
        @inbounds for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
            k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2)
            coeff[ix,iy,iz] = k <= km ?
                coeff[ix,iy,iz] * TKM3D.truncated_laplace3d_hat(k, g.L) :
                zero(eltype(coeff))
        end
        txn = Float64[dk[1] * (trg[1,i] - g.center[1]) for i in bidx]
        tyn = Float64[dk[2] * (trg[2,i] - g.center[2]) for i in bidx]
        tzn = Float64[dk[3] * (trg[3,i] - g.center[3]) for i in bidx]
        vals = TKM3D._finufft_type2_eval_3d(txn, tyn, tzn, 1, eps, coeff)
        pref = dk[1] * dk[2] * dk[3] / (2π)^3
        for (m, i) in enumerate(bidx)
            out[i] = pref * real(vals[m])
        end
    end
    if !isempty(oidx)
        vals = lfmm3d(eps, src; charges = q, targets = trg[:, oidx], pgt = 1)
        for (m, i) in enumerate(oidx)
            out[i] = vals.pottarg[m] / (4π)
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Panel (a): convergence vs n — hybrid and pure particle sum, several tolerances
# ---------------------------------------------------------------------------
err_tkm_n  = Dict{Float64,Vector{Float64}}()
err_fmm_n  = Dict{Float64,Vector{Float64}}()
err_near_n = Dict{Float64,Vector{Float64}}()
err_far_n  = Dict{Float64,Vector{Float64}}()
for eps in eps_list
    err_tkm_n[eps]  = Float64[]
    err_fmm_n[eps]  = Float64[]
    err_near_n[eps] = Float64[]
    err_far_n[eps]  = Float64[]
end
hn_by_n = Float64[]

@info "Panel (a) — convergence vs n, c_pad = $C_PAD, $N_near near + $N_far far targets"
for n in n_list
    vs = grid_source(n)
    g  = BI.near_field_geometry(vs; c_pad = C_PAD)
    push!(hn_by_n, g.hn)

    # the target set must classify as designed at every n
    nn = count(i -> BI.in_near_region(g, targets, i), 1:N_targets)
    nn == N_near || error("n = $n: $nn targets classified near, expected $N_near " *
                          "(B_pad upper corner = $(g.hi[1]))")

    src, q = BI._volume_source_fmm_sources(vs)
    for eps in eps_list
        field = PrecomputedVolumeField(vs; tol = eps, c_pad = C_PAD, compute_grad = false)
        u_h   = volume_field_potential(field, targets)
        push!(err_tkm_n[eps],  relerr(u_h))
        push!(err_near_n[eps], relerr(u_h, near_rng))
        push!(err_far_n[eps],  relerr(u_h, far_rng))

        # pure particle sum at every target (Eq. 3.13); tau does not enter it,
        # so the four curves coincide
        u_p = lfmm3d(eps, src; charges = q, targets = targets, pgt = 1).pottarg ./ (4π)
        push!(err_fmm_n[eps], relerr(u_p))

        @info @sprintf("  n = %3d  tau = %.0e  E_hyb = %.3e (near %.3e, far %.3e)  E_part = %.3e",
                       n, eps, err_tkm_n[eps][end], err_near_n[eps][end],
                       err_far_n[eps][end], err_fmm_n[eps][end])
    end
end

# ---------------------------------------------------------------------------
# Panel (b): periodization threshold — vary eta at fixed n
# ---------------------------------------------------------------------------
err_tkm_eta = Dict{Float64,Vector{Float64}}()
for eps in eps_list
    err_tkm_eta[eps] = Float64[]
end

@info "Panel (b) — eta sweep at n = $n_for_eta"
vs_c = grid_source(n_for_eta)

# cross-check: at eta = 1 the standalone evaluator must reproduce the production path
let g = BI.near_field_geometry(vs_c; c_pad = C_PAD)
    f  = PrecomputedVolumeField(vs_c; tol = 1e-12, c_pad = C_PAD, compute_grad = false)
    u_prod = volume_field_potential(f, targets)
    u_std  = hybrid_eval_eta(vs_c, targets, 1.0; eps = 1e-12)
    d = maximum(abs.(u_prod .- u_std)) / maximum(abs.(u_prod))
    @info @sprintf("eta = 1 cross-check vs PrecomputedVolumeField: max rel diff = %.3e", d)
    d < 1e-10 || error("standalone eta evaluator disagrees with the production path ($d)")
end

for eps in eps_list
    for eta in eta_list
        u = hybrid_eval_eta(vs_c, targets, eta; eps = eps)
        push!(err_tkm_eta[eps], relerr(u))
    end
    @info @sprintf("  tau = %.0e :  E(0.5) = %.3e   E(~1) = %.3e   E(1.6) = %.3e",
                   eps, err_tkm_eta[eps][1],
                   err_tkm_eta[eps][argmin(abs.(eta_list .- 1.0))],
                   err_tkm_eta[eps][end])
end

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    s_gauss     = s_gauss,
    y0          = y0,
    box_half    = box_half,
    l_box       = l_box,
    c_pad       = C_PAD,
    hn_by_n     = hn_by_n,
    eps_list    = collect(eps_list),
    n_list      = n_list,
    err_tkm_n   = err_tkm_n,
    err_fmm_n   = err_fmm_n,
    err_near_n  = err_near_n,
    err_far_n   = err_far_n,
    n_for_eta   = n_for_eta,
    eta_list    = eta_list,
    err_tkm_eta = err_tkm_eta,
    N_targets   = N_targets,
    n_near      = N_near,
    n_far       = N_far,
    tail_mass   = tail_mass(box_half, s_gauss),
)

datapath = joinpath(@__DIR__, "fig5_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes=stat(datapath).size
```

`L_diag` and the old `grid_sources`/`tkm_eval_eta` are gone; `FINUFFT` is no longer imported directly (reached via `TKM3D.FINUFFT`).

- [ ] **Step 4: Smoke run before the full sweep**

Back up the existing data first, then run a reduced sweep by temporarily editing the two list constants:

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
cp fig5_data.jls fig5_data.jls.pre_cpad
sed -i.smokebak \
  -e 's/^const n_list   = collect(8:8:64)$/const n_list   = collect(8:8:16)/' \
  -e 's/^const eps_list = (1e-3, 1e-6, 1e-9, 1e-12)$/const eps_list = (1e-6,)/' \
  fig5_data.jl
julia --project=. fig5_data.jl 2>&1 | tail -30
```

Expected: the classification check passes at both `n = 8` and `n = 16` (no `error(...)`), the `eta = 1` cross-check reports `< 1e-10`, and the hybrid error is well below the particle-sum error.

The classification check is the important one. If it errors, the far shell and `B_pad` overlap at that `n` and the target-set bounds need adjusting — do not delete the check.

Restore the full sweep:

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
mv fig5_data.jl.smokebak fig5_data.jl
```

- [ ] **Step 5: Commit the script**

```bash
cd ~/work/four_index_integral_solver/codes
git add fig_gen/fig5_data.jl
git commit -m "$(cat <<'EOF'
fig5: hybrid evaluator at c_pad = 5 with mixed near+far targets

The experiment set h_n = 0 and evaluated only targets inside B, so it never
exercised the near/far switch that Section 3 specifies. It now uses
PrecomputedVolumeField at c_pad = 5 over 200 near plus 200 far targets, with
an assertion that the classification is the intended one at every n.

Panel (b)'s eta sweep keeps a standalone tunable-dk evaluator, since the
production path now sits exactly at eta = 1, and cross-checks against it there.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Regenerate Figure 5

**Files:**
- Modify: `~/work/four_index_integral_solver/codes/fig_gen/fig5_plot.jl:72-100`
- Regenerate: `fig_gen/figs/fig5_tkm_validation.pdf`

**Interfaces:**
- Consumes: `fig5_data.jls` from Task 7 (fields `eta_list`, `err_tkm_eta`, `err_tkm_n`, `err_fmm_n`, `n_list`, `eps_list`)
- Produces: the updated PDF at `fig_gen/figs/fig5_tkm_validation.pdf`

- [ ] **Step 1: Run the full data sweep**

`n = 64` with 4 tolerances plus a 20-point `η` sweep is the expensive part; `η = 1.6` needs roughly `(1.6)³ ≈ 4×` the modes of `η = 1`. Run it detached with a log:

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
nohup julia --project=. fig5_data.jl > fig5_data.log 2>&1 &
echo $!
```

Poll the log (do not block on it):

```bash
tail -20 ~/work/four_index_integral_solver/codes/fig_gen/fig5_data.log
```

Expected in the log: the classification check silent at all 8 values of `n`, the `eta = 1` cross-check `< 1e-10`, `E_hyb` dropping to the `τ` plateau, `E_part` stalling near `10⁻⁴`, and finally `Saved data`.

- [ ] **Step 2: Fix the panel (b) axis label and drop the hardcoded slice**

In `fig_gen/fig5_plot.jl`, line 79, change the x-label to include `h_n`:

```julia
                xlabel = L"\eta = \min_{\alpha}\, L_\alpha / (l_\alpha + h_n + L)",
```

Lines 86-91: `eta_list` now spans `[0.5, 1.6]` with 20 points, so plot all of it instead of slicing:

```julia
    for eps in eps_list
        scatterlines!(ax_b, eta_list, clip.(err_tkm_eta[eps]);
                      color = eps_colors[eps], marker = eps_markers[eps],
                      markersize = MS, linewidth = LW_DATA,
                      label = eps_label(eps))
    end
```

Line 99, widen the x limits to the new range:

```julia
    xlims!(ax_b, 0.48, 1.62)
```

Also update the panel (a) `ylabel`: it reads `\mathcal{E}_u`, which is still correct, so leave it.

- [ ] **Step 3: Render and inspect**

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
julia --project=. fig5_plot.jl 2>&1 | tail -10
```

Then look at the result:

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
julia --project=. -e '
using Serialization
d = open(deserialize, "fig5_data.jls", "r")
println("c_pad = ", d.c_pad, "  hn(n=8) = ", d.hn_by_n[1], "  hn(n=64) = ", d.hn_by_n[end])
for eps in d.eps_list
    println("tau = ", eps, "  E(n=64) = ", d.err_tkm_n[eps][end],
            "  near = ", d.err_near_n[eps][end], "  far = ", d.err_far_n[eps][end])
end
e = d.err_tkm_eta[d.eps_list[end]]
i = findfirst(j -> e[j] < 1e-6, eachindex(e))
println("tightest-tau eta crossing below 1e-6 at eta = ", d.eta_list[i])'
```

Read the PDF with the Read tool to confirm both panels render sensibly.

Expected: the `η` crossing lands somewhat **below** 1 — near `0.91` — because the test Gaussian's `s = 0.10` gives it a numerical support radius of about `0.6` rather than the full `B`, so the periodic image clears the target region early. Eq. (3.16) is a sufficient condition, not a sharp one. **Do not tune anything to move the crossing to 1.** Record the observed value; Task 9 words the caption to match.

- [ ] **Step 4: Commit**

```bash
cd ~/work/four_index_integral_solver/codes
git add fig_gen/fig5_plot.jl fig_gen/fig5_data.jls figs/fig5_tkm_validation.pdf fig_gen/figs/fig5_tkm_validation.pdf
git commit -m "$(cat <<'EOF'
fig5: regenerate at c_pad = 5 with the near+far target set

Panel (b)'s x-axis is now eta = min L_alpha/(l_alpha + h_n + L); the whole
sweep is plotted rather than a hardcoded slice.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

Adjust the paths in `git add` to whichever of the two `figs/` directories the plot script actually writes to (`fig5_plot.jl:109` writes `joinpath(@__DIR__, "figs/...")`), and copy it to `~/Articles/four_indices_bie/figs/` if that is a separate file rather than a symlink — check with `ls -l ~/Articles/four_indices_bie/figs/fig5_tkm_validation.pdf`.

---

### Task 9: Update Section 3.4 and the caption

**Files:**
- Modify: `~/Articles/four_indices_bie/main.tex:1227-1343`

**Interfaces:**
- Consumes: the observed numbers from Task 8 (the `η` crossing, the `n` at which each `τ` plateaus, the particle-sum plateau).

- [ ] **Step 1: Rewrite the subsection opening**

Replace lines 1229-1233 (the "We validate two components separately..." paragraph, whose last sentence disclaims the switch) with:

```latex
We validate three things: convergence of the hybrid evaluator, the
Fourier-spacing condition of Lemma~\ref{lem:padding}, and the near--far switch
itself. All experiments use $c_{\mathrm{pad}}=5$, so $h_n=5h$ throughout, and
every error is measured on a target set that straddles
$\partial\mathcal B_{\mathrm{pad}}$.
```

- [ ] **Step 2: Update the target-set description**

Replace lines 1274-1275 ("We measure the relative $L^2$ error at a fixed set of $200$ off-grid random targets inside $\mathcal B$.") with:

```latex
We measure the relative $L^2$ error on a fixed set of $400$ off-grid targets:
$200$ drawn uniformly from $[-0.9,0.9]^3$, which lie in
$\mathcal B_{\mathrm{pad}}$ for every resolution in the sweep, and $200$ with
$\max_\alpha\abs{x_\alpha}\in[2.5,4]$, which lie outside
$\mathcal B_{\mathrm{pad}}$ for every resolution. Because $h_n=5h$ shrinks with
$n$, holding both groups clear of the boundary is what makes the classification,
and hence the curves, independent of $n$: the near targets are evaluated by the
TKM and the far targets by the particle/FMM representation at every $n$.
```

- [ ] **Step 3: Rewrite the TKM convergence paragraph**

Replace lines 1279-1285 with:

```latex
The first experiment fixes $c_{\mathrm{pad}}=5$ and evaluates the complete
hybrid \eqref{eq:near_far_switch} at all $400$ targets.
Figure~\ref{fig:tkm}(a) shows rapid convergence as the source-grid resolution
$n$ increases, followed by a plateau set by the selected Fourier-tail and NUFFT
tolerances $\tau\in\{10^{-3},10^{-6},10^{-9},10^{-12}\}$. Both branches
converge spectrally: the near branch by the truncated-kernel argument above,
and the far branch because the integrand $\rho(y)G(x-y)$ is smooth at a far
target, so the uniform-lattice rule inherits the exponential accuracy of the
trapezoidal rule for a rapidly decaying integrand.
```

Insert the observed plateau resolution from Task 8 in place of the old "At $\tau=10^{-12}$ the error reaches machine precision by approximately $n=48$" sentence, using the measured value.

Delete lines 1287-1288 (the `\jk`/`\xg` pair about how `\tau` enters), since its `\xg` reply explicitly says "because the test uses only near targets with $h_n = 0$" and is now false.

- [ ] **Step 4: Fix the particle-sum comparison paragraph**

In lines 1290-1297, replace "at the same near targets" with "at the same $400$ targets", and replace the final sentence (lines 1295-1297, "This comparison illustrates why the particle representation is not used for near targets; it does not test the far-target FMM acceleration.") with:

```latex
Its error is dominated by the near targets, where the kernel is singular or
nearly singular on the scale of a source-lattice cell; at the far targets it
agrees with the hybrid, which is what the switch is designed to exploit.
```

- [ ] **Step 5: Rewrite the periodization paragraph**

Replace lines 1307-1318 with:

```latex
The second experiment fixes $n=64$ and $c_{\mathrm{pad}}=5$, and varies the
Fourier spacings through
\begin{equation}
\label{eq:padding_ratio}
    \eta
    =
    \min_\alpha
    \frac{L_\alpha}{l_\alpha+h_n+L}.
\end{equation}
```

Then reword lines 1319-1324 so the claim matches the data. Using the crossing value measured in Task 8:

```latex
For $\eta$ well below $1$, periodic images overlap the target region and
introduce $\mathcal O(1)$ contamination; once $\eta\geq1$ the periodization error
is absent and the total error drops to the level set by the remaining
discretization parameters. Equation~\eqref{eq:padding} is sufficient rather than
sharp: the observed transition lies slightly below $\eta=1$ because the
numerical support of the test Gaussian is a ball of radius $\approx0.6$ rather
than all of $\mathcal B$, so the nearest periodic image clears the target region
a little earlier than the bound requires. Figure~\ref{fig:tkm}(b) therefore
confirms the sufficiency of \eqref{eq:padding} and the plateau beyond it.
```

Update the `\xg` reply at lines 1326-1327 to match, since it currently asserts that $\eta=1$ is "the transition threshold".

- [ ] **Step 6: Rewrite the caption**

Replace the caption body at lines 1332-1342 with:

```latex
    \caption{Component-level validation of the hybrid incident-potential
    evaluator for the Gaussian density \eqref{eq:gaussian_source}, at
    $c_{\mathrm{pad}}=5$. Errors are measured on $400$ targets: $200$ in
    $[-0.9,0.9]^3$, inside $\mathcal B_{\mathrm{pad}}$ at every resolution, and
    $200$ with $\max_\alpha\abs{x_\alpha}\in[2.5,4]$, outside it at every
    resolution.
    (a) Relative $L^2$ error versus source-grid resolution $n$ at
    tolerances $\tau\in\{10^{-3},10^{-6},10^{-9},10^{-12}\}$, compared
    with the particle midpoint sum \eqref{eq:fmm_far_field} evaluated at all
    $400$ targets. (b) Relative $L^2$ error versus the period ratio
    $\eta=\min_\alpha L_\alpha/(l_\alpha+h_n+L)$ at $n=64$. Beyond the
    transition the remaining discretization errors set the plateau.}
```

- [ ] **Step 7: Check the Notation table and the error-budget text**

Line 154 defines `$h_n,\ c_{\mathrm{pad}}$`; line 1084 attributes $\mathcal E_{\mathrm{far}}$ to $c_{\mathrm{pad}}$. Both remain correct — read them and confirm, but do not change §3.1-3.3: the code moved to the paper, not the reverse.

Line 853 says `$c_{\mathrm{pad}}=c_{\mathrm{pad}}(\tau)$` is selected from the particle-quadrature error. The implementation uses a fixed `5.0` independent of `τ`. Add one sentence after line 857 recording that:

```latex
In the present implementation $c_{\mathrm{pad}}=5$ is used for all tolerances;
Figure~\ref{fig:tkm}(a) shows that this choice keeps
$\mathcal{E}_{\mathrm{far}}$ below the plateau set by $\tau$ over the full range
tested.
```

Only include this sentence if Task 8's `err_far_n` data actually supports it — check that the far-target error sits at or below the overall plateau for each `τ`. If it does not, report the discrepancy instead of writing the sentence.

- [ ] **Step 8: Build the paper**

```bash
cd ~/Articles/four_indices_bie
latexmk -pdf main.tex 2>&1 | tail -25
```

Expected: no new warnings about undefined references, and `eq:padding_ratio` still resolves. Confirm the Fig. 5 caption and §3.4 read correctly in `main.pdf`.

- [ ] **Step 9: Commit**

```bash
cd ~/Articles/four_indices_bie
git add main.tex figs/fig5_tkm_validation.pdf
git commit -m "$(cat <<'EOF'
sec3.4: validate the hybrid at c_pad = 5 with near and far targets

The experiment set h_n = 0 and used only targets inside B, so it did not
exercise the near/far switch, and the text said so. It now runs the complete
hybrid at c_pad = 5 over 400 targets straddling dB_pad, so those disclaimers
are gone.

Panel (b)'s eta now carries h_n, and the text no longer claims the transition
is exactly at eta = 1: Eq. (3.16) is sufficient, not sharp, and the test
Gaussian's numerical support is smaller than B, so the crossing sits slightly
below 1.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Merge the worktree

**Files:** none

- [ ] **Step 1: Re-run the full BoundaryIntegral suite one last time**

```bash
cd ~/codes/bi-cpad-align
julia --project=. test/runtests.jl 2>&1 | tail -30
```

Expected: green.

- [ ] **Step 2: Merge and clean up**

```bash
cd ~/codes/BoundaryIntegral.jl
git merge --no-ff cpad-paper-alignment
julia --project=. test/runtests.jl 2>&1 | tail -20
```

Only after the suite passes in the live checkout:

```bash
cd ~/codes/BoundaryIntegral.jl
git worktree remove ~/codes/bi-cpad-align
```

- [ ] **Step 3: Point `fig_gen` back at the live checkout**

```bash
cd ~/work/four_index_integral_solver/codes/fig_gen
julia --project=. -e 'using Pkg; Pkg.develop(path = "/mnt/home/xgao1/codes/BoundaryIntegral.jl"); Pkg.resolve()' 2>&1 | tail -10
```

- [ ] **Step 4: Commit the Manifest change if `fig_gen/Manifest.toml` is tracked**

```bash
cd ~/work/four_index_integral_solver/codes
git status --short fig_gen/
```

Commit only if the Manifest is tracked and changed.

---

## Self-Review

**Spec coverage:**

| spec item | task |
|---|---|
| Mismatch 1: `h` = `‖A_ρ‖₂` | 1 (basis field), 2 (`lattice_spacing`) |
| Mismatch 2: two `c_pad` values | 3, 4, 5 |
| Mismatch 3: ball → box | 4 |
| Mismatch 4: `L` and period over-estimated | 2 (`near_field_geometry`), 3 |
| Component 1: `VolumeSource` basis + propagation audit | 1 (steps 4-5), 5 (step 3) |
| Component 2: shared geometry helper | 2 |
| Component 3: three call sites + campaign plumbing | 3, 4, 5 |
| Component 4: Fig. 5 data script | 7 |
| Component 5: paper edits | 9 |
| Verification 1: full suite | 0 (baseline), 3, 4, 5, 6, 10 |
| Verification 2: cross-path equivalence | 6 |
| Verification 3: cubic-lattice unchanged | 6 |
| Verification 4: Fig. 5 smoke run | 7 step 4 |
| Verification 5: full sweep + `η` check | 8 |
| Bounding the un-rerun §6.4 change | 5 (step 1 test) |
| `B` = cell box | 2 (`source_box`), verified in 7 step 2 |
| No Stage C | plan has no rerun task |

**Type consistency:** `NearFieldGeometry` field names (`lo, hi, center, l, h, hn, L, dks`) are used identically in Tasks 2, 3, 5, 7. `c_pad` is a keyword everywhere (`PrecomputedVolumeField`, `near_field_geometry`, `_classify_near_far_targets`, `evaluate_batch_potential`, `rhs_dielectric_box3d_hybrid`). `with_density` and `_primitive_basis` are defined in Task 1 and consumed in Tasks 1 and 5. `in_near_region(g, targets, i)` has one signature throughout.

**Known open item, deliberately left to the implementer:** Task 4 Step 3 depends on the current body of `_classify_near_far_panels`, which builds per-panel representative points in a way this plan does not reproduce verbatim. The step says to read lines 85-103 and preserve the existing representative-point choice, extracting it into `_panel_representative_points` if inlined. This is a read-then-preserve instruction rather than a placeholder, but it is the one place the implementer must consult the source to write the code.

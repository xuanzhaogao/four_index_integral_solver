# PrecomputedVolumeField Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a target-independent precomputed spectral representation of a volume source's free-space potential/gradient to BoundaryIntegral.jl, so the RHS-adaptive mesh builder, RHS assembly, and volume-potential evaluation stop redoing a ~23 s type-1 NUFFT + FFT-replanning per call (measured 16x speedup of the adaptive build on the production monolayer problem; max rel deviation 1.4e-5; 0 refinement-decision flips — see `numerical_results/exp65_orbital_bench/REPORT.md` addendum and `scripts/bench_meshgen.jl` for the validated prototype this code is lifted from).

**Architecture:** New struct `PrecomputedVolumeField` (built once per source density): fixed evaluation box B = source bbox + `margin_h * h`; on B's target-independent Fourier box, type-1 NUFFT of all source charges + truncated-kernel diagonal scaling (+ spectral gradient coefficients) are computed once and stored. Every later evaluation routes targets inside B through a type-2 NUFFT on the stored coefficients (fixed FFT dims => FFTW wisdom reuse) and targets outside B through an exact direct threaded sum (measured 7.3e10 pairs/s — no FMM needed at these batch sizes). Consumers: a field overload of `single_dielectric_box3d_rhs_adaptive`, a field RHS assembly `rhs_dielectric_box3d_field`, and exported `volume_field_potential` / `volume_field_gradient`.

**Tech Stack:** Julia, TKM3D.jl internals (`combined_box_geometry_3xn`, `centered_mode_axis`, `truncated_laplace3d_hat`, `_spectral_gradient_coeffs_3d`, `_finufft_make_type2_plan_3d`), FINUFFT, Threads.

**CRITICAL — worktree rule:** All changes happen in a dedicated git worktree at
`$WT = /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf` (branch `feature/precomputed-volume-field`).
NEVER modify `/mnt/home/xgao1/codes/BoundaryIntegral.jl` (the user's live checkout) directly.
Julia binary: `/mnt/home/xgao1/.juliaup/bin/julia` (never `module load julia`).

**Normalization convention:** everything in this feature uses the free-space kernel `1/(4π r)`
(TKM3D convention). The FMM3D `1/r` convention appears nowhere in the new code.

---

## File structure

- Create: `$WT/src/shape/volume_field.jl` — struct, constructor, eval functions, direct-sum kernels, `_rhs_volume_targets_field`. (Lives in `shape/` because it depends on `_volume_source_fmm_sources` / `_estimate_source_spacing` / `_estimate_tkm3dc_kmax` from `shape/box3d_fmm_helpers.jl`.)
- Modify: `$WT/src/BoundaryIntegral.jl` — `include("shape/volume_field.jl")` after line 107 (`include("shape/box3d_fmm_helpers.jl")`), plus exports.
- Modify: `$WT/src/shape/box3d_rhs_adaptive.jl` — `_rhs_panel3d_resolved_volume_field` + field overload of `single_dielectric_box3d_rhs_adaptive` (append at end of file).
- Modify: `$WT/src/solver/dielectric_box3d.jl` — `rhs_dielectric_box3d_field` (append after `rhs_dielectric_box3d_hybrid`, i.e. after line ~216).
- Create: `$WT/test/shape/volume_field.jl` — all tests for this feature (one file, testsets appended task by task).
- Modify: `$WT/test/runtests.jl` — register the test file after `include("core/sources.jl")`.

---

### Task 1: Worktree setup and baseline

**Files:** none modified (setup only)

- [ ] **Step 1: Create the worktree and branch**

```bash
git -C /mnt/home/xgao1/codes/BoundaryIntegral.jl worktree add \
    /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf -b feature/precomputed-volume-field
```

Expected: `Preparing worktree (new branch 'feature/precomputed-volume-field')`.
If the worktree already exists from a previous attempt, reuse it (`git -C $WT status` must be clean).

- [ ] **Step 2: Provide a Manifest (gitignored in BI, so the worktree has none)**

```bash
cp /mnt/home/xgao1/codes/BoundaryIntegral.jl/Manifest.toml /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf/
```

The copied Manifest dev-paths TKM3D to `/mnt/home/xgao1/codes/TKM3D.jl/` (absolute), which is correct for the worktree too.

- [ ] **Step 3: Verify the package loads from the worktree**

```bash
cd /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf && \
/mnt/home/xgao1/.juliaup/bin/julia --project=. -e 'using Pkg; Pkg.instantiate(); using BoundaryIntegral; println("baseline OK")'
```

Expected: `baseline OK` (precompilation may take a few minutes the first time).

---

### Task 2: `PrecomputedVolumeField` struct, constructor, potential/gradient evaluation

**Files:**
- Create: `$WT/src/shape/volume_field.jl`
- Modify: `$WT/src/BoundaryIntegral.jl` (include at line 107, exports near line 24)
- Create: `$WT/test/shape/volume_field.jl`
- Modify: `$WT/test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Create `$WT/test/shape/volume_field.jl`:

```julia
# Reference O(N^2) sums in the free-space 1/(4π r) convention.
function _ref_potential(sources, charges, targets)
    n = size(targets, 2)
    out = zeros(n)
    for i in 1:n, j in eachindex(charges)
        dx = targets[1, i] - sources[1, j]
        dy = targets[2, i] - sources[2, j]
        dz = targets[3, i] - sources[3, j]
        out[i] += charges[j] / sqrt(dx^2 + dy^2 + dz^2)
    end
    return out ./ (4π)
end

function _ref_gradient(sources, charges, targets)
    n = size(targets, 2)
    out = zeros(3, n)
    for i in 1:n, j in eachindex(charges)
        dx = targets[1, i] - sources[1, j]
        dy = targets[2, i] - sources[2, j]
        dz = targets[3, i] - sources[3, j]
        r2 = dx^2 + dy^2 + dz^2
        s = charges[j] / (r2 * sqrt(r2))
        out[1, i] -= s * dx; out[2, i] -= s * dy; out[3, i] -= s * dz
    end
    return out ./ (4π)
end

@testset "PrecomputedVolumeField potential/gradient" begin
    gsrc = BoundaryIntegral.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 12, 1e-6)
    field = PrecomputedVolumeField(gsrc; tol = 1e-6)
    src, q = BoundaryIntegral._volume_source_fmm_sources(gsrc)

    # mixed batch: 4 targets inside the field box, 4 well outside it
    targets = [0.05  0.0   0.2  -0.15  3.0  0.0  -3.0  4.0;
               0.0   0.1  -0.1   0.05  0.0  3.5   2.0  4.0;
               0.0   0.02  0.1  -0.1   0.0  1.0  -2.0  4.0]
    @test count(i -> BoundaryIntegral.in_field_box(field, targets, i), 1:8) == 4

    pot = volume_field_potential(field, targets)
    grad = volume_field_gradient(field, targets)
    pref = _ref_potential(src, q, targets)
    gref = _ref_gradient(src, q, targets)

    # in-box: spectral path (kmax = grid Nyquist, tol 1e-6)
    for i in 1:4
        @test isapprox(pot[i], pref[i]; rtol = 1e-4)
        @test isapprox(grad[:, i], gref[:, i]; rtol = 1e-3)
    end
    # out-of-box: direct sum, identical arithmetic to the reference
    for i in 5:8
        @test isapprox(pot[i], pref[i]; rtol = 1e-12)
        @test isapprox(grad[:, i], gref[:, i]; rtol = 1e-12)
    end

    # gradient-only / potential-only construction guards
    fg = PrecomputedVolumeField(gsrc; tol = 1e-6, compute_pot = false)
    @test_throws ArgumentError volume_field_potential(fg, targets)
    fp = PrecomputedVolumeField(gsrc; tol = 1e-6, compute_grad = false)
    @test_throws ArgumentError volume_field_gradient(fp, targets)
end
```

Register in `$WT/test/runtests.jl` — directly after the line `include("core/sources.jl")`, add:

```julia
    include("shape/volume_field.jl")
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf && \
/mnt/home/xgao1/.juliaup/bin/julia --project=. -e \
  'using BoundaryIntegral, LinearAlgebra, Test; include("test/shape/volume_field.jl")'
```

Expected: FAIL with `UndefVarError: PrecomputedVolumeField not defined`.

- [ ] **Step 3: Write the implementation**

Create `$WT/src/shape/volume_field.jl`:

```julia
"""
    PrecomputedVolumeField(vs; tol, kmax = nothing, margin_h = 5.0,
                           compute_pot = true, compute_grad = true)

Target-independent precomputed spectral representation of the free-space
Laplace potential (`1/(4π r)` convention) of a `VolumeSource`, for evaluation
at many target batches without redoing per-call setup.

Construction (once): fix the evaluation box B = source bounding box extended
by `margin_h * h` (`h` = mean source spacing); on the Fourier box determined
by B alone, run the type-1 NUFFT of all source charges, apply the truncated
Laplace kernel `TKM3D.truncated_laplace3d_hat`, and (optionally) materialize
the spectral gradient coefficients.

Evaluation (per batch): targets inside B use a type-2 NUFFT on the stored
coefficients — the FFT dimensions never change, so FFTW planning is reused;
targets outside B use an exact direct threaded sum (no FMM setup).

This removes the dominant cost of the RHS-adaptive mesh build, where the
target-dependent Fourier box previously forced a fresh type-1 NUFFT + FFT
plan at every refinement depth (~23 s per call on the production monolayer
density; measured 16x build speedup with max rel deviation 1.4e-5 and no
refinement-decision changes).
"""
struct PrecomputedVolumeField{T <: AbstractFloat}
    sources::Matrix{T}
    charges::Vector{T}
    lo::NTuple{3, T}
    hi::NTuple{3, T}
    center::NTuple{3, T}
    dks::NTuple{3, T}
    nmodes::NTuple{3, Int}
    kmax::T
    tol::T
    prefactor::T
    coeff::Union{Nothing, Array{Complex{T}, 3}}
    grad_coeff::Union{Nothing, Array{Complex{T}, 4}}
end

function PrecomputedVolumeField(
    vs::VolumeSource{T, 3};
    tol::Real,
    kmax::Union{Nothing, Real} = nothing,
    margin_h::Real = 5.0,
    compute_pot::Bool = true,
    compute_grad::Bool = true,
) where {T <: AbstractFloat}
    (compute_pot || compute_grad) || throw(ArgumentError("nothing to precompute"))
    tolT = T(tol)
    tolT > zero(T) || throw(ArgumentError("tol must be positive"))
    sources, charges = _volume_source_fmm_sources(vs)
    h = _estimate_source_spacing(vs)
    km = isnothing(kmax) ? T(_estimate_tkm3dc_kmax(h)) : T(kmax)
    km > zero(T) || throw(ArgumentError("kmax must be positive"))

    m = T(margin_h) * h
    lo = ntuple(d -> minimum(view(sources, d, :)) - m, 3)
    hi = ntuple(d -> maximum(view(sources, d, :)) + m, 3)
    corners = T[lo[1] hi[1]; lo[2] hi[2]; lo[3] hi[3]]
    lengths, center = TKM3D.combined_box_geometry_3xn(sources, corners)
    Lbig = sqrt(sum(abs2, lengths))
    dks = ntuple(d -> prevfloat(T(2π) / (lengths[d] + Lbig)), 3)
    kx = TKM3D.centered_mode_axis(dks[1], km)
    ky = TKM3D.centered_mode_axis(dks[2], km)
    kz = TKM3D.centered_mode_axis(dks[3], km)
    nmodes = (length(kx), length(ky), length(kz))

    srcx = dks[1] .* (vec(view(sources, 1, :)) .- center[1])
    srcy = dks[2] .* (vec(view(sources, 2, :)) .- center[2])
    srcz = dks[3] .* (vec(view(sources, 3, :)) .- center[3])
    coeff0 = TKM3D.FINUFFT.nufft3d1(srcx, srcy, srcz, complex.(charges), -1, tolT, nmodes...)
    coeff = ndims(coeff0) == 4 ? dropdims(coeff0; dims = 4) : coeff0
    @inbounds for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2)
        coeff[ix, iy, iz] = k <= km ?
            coeff[ix, iy, iz] * TKM3D.truncated_laplace3d_hat(k, Lbig) :
            zero(eltype(coeff))
    end
    grad_coeff = compute_grad ? TKM3D._spectral_gradient_coeffs_3d(coeff, kx, ky, kz) : nothing
    prefactor = dks[1] * dks[2] * dks[3] / T(2π)^3
    return PrecomputedVolumeField{T}(
        sources, charges, lo, hi, (center[1], center[2], center[3]),
        dks, nmodes, km, tolT, prefactor,
        compute_pot ? coeff : nothing, grad_coeff)
end

@inline in_field_box(f::PrecomputedVolumeField, targets::AbstractMatrix, i::Integer) =
    (f.lo[1] <= targets[1, i] <= f.hi[1]) &&
    (f.lo[2] <= targets[2, i] <= f.hi[2]) &&
    (f.lo[3] <= targets[3, i] <= f.hi[3])

function _field_scaled_targets(f::PrecomputedVolumeField{T}, targets, idxs) where {T}
    txn = T[f.dks[1] * (targets[1, i] - f.center[1]) for i in idxs]
    tyn = T[f.dks[2] * (targets[2, i] - f.center[2]) for i in idxs]
    tzn = T[f.dks[3] * (targets[3, i] - f.center[3]) for i in idxs]
    return txn, tyn, tzn
end

function _direct_potential!(out::AbstractVector{T}, sources, charges, targets, idxs) where {T}
    isempty(idxs) && return out
    sx = sources[1, :]; sy = sources[2, :]; sz = sources[3, :]
    Threads.@threads for ii in eachindex(idxs)
        i = idxs[ii]
        x, y, z = targets[1, i], targets[2, i], targets[3, i]
        acc = zero(T)
        @inbounds @simd for j in eachindex(charges)
            dx = x - sx[j]; dy = y - sy[j]; dz = z - sz[j]
            acc += charges[j] / sqrt(dx * dx + dy * dy + dz * dz)
        end
        out[i] = acc / (4 * T(π))
    end
    return out
end

function _direct_gradient!(out::AbstractMatrix{T}, sources, charges, targets, idxs) where {T}
    isempty(idxs) && return out
    sx = sources[1, :]; sy = sources[2, :]; sz = sources[3, :]
    Threads.@threads for ii in eachindex(idxs)
        i = idxs[ii]
        x, y, z = targets[1, i], targets[2, i], targets[3, i]
        gx = zero(T); gy = zero(T); gz = zero(T)
        @inbounds @simd for j in eachindex(charges)
            dx = x - sx[j]; dy = y - sy[j]; dz = z - sz[j]
            r2 = dx * dx + dy * dy + dz * dz
            s = charges[j] / (r2 * sqrt(r2))
            gx -= s * dx; gy -= s * dy; gz -= s * dz
        end
        out[1, i] = gx / (4 * T(π))
        out[2, i] = gy / (4 * T(π))
        out[3, i] = gz / (4 * T(π))
    end
    return out
end

"""
    volume_field_potential(field, targets) -> Vector

Potential of the precomputed field at `targets` (3 x n), free-space
`1/(4π r)` normalization. In-box targets via type-2 NUFFT, out-of-box via
direct threaded summation.
"""
function volume_field_potential(f::PrecomputedVolumeField{T}, targets::AbstractMatrix{<:Real}) where {T}
    f.coeff === nothing && throw(ArgumentError("field was built with compute_pot = false"))
    trg = Matrix{T}(targets)
    n = size(trg, 2)
    out = Vector{T}(undef, n)
    inb = [in_field_box(f, trg, i) for i in 1:n]
    bidx = findall(inb); oidx = findall(!, inb)
    if !isempty(bidx)
        txn, tyn, tzn = _field_scaled_targets(f, trg, bidx)
        plan = TKM3D._finufft_make_type2_plan_3d(txn, tyn, tzn, 1, f.tol, f.nmodes, 1, T)
        vals = try
            TKM3D._finufft_type2_exec_3d(plan, f.coeff)
        finally
            TKM3D.FINUFFT.finufft_destroy!(plan)
        end
        for (k, i) in enumerate(bidx)
            out[i] = f.prefactor * real(vals[k])
        end
    end
    _direct_potential!(out, f.sources, f.charges, trg, oidx)
    return out
end

"""
    volume_field_gradient(field, targets) -> 3 x n Matrix

Gradient of the precomputed field at `targets` (3 x n), free-space
`1/(4π r)` normalization.
"""
function volume_field_gradient(f::PrecomputedVolumeField{T}, targets::AbstractMatrix{<:Real}) where {T}
    f.grad_coeff === nothing && throw(ArgumentError("field was built with compute_grad = false"))
    trg = Matrix{T}(targets)
    n = size(trg, 2)
    out = Matrix{T}(undef, 3, n)
    inb = [in_field_box(f, trg, i) for i in 1:n]
    bidx = findall(inb); oidx = findall(!, inb)
    if !isempty(bidx)
        txn, tyn, tzn = _field_scaled_targets(f, trg, bidx)
        plan = TKM3D._finufft_make_type2_plan_3d(txn, tyn, tzn, 1, f.tol, f.nmodes, 3, T)
        vals = try
            TKM3D._finufft_type2_exec_3d(plan, f.grad_coeff)   # n_in x 3
        finally
            TKM3D.FINUFFT.finufft_destroy!(plan)
        end
        for (k, i) in enumerate(bidx)
            out[1, i] = f.prefactor * real(vals[k, 1])
            out[2, i] = f.prefactor * real(vals[k, 2])
            out[3, i] = f.prefactor * real(vals[k, 3])
        end
    end
    _direct_gradient!(out, f.sources, f.charges, trg, oidx)
    return out
end
```

Wire into `$WT/src/BoundaryIntegral.jl`:
1. After the line `include("shape/box3d_fmm_helpers.jl") # FMM/TKM helpers (needed by rhs_adaptive)` (line 107), add:

```julia
include("shape/volume_field.jl") # precomputed spectral volume field
```

2. Next to the existing `export PointSource, VolumeSource` block (line 22), add:

```julia
export PrecomputedVolumeField, volume_field_potential, volume_field_gradient
```

- [ ] **Step 4: Run test to verify it passes**

Same command as Step 2. Expected: all `PrecomputedVolumeField potential/gradient` tests PASS.

- [ ] **Step 5: Commit**

```bash
cd /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf && \
git add src/shape/volume_field.jl src/BoundaryIntegral.jl test/shape/volume_field.jl test/runtests.jl && \
git commit -m "feat: PrecomputedVolumeField — precompute-once spectral volume field

Type-1 NUFFT + truncated-kernel scaling (+ spectral gradient coeffs) on a
fixed density-defined box, evaluated per batch via type-2 NUFFT (in box) or
direct threaded sum (out of box). Validated against the exp65 meshgen
prototype (16x adaptive-build speedup, rel dev 1.4e-5)."
```

---

### Task 3: `_rhs_volume_targets_field` (RHS values from a field)

**Files:**
- Modify: `$WT/src/shape/volume_field.jl` (append)
- Modify: `$WT/test/shape/volume_field.jl` (append)

- [ ] **Step 1: Write the failing test** — append to `$WT/test/shape/volume_field.jl`:

```julia
@testset "_rhs_volume_targets_field vs hybrid" begin
    gsrc = BoundaryIntegral.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 12, 1e-6)
    field = PrecomputedVolumeField(gsrc; tol = 1e-6, compute_pot = false)
    src, q = BoundaryIntegral._volume_source_fmm_sources(gsrc)
    h = BoundaryIntegral._estimate_source_spacing(gsrc)
    kmax = BoundaryIntegral._estimate_tkm3dc_kmax(h)

    targets = [0.05  0.2  -0.15  3.0  -3.0  4.0;
               0.0  -0.1   0.05  0.0   2.0  4.0;
               0.0   0.1  -0.1   0.0  -2.0  4.0]
    normals = [1.0  0.0  0.0  0.0  1.0  0.0;
               0.0  1.0  0.0  0.0  0.0  1.0;
               0.0  0.0  1.0  1.0  0.0  0.0]

    rhs_field = BoundaryIntegral._rhs_volume_targets_field(field, targets, normals, 1.0)

    is_near = BoundaryIntegral._classify_near_far_targets(targets, gsrc, h)
    rhs_hyb, _, _ = BoundaryIntegral._rhs_volume_targets_hybrid(
        src, q, targets, normals, 1.0, 1e-6, kmax, is_near)

    scale = maximum(abs, rhs_hyb)
    for i in 1:6
        @test abs(rhs_field[i] - rhs_hyb[i]) <= 1e-4 * scale
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Same command as Task 2 Step 2. Expected: FAIL with `UndefVarError: _rhs_volume_targets_field not defined`.

- [ ] **Step 3: Implement** — append to `$WT/src/shape/volume_field.jl`:

```julia
"""
    _rhs_volume_targets_field(field, targets, normals, eps_src) -> Vector

Dielectric-interface RHS values `-(n · ∇φ) / eps_src` at `targets` from a
precomputed field (free-space normalization, matching the TKM branch of
`_rhs_volume_targets_hybrid`). Replaces the per-call KDTree classify +
lfmm3d/ltkm3dc pair with one field evaluation.
"""
function _rhs_volume_targets_field(
    field::PrecomputedVolumeField{T},
    targets::Matrix{T},
    normals::Matrix{T},
    eps_src::T,
) where {T}
    n = size(targets, 2)
    @assert size(normals, 2) == n
    grad = volume_field_gradient(field, targets)
    rhs_vals = Vector{T}(undef, n)
    @inbounds for i in 1:n
        rhs_vals[i] = -(normals[1, i] * grad[1, i] + normals[2, i] * grad[2, i] +
                        normals[3, i] * grad[3, i]) / eps_src
    end
    return rhs_vals
end
```

- [ ] **Step 4: Run test to verify it passes** — same command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf && \
git add src/shape/volume_field.jl test/shape/volume_field.jl && \
git commit -m "feat: _rhs_volume_targets_field — interface RHS from precomputed field"
```

---

### Task 4: Field overload of the RHS-adaptive builder

**Files:**
- Modify: `$WT/src/shape/box3d_rhs_adaptive.jl` (append both functions at end of file)
- Modify: `$WT/test/shape/volume_field.jl` (append)

- [ ] **Step 1: Write the failing test** — append to `$WT/test/shape/volume_field.jl`:

```julia
@testset "adaptive builder: field overload reproduces vs path" begin
    gsrc = BoundaryIntegral.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 12, 1e-6)
    field = PrecomputedVolumeField(gsrc; tol = 1e-4, compute_pot = false)

    iface_vs = BoundaryIntegral.single_dielectric_box3d_rhs_adaptive(
        4.0, 4.0, 1.0, 4, gsrc, 1.0, 0.26, 1e-3, 3.0, 1.0, Float64; max_depth = 8)
    iface_f = BoundaryIntegral.single_dielectric_box3d_rhs_adaptive(
        4.0, 4.0, 1.0, 4, field, 1.0, 0.26, 1e-3, 3.0, 1.0, Float64; max_depth = 8)

    @test length(iface_f.panels) == length(iface_vs.panels)
    @test BoundaryIntegral.num_points(iface_f) == BoundaryIntegral.num_points(iface_vs)
    # identical refinement path => identical panel geometry
    c_vs = sort([sum(p.corners[1]) + 3 * sum(p.corners[3]) for p in iface_vs.panels])
    c_f  = sort([sum(p.corners[1]) + 3 * sum(p.corners[3]) for p in iface_f.panels])
    @test isapprox(c_vs, c_f; rtol = 1e-12)
end
```

- [ ] **Step 2: Run test to verify it fails**

Same command. Expected: FAIL with `MethodError: no method matching single_dielectric_box3d_rhs_adaptive(::Float64, ::Float64, ::Float64, ::Int64, ::PrecomputedVolumeField{Float64}, ...)`.

- [ ] **Step 3: Implement** — append to `$WT/src/shape/box3d_rhs_adaptive.jl`:

```julia
# --- PrecomputedVolumeField path: same refinement logic, but the per-depth
# --- RHS evaluation uses the stored spectral field (no per-call type-1
# --- NUFFT / FFT replanning / KDTree classify / FMM setup).

function _rhs_panel3d_resolved_volume_field(
    panels::Vector{TempPanel3D{T}},
    field::PrecomputedVolumeField{T},
    eps_src::T,
    ns::Vector{T},
    ws::Vector{T},
    atol::T,
) where {T}
    n_panels = length(panels)
    if n_panels == 0
        return Bool[]
    end
    resolved = fill(false, n_panels)
    n_quad = length(ns)
    λ = gl_barycentric_weights(ns, ws)
    n_pts = 10
    xs = range(-one(T), one(T); length = n_pts)
    ys = range(-one(T), one(T); length = n_pts)
    targets, normals, n_per_panel = _rhs_panel3d_refinement_targets(panels, ns, ws; n_pts = n_pts)

    rhs_vals = _rhs_volume_targets_field(field, targets, normals, eps_src)
    @info "    rhs panel field evaluation, targets: $(size(targets, 2))"

    idx = 0
    for p in 1:n_panels
        quad_vals = Matrix{T}(undef, n_quad, n_quad)
        for i in 1:n_quad
            for j in 1:n_quad
                idx += 1
                quad_vals[i, j] = rhs_vals[idx]
            end
        end
        err = zero(T)
        for u in xs
            rx = T.(barycentric_row(ns, λ, u))
            for v in ys
                ry = T.(barycentric_row(ns, λ, v))
                approx = zero(T)
                for i in 1:n_quad
                    for j in 1:n_quad
                        approx += quad_vals[i, j] * rx[i] * ry[j]
                    end
                end
                idx += 1
                err = max(err, abs(rhs_vals[idx] - approx))
            end
        end
        resolved[p] = err <= atol
    end
    return resolved
end

function single_dielectric_box3d_rhs_adaptive(
    Lx::T,
    Ly::T,
    Lz::T,
    n_quad::Int,
    field::PrecomputedVolumeField{T},
    eps_src::T,
    l_ec::T,
    rhs_atol::T,
    eps_in::T,
    eps_out::T,
    ::Type{T} = Float64;
    max_depth::Int = 128,
    alpha::T = sqrt(T(2)),
) where {T}
    ns, ws = gausslegendre(n_quad)
    @info "box3d volume field rhs adaptive panel generation, source points: $(length(field.charges))"

    solved = TempPanel3D{T}[]
    unsolved = _box3d_rhs_adaptive_initial_panels(Lx, Ly, Lz, alpha)
    depth = 0
    while !isempty(unsolved) && depth < max_depth
        @info "  depth $depth, unsolved panels: $(length(unsolved))"
        resolved = _rhs_panel3d_resolved_volume_field(unsolved, field, eps_src, ns, ws, rhs_atol)
        next_unsolved = TempPanel3D{T}[]
        for i in eachindex(unsolved)
            tpl = unsolved[i]
            if resolved[i]
                push!(solved, tpl)
            else
                append!(next_unsolved, divide_temp_panel3d(tpl, 2, 2))
            end
        end
        unsolved = next_unsolved
        depth += 1
    end
    append!(solved, unsolved)

    rough_ec = copy(solved)
    refined = TempPanel3D{T}[]
    while !isempty(rough_ec)
        tpl = popfirst!(rough_ec)
        has_ec = tpl.is_a_corner || tpl.is_b_corner || tpl.is_c_corner || tpl.is_d_corner ||
            tpl.is_ab_edge || tpl.is_bc_edge || tpl.is_cd_edge || tpl.is_da_edge
        L_ab = norm(tpl.b .- tpl.a)
        L_da = norm(tpl.a .- tpl.d)
        if has_ec && max(L_ab, L_da) > l_ec
            append!(rough_ec, divide_temp_panel3d(tpl, 2, 2))
        else
            push!(refined, tpl)
        end
    end

    panels = Vector{FlatPanel{T, 3}}()
    for tpl in refined
        is_edge = tpl.is_ab_edge || tpl.is_bc_edge || tpl.is_cd_edge || tpl.is_da_edge ||
            tpl.is_a_corner || tpl.is_b_corner || tpl.is_c_corner || tpl.is_d_corner
        push!(panels, rect_panel3d_discretize(tpl.a, tpl.b, tpl.c, tpl.d, ns, ws, tpl.normal; is_edge = is_edge))
    end

    return DielectricInterface(panels, fill(eps_in, length(panels)), fill(eps_out, length(panels)))
end
```

- [ ] **Step 4: Run test to verify it passes** — same command. Expected: PASS.
(If the panel-count assertion fails by a small number, a refinement decision flipped on a borderline panel; report this rather than loosening the test — at production parameters the prototype measured 0 flips.)

- [ ] **Step 5: Commit**

```bash
cd /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf && \
git add src/shape/box3d_rhs_adaptive.jl test/shape/volume_field.jl && \
git commit -m "feat: single_dielectric_box3d_rhs_adaptive overload on PrecomputedVolumeField"
```

---

### Task 5: Field RHS assembly for the solver

**Files:**
- Modify: `$WT/src/solver/dielectric_box3d.jl` (append after `rhs_dielectric_box3d_hybrid`, ~line 216)
- Modify: `$WT/src/BoundaryIntegral.jl` (export)
- Modify: `$WT/test/shape/volume_field.jl` (append)

- [ ] **Step 1: Write the failing test** — append to `$WT/test/shape/volume_field.jl`:

```julia
@testset "rhs_dielectric_box3d_field vs hybrid" begin
    gsrc = BoundaryIntegral.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 12, 1e-6)
    field = PrecomputedVolumeField(gsrc; tol = 1e-6, compute_pot = false)
    iface = BoundaryIntegral.single_dielectric_box3d_rhs_adaptive(
        4.0, 4.0, 1.0, 4, gsrc, 1.0, 0.26, 1e-3, 3.0, 1.0, Float64; max_depth = 8)

    rhs_f = rhs_dielectric_box3d_field(iface, field, 1.0)
    rhs_h = BoundaryIntegral.rhs_dielectric_box3d_hybrid(iface, gsrc, 1.0, 1e-6)

    @test length(rhs_f) == BoundaryIntegral.num_points(iface)
    @test maximum(abs, rhs_f .- rhs_h) <= 1e-4 * maximum(abs, rhs_h)
end
```

- [ ] **Step 2: Run test to verify it fails**

Same command. Expected: FAIL with `UndefVarError: rhs_dielectric_box3d_field not defined`.

- [ ] **Step 3: Implement** — append to `$WT/src/solver/dielectric_box3d.jl` (after the `Rhs_dielectric_box3d_hybrid` backward-compat alias, line ~216):

```julia
"""
    rhs_dielectric_box3d_field(interface, field, eps_src) -> Vector

Dielectric-interface RHS assembled from a `PrecomputedVolumeField`
(in-box targets: type-2 NUFFT on stored coefficients; out-of-box: direct
threaded summation). Drop-in alternative to `rhs_dielectric_box3d_hybrid`
when the same source density is used for many assemblies/evaluations.
"""
function rhs_dielectric_box3d_field(
    interface::DielectricInterface{P, Float64},
    field::PrecomputedVolumeField{Float64},
    eps_src::Float64,
) where {P <: AbstractPanel}
    n_points = num_points(interface)
    targets = Matrix{Float64}(undef, 3, n_points)
    normals = Matrix{Float64}(undef, 3, n_points)
    for (i, point) in enumerate(eachpoint(interface))
        targets[1, i] = point.panel_point.point[1]
        targets[2, i] = point.panel_point.point[2]
        targets[3, i] = point.panel_point.point[3]
        normals[1, i] = point.panel_point.normal[1]
        normals[2, i] = point.panel_point.normal[2]
        normals[3, i] = point.panel_point.normal[3]
    end
    return _rhs_volume_targets_field(field, targets, normals, eps_src)
end
```

In `$WT/src/BoundaryIntegral.jl`, extend the export added in Task 2 to:

```julia
export PrecomputedVolumeField, volume_field_potential, volume_field_gradient,
       rhs_dielectric_box3d_field
```

- [ ] **Step 4: Run test to verify it passes** — same command. Expected: PASS (all four testsets).

- [ ] **Step 5: Commit**

```bash
cd /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf && \
git add src/solver/dielectric_box3d.jl src/BoundaryIntegral.jl test/shape/volume_field.jl && \
git commit -m "feat: rhs_dielectric_box3d_field — RHS assembly from precomputed field"
```

---

### Task 6: Regression check of touched subsystems

**Files:** none modified

- [ ] **Step 1: Run the existing test files for every modified subsystem plus the new file**

```bash
cd /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf && \
/mnt/home/xgao1/.juliaup/bin/julia --project=. -e '
using BoundaryIntegral, LinearAlgebra, Test
@testset "regression" begin
    include("test/core/sources.jl")
    include("test/shape/volume_field.jl")
end'
```

Expected: all PASS. (The full `test/runtests.jl` suite is long; the user runs it before merging. If `test/shape/` contains files exercising `box3d_rhs_adaptive` or `test/solver/` files exercising `dielectric_box3d`, include those too — check with `ls test/shape test/solver` and add the relevant `include`s to the command.)

- [ ] **Step 2: Verify the live checkout is untouched and summarize**

```bash
git -C /mnt/home/xgao1/codes/BoundaryIntegral.jl status --short    # must be empty
git -C /mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf log --oneline main..HEAD
```

Expected: empty status; 4 feature commits listed. Report the commit list — merging is the user's call.

---

## After merge (not part of this plan; orchestrator follow-up)

- Production-scale validation: rerun `numerical_results/exp65_orbital_bench/scripts/run_single_rhs.jl` with a field-based variant and confirm the projected ~70 s end-to-end (vs 296 s) and unchanged u_onsite to ~1e-5.
- Optional: multi-box (`multi_dielectric_box3d_rhs_adaptive`) field overload — same pattern (YAGNI for the monolayer production path).
- Optional: profile/remove the residual per-call cost in `TKM3D.ltkm3dc` (unused ntrans=1 plan creation on the gradient-only path).

# Screened Monolayer Graphene at k_323201 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce CoQui's screened `U_ijkl` (4 distinct values: onsite, nn, hund_sf, hund_ph) for the monolayer graphene k_323201 Wannier dataset by running the bilayer slab GMRES dielectric solver as a black box, with monolayer geometry (`LZ = 3.35 Å`, `eps_in = 2.4`, `eps_out = 1.0`, orbital at slab midplane).

**Architecture:** Two GMRES dielectric solves total — one density-source (covers onsite + nn over a Sharp + 4-bandwidth sweep) and one Hund's-source (signed `phi1·phi2`, Sharp mode only). The new `MonolayerScreenedSolve` wrapper composes `bilayer_slab/src/ScreenedOrbitalSolve.jl`'s `solve_screened_mode` with monolayer-specific source/target construction; no changes to the bilayer solver.

**Tech Stack:** Julia 1.11 via juliaup, shared `Project.toml` (`BoundaryIntegral.jl`, `Krylov.jl`, `CSV.jl`, `DataFrames.jl`, `Test.jl`). XSFs and CoQui reference: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/`.

---

## File Structure

Paths relative to `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/`.

**New files:**

- `monolayer/src/MonolayerScreenedSolve.jl` — wrapper module providing:
  - `parse_coqui_loc(path)` — CoQui `_coqui_crpa_loc.out` parser; returns `Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}` for the 4 channels.
  - `monolayer_screened_sources(; orbital_1, orbital_2, source_tol, z_center)` — wraps the existing `MonolayerOrbitalLoader.centered_monolayer_sources_padded` and applies a z-only shift so the orbital sits at slab midplane; also builds the `phi1·phi2` product VolumeSource.
  - `density_target_specs(vs1, vs2)` — produces a 2-entry `target_specs` list (`:onsite` → target=vs1, `:nn` → target=vs2) with `Na = Nb = Nphi1` for onsite and `Na = Nphi1, Nb = Nphi2` for nn (matching the bare-channel Hund's-safe normalization).
  - `hund_target_specs(vs_product, vs1, vs2)` — produces a 1-entry `target_specs` list (`:hund` → target=vs_product) with `Na = Nphi1, Nb = Nphi2`.
  - `screened_run_record(...)` — packs one `pair_result` from `solve_screened_mode` + slab/mode metadata into a flat NamedTuple ready for CSV.
- `monolayer/scripts/screened_monolayer_k323201.jl` — driver: load sources, run density sweep + Hund's Sharp, parse CoQui, write CSV, print summary.
- `monolayer/data/screened_monolayer_k323201.csv` — output (11 rows: 2 density channels × 5 modes + 1 Hund's channel × 1 mode).
- `monolayer/results/2026-05-13-screened-monolayer-k323201-report.md` — dated report.
- `monolayer/test/test_monolayer_screened_solve.jl` — TDD tests for the parser + source builders. Cheap to run.

**Existing test wiring updated:**

- `monolayer/test/runtests.jl` — append `include("test_monolayer_screened_solve.jl")`.

**Reused unchanged (black box):**

- `bilayer_slab/src/ScreenedOrbitalSolve.jl` — `solve_screened_mode` plus `BI.SharpScreening` / `BI.SoftMixInversePermittivity` types pulled in via the include chain.
- `monolayer/src/MonolayerOrbitalLoader.jl` — `centered_monolayer_sources_padded`, `load_signed_xsf`, `load_squared_xsf`, `shift_datagrid`.
- `monolayer/src/MonolayerBareIntegrals.jl` — not consumed at runtime; the inline product-datagrid helper in `MonolayerScreenedSolve` duplicates its 4-line `_product_datagrid` body to avoid a cross-module private-symbol dependency.

---

## CoQui reference table (verified during planning)

From `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out`:

```
bare (v) and cRPA screened [U(iw=0)] interactions in eV
  i  j  k  l    v_ijkl    U_ijkl
--------------------------------
  0  0  0  0   17.4029    9.7824
  0  0  0  1    0.1277    0.0793
  0  0  1  0    0.1277    0.0793
  0  0  1  1    8.8319    5.1678
  0  1  0  1    0.1316    0.0937
  0  1  1  0    0.1316    0.0937
  1  1  1  1   17.4043    9.7829
  ...
```

The 4 distinct (i,j,k,l) tuples we parse:

- `(0,0,0,0)` → `:onsite`
- `(0,0,1,1)` → `:nn`
- `(0,1,0,1)` → `:hund_ph`
- `(0,1,1,0)` → `:hund_sf`

---

## Bilayer solver interface (consumed black-box)

Verified at planning time from `bilayer_slab/src/ScreenedOrbitalSolve.jl`:

```julia
function solve_screened_mode(
    source_vs::BI.VolumeSource{Float64, 3},
    target_specs,                          # vector of NamedTuples with fields :pair, :target_vs, :Na, :Nb
    Lx::Real, Ly::Real, Lz::Real,
    eps_in::Real, eps_out::Real,
    mode;                                  # BI.SharpScreening() or BI.SoftMixInversePermittivity(bw)
    n_quad::Integer = 6, edge_refine_level::Integer = 4,
    rhs_tol::Real = 1e-4, lhs_tol::Real = 1e-6,
    gmres_atol::Real = 1e-6, gmres_rtol::Real = 1e-6,
    max_order::Integer = 128, max_depth::Integer = 100,
    volume_tol::Real = 1e-6, scatter_range_factor::Real = 5.0, tkm_kmax = nothing,
)
```

Returns a NamedTuple with fields `interface, screened_vs, sigma, gmres_status, residual, tkm_kmax, n_interface_points, pair_results`. Each `pair_result` NamedTuple carries `pair (String), target_vs, Na, Nb, u_int_raw, u_scatter_raw, u_total_raw, u_int_ev, u_scatter_ev, u_total_ev`.

---

## Task 1: Inspect reference files and confirm format

**Files:** none (read-only)

- [ ] **Step 1: Confirm input files are readable**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
DATASET_DIR=/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15
for f in graphene_00001.xsf graphene_00002.xsf _coqui_crpa_loc.out; do
  p=$DATASET_DIR/$f
  [ -r "$p" ] && echo "OK   $f" || echo "MISS $f"
done
```

Expected: all three lines `OK`. If any file is unreadable, stop and report BLOCKED.

- [ ] **Step 2: Confirm the 4 reference U values**

```bash
grep -E "^\s*0\s+(0\s+0\s+0|0\s+1\s+1|1\s+0\s+1|1\s+1\s+0)" /mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out
```

Expected output (verified at planning time):

```
  0  0  0  0   17.4029    9.7824
  0  0  1  1    8.8319    5.1678
  0  1  0  1    0.1316    0.0937
  0  1  1  0    0.1316    0.0937
```

If any value differs, stop and report.

---

## Task 2: CoQui loc.out parser (TDD)

**Files:**
- Create: `monolayer/src/MonolayerScreenedSolve.jl` (parser only at this stage)
- Create: `monolayer/test/test_monolayer_screened_solve.jl`
- Modify: `monolayer/test/runtests.jl` (append one include)

- [ ] **Step 1: Write the failing test**

Create `monolayer/test/test_monolayer_screened_solve.jl`:

```julia
using Test

include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const COQUI_OUT = joinpath(REF_DIR, "_coqui_crpa_loc.out")
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")

@testset "parse_coqui_loc" begin
    @test isfile(COQUI_OUT)
    out = parse_coqui_loc(COQUI_OUT)
    @test sort(collect(keys(out))) == [:hund_ph, :hund_sf, :nn, :onsite]
    @test out[:onsite].v_ijkl ≈ 17.4029 atol = 1e-4
    @test out[:onsite].U_ijkl ≈ 9.7824  atol = 1e-4
    @test out[:nn].v_ijkl     ≈ 8.8319  atol = 1e-4
    @test out[:nn].U_ijkl     ≈ 5.1678  atol = 1e-4
    @test out[:hund_ph].v_ijkl ≈ 0.1316 atol = 1e-4
    @test out[:hund_ph].U_ijkl ≈ 0.0937 atol = 1e-4
    @test out[:hund_sf].v_ijkl ≈ 0.1316 atol = 1e-4
    @test out[:hund_sf].U_ijkl ≈ 0.0937 atol = 1e-4
end
```

Append to `monolayer/test/runtests.jl`:

```julia
include("test_monolayer_screened_solve.jl")
```

- [ ] **Step 2: Run and confirm FAIL**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: FAIL — `MonolayerScreenedSolve` module not found.

- [ ] **Step 3: Implement the parser**

Create `monolayer/src/MonolayerScreenedSolve.jl`:

```julia
module MonolayerScreenedSolve

export parse_coqui_loc

"""
    parse_coqui_loc(path)

Parse a CoQui `_coqui_crpa_loc.out` file's `i j k l v_ijkl U_ijkl` table and
return a Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}
for the four distinct channels:

- (0,0,0,0) → :onsite
- (0,0,1,1) → :nn
- (0,1,0,1) → :hund_ph
- (0,1,1,0) → :hund_sf
"""
function parse_coqui_loc(path::AbstractString)
    isfile(path) || throw(ArgumentError("CoQui loc file not found: $path"))
    targets = Dict{NTuple{4, Int}, Symbol}(
        (0, 0, 0, 0) => :onsite,
        (0, 0, 1, 1) => :nn,
        (0, 1, 0, 1) => :hund_ph,
        (0, 1, 1, 0) => :hund_sf,
    )
    result = Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}()
    pat = r"^\s*(\d)\s+(\d)\s+(\d)\s+(\d)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
    for line in eachline(path)
        m = match(pat, line)
        m === nothing && continue
        ijkl = (parse(Int, m.captures[1]), parse(Int, m.captures[2]),
                parse(Int, m.captures[3]), parse(Int, m.captures[4]))
        sym = get(targets, ijkl, nothing)
        sym === nothing && continue
        v = parse(Float64, m.captures[5])
        U = parse(Float64, m.captures[6])
        result[sym] = (v_ijkl = v, U_ijkl = U)
    end
    for k in (:onsite, :nn, :hund_ph, :hund_sf)
        haskey(result, k) || error("parse_coqui_loc: missing channel $k in $path")
    end
    return result
end

end # module
```

- [ ] **Step 4: Run and confirm PASS**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: all tests pass (including pre-existing tests).

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerScreenedSolve.jl monolayer/test/test_monolayer_screened_solve.jl monolayer/test/runtests.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): CoQui loc.out parser for U_ijkl channels

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Source construction (z-shifted vs1/vs2 + product VS)

**Files:**
- Modify: `monolayer/src/MonolayerScreenedSolve.jl`
- Modify: `monolayer/test/test_monolayer_screened_solve.jl`

- [ ] **Step 1: Write the failing test**

Append to `monolayer/test/test_monolayer_screened_solve.jl`:

```julia
@testset "monolayer_screened_sources — z-centering and product source" begin
    out = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = 1e-3, z_center = 1.675,
    )
    @test hasproperty(out, :vs1)
    @test hasproperty(out, :vs2)
    @test hasproperty(out, :vs_product)
    @test hasproperty(out, :Nphi1)
    @test hasproperty(out, :Nphi2)

    # vs1 density centroid should be at (0, 0, z_center) since the loader
    # zeros out the in-plane centroid and we add z_center on top.
    weights = out.vs1.weights .* out.vs1.density
    total = sum(weights)
    centroid = (
        sum(out.vs1.positions[1, :] .* weights) / total,
        sum(out.vs1.positions[2, :] .* weights) / total,
        sum(out.vs1.positions[3, :] .* weights) / total,
    )
    @test abs(centroid[1]) < 1e-3
    @test abs(centroid[2]) < 1e-3
    @test abs(centroid[3] - 1.675) < 1e-3

    # vs_product integrates to ~0 (orthogonal Wannier orbitals).
    n_prod = sum(out.vs_product.weights .* out.vs_product.density)
    @test abs(n_prod) < 0.1

    # Nphi1, Nphi2 are the squared-orbital norms (~1 each for a normalized Wannier function).
    @test 0.8 < out.Nphi1 < 1.2
    @test 0.8 < out.Nphi2 < 1.2
end
```

- [ ] **Step 2: Run and confirm FAIL**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: FAIL — `monolayer_screened_sources` not defined.

- [ ] **Step 3: Implement source construction**

Append to `monolayer/src/MonolayerScreenedSolve.jl` (before `end # module`):

```julia
using BoundaryIntegral
import BoundaryIntegral as BI

include(joinpath(@__DIR__, "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

export monolayer_screened_sources

function _product_datagrid(a, b)
    size(a.values) == size(b.values) || throw(ArgumentError("grid shapes differ"))
    a.origin == b.origin || throw(ArgumentError("grid origins differ"))
    values = a.values .* b.values
    return merge(a, (; values = values))
end

function _shift_z(datagrid, dz::Real)
    origin = (Float64(datagrid.origin[1]),
              Float64(datagrid.origin[2]),
              Float64(datagrid.origin[3]) + Float64(dz))
    return merge(datagrid, (; origin = origin))
end

"""
    monolayer_screened_sources(; orbital_1, orbital_2, source_tol, z_center)

Load and center the two monolayer Wannier orbitals, then apply an additional
z-only shift of `+z_center` so the squared-orbital-1 centroid lands at
`(0, 0, z_center)` (orbital sits at slab midplane). Returns the three
`VolumeSource` objects (`vs1`, `vs2`, `vs_product = phi1*phi2`) plus the
two squared-orbital norms (`Nphi1`, `Nphi2`) used for the eV conversion.
"""
function monolayer_screened_sources(;
    orbital_1::AbstractString,
    orbital_2::AbstractString,
    source_tol::Real = 1e-3,
    z_center::Real = 0.0,
)
    sq = MonolayerOrbitalLoader.centered_monolayer_sources_padded(;
        orbital_1 = orbital_1, orbital_2 = orbital_2,
        source_tol = source_tol, square = true, mirror_pad_level = 0,
    )
    sg = MonolayerOrbitalLoader.centered_monolayer_sources_padded(;
        orbital_1 = orbital_1, orbital_2 = orbital_2,
        source_tol = source_tol, square = false, mirror_pad_level = 0,
    )

    dg1_sq_z = _shift_z(sq.datagrid_1, z_center)
    dg2_sq_z = _shift_z(sq.datagrid_2, z_center)
    dg1_sg_z = _shift_z(sg.datagrid_1, z_center)
    dg2_sg_z = _shift_z(sg.datagrid_2, z_center)

    vs1 = BI.VolumeSource(dg1_sq_z, tol = source_tol)
    vs2 = BI.VolumeSource(dg2_sq_z, tol = source_tol)

    product_dg = _product_datagrid(dg1_sg_z, dg2_sg_z)
    vs_product = BI.VolumeSource(product_dg, tol = source_tol)

    Nphi1 = sum(vs1.weights .* vs1.density)
    Nphi2 = sum(vs2.weights .* vs2.density)

    return (
        vs1 = vs1,
        vs2 = vs2,
        vs_product = vs_product,
        Nphi1 = Nphi1,
        Nphi2 = Nphi2,
        shared_shift = sq.shared_shift,
        z_center = Float64(z_center),
    )
end
```

Note: the `using BoundaryIntegral` / `import BoundaryIntegral as BI` lines and the `include`/`using .MonolayerOrbitalLoader` must appear in the module body. Move them to the top of the module (just below `module MonolayerScreenedSolve`) when you implement, so all functions see them. The `export parse_coqui_loc` block grows to `export parse_coqui_loc, monolayer_screened_sources`.

- [ ] **Step 4: Run and confirm PASS**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerScreenedSolve.jl monolayer/test/test_monolayer_screened_solve.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): z-centered VS + product-source builder for screened solve

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Target-specs builders for density and Hund's channels

**Files:**
- Modify: `monolayer/src/MonolayerScreenedSolve.jl`
- Modify: `monolayer/test/test_monolayer_screened_solve.jl`

The bilayer solver's `solve_screened_mode` consumes a `target_specs` vector of NamedTuples with fields `pair, target_vs, Na, Nb`. We build two such lists.

- [ ] **Step 1: Write the failing test**

Append to `monolayer/test/test_monolayer_screened_solve.jl`:

```julia
@testset "density_target_specs and hund_target_specs" begin
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = 1e-3, z_center = 1.675,
    )
    density_specs = density_target_specs(src)
    @test length(density_specs) == 2
    @test density_specs[1].pair == "onsite"
    @test density_specs[2].pair == "nn"
    @test density_specs[1].Na ≈ src.Nphi1
    @test density_specs[1].Nb ≈ src.Nphi1
    @test density_specs[2].Na ≈ src.Nphi1
    @test density_specs[2].Nb ≈ src.Nphi2
    @test density_specs[1].target_vs === src.vs1
    @test density_specs[2].target_vs === src.vs2

    hund_specs = hund_target_specs(src)
    @test length(hund_specs) == 1
    @test hund_specs[1].pair == "hund"
    @test hund_specs[1].Na ≈ src.Nphi1
    @test hund_specs[1].Nb ≈ src.Nphi2
    @test hund_specs[1].target_vs === src.vs_product
end
```

- [ ] **Step 2: Run and confirm FAIL**

```bash
julia --project=. monolayer/test/runtests.jl
```

- [ ] **Step 3: Implement**

Append to `monolayer/src/MonolayerScreenedSolve.jl` (before `end # module`), and extend the export list to include `density_target_specs, hund_target_specs`:

```julia
export density_target_specs, hund_target_specs

"""
    density_target_specs(sources)

Build a 2-entry target_specs vector consumed by `solve_screened_mode` for the
density-source solve. Source is `sources.vs1` (squared orbital 1); targets are
`:onsite` → vs1 and `:nn` → vs2. Normalizations follow the bare-channel
Hund's-safe convention: divide by orbital norms (Nphi1, Nphi2), not by the
target-density integral.
"""
function density_target_specs(sources)
    return [
        (pair = "onsite", target_vs = sources.vs1, Na = sources.Nphi1, Nb = sources.Nphi1),
        (pair = "nn",     target_vs = sources.vs2, Na = sources.Nphi1, Nb = sources.Nphi2),
    ]
end

"""
    hund_target_specs(sources)

Build a 1-entry target_specs vector for the Hund's-source solve. Source and
target are both `sources.vs_product = phi1*phi2`. Normalizations follow the
same orbital-norm convention.
"""
function hund_target_specs(sources)
    return [
        (pair = "hund", target_vs = sources.vs_product, Na = sources.Nphi1, Nb = sources.Nphi2),
    ]
end
```

- [ ] **Step 4: Run and confirm PASS**

```bash
julia --project=. monolayer/test/runtests.jl
```

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerScreenedSolve.jl monolayer/test/test_monolayer_screened_solve.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): density and Hund's target_specs builders

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: `screened_run_record` flattener

**Files:**
- Modify: `monolayer/src/MonolayerScreenedSolve.jl`
- Modify: `monolayer/test/test_monolayer_screened_solve.jl`

Each `solve_screened_mode` returns a NamedTuple containing a `pair_results` vector. We flatten one pair_result into a CSV-ready row carrying the slab and mode metadata.

- [ ] **Step 1: Write the failing test** (pure data-massaging test using a hand-built fake input)

Append to `monolayer/test/test_monolayer_screened_solve.jl`:

```julia
@testset "screened_run_record packs pair_result + metadata" begin
    fake_pair_result = (
        pair = "onsite",
        target_vs = nothing,
        Na = 1.0, Nb = 1.0,
        u_int_raw = 0.1, u_scatter_raw = -0.02, u_total_raw = 0.08,
        u_int_ev = 17.4, u_scatter_ev = -3.5, u_total_ev = 13.9,
    )
    fake_solve_result = (
        residual = 1e-6,
        n_interface_points = 1234,
        tkm_kmax = 40.0,
        gmres_status = nothing,
        pair_results = [fake_pair_result],
    )
    row = screened_run_record(fake_pair_result, fake_solve_result;
        mode_label = "Sharp", bandwidth = missing,
        eps_in = 2.4, eps_out = 1.0, Lz = 3.35, L = 90.0, z_center = 1.675,
        source_tol = 1e-3, rhs_tol = 1e-3, lhs_tol = 1e-5,
        gmres_atol = 1e-5, gmres_rtol = 1e-5,
        n_quad = 6, edge_refine_level = 4, max_order = 64, max_depth = 12,
    )
    @test row.channel == "onsite"
    @test row.mode == "Sharp"
    @test ismissing(row.bandwidth)
    @test row.u_total_ev ≈ 13.9
    @test row.u_int_ev ≈ 17.4
    @test row.u_scatter_ev ≈ -3.5
    @test row.eps_in == 2.4
    @test row.Lz == 3.35
    @test row.n_interface_points == 1234
    @test row.sigma_residual ≈ 1e-6
end
```

- [ ] **Step 2: Run and confirm FAIL**

```bash
julia --project=. monolayer/test/runtests.jl
```

- [ ] **Step 3: Implement**

Append to `monolayer/src/MonolayerScreenedSolve.jl`, extend exports with `, screened_run_record`:

```julia
export screened_run_record

"""
    screened_run_record(pair_result, solve_result; mode_label, bandwidth, ...slab params...)

Flatten one `pair_result` from `solve_screened_mode` plus the surrounding solve
NamedTuple and the slab/mode parameters into a single NamedTuple ready for
DataFrame / CSV serialization.
"""
function screened_run_record(pair_result, solve_result;
    mode_label::AbstractString,
    bandwidth::Union{Real, Missing},
    eps_in::Real, eps_out::Real, Lz::Real, L::Real, z_center::Real,
    source_tol::Real, rhs_tol::Real, lhs_tol::Real,
    gmres_atol::Real, gmres_rtol::Real,
    n_quad::Integer, edge_refine_level::Integer,
    max_order::Integer, max_depth::Integer,
)
    return (
        channel = pair_result.pair,
        mode = String(mode_label),
        bandwidth = bandwidth,
        u_int_raw = pair_result.u_int_raw,
        u_scatter_raw = pair_result.u_scatter_raw,
        u_total_raw = pair_result.u_total_raw,
        u_int_ev = pair_result.u_int_ev,
        u_scatter_ev = pair_result.u_scatter_ev,
        u_total_ev = pair_result.u_total_ev,
        Na = pair_result.Na,
        Nb = pair_result.Nb,
        sigma_residual = solve_result.residual,
        n_interface_points = solve_result.n_interface_points,
        tkm_kmax = solve_result.tkm_kmax,
        eps_in = Float64(eps_in),
        eps_out = Float64(eps_out),
        Lz = Float64(Lz),
        L = Float64(L),
        z_center = Float64(z_center),
        source_tol = Float64(source_tol),
        rhs_tol = Float64(rhs_tol),
        lhs_tol = Float64(lhs_tol),
        gmres_atol = Float64(gmres_atol),
        gmres_rtol = Float64(gmres_rtol),
        n_quad = Int(n_quad),
        edge_refine_level = Int(edge_refine_level),
        max_order = Int(max_order),
        max_depth = Int(max_depth),
    )
end
```

- [ ] **Step 4: Run and confirm PASS**

```bash
julia --project=. monolayer/test/runtests.jl
```

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerScreenedSolve.jl monolayer/test/test_monolayer_screened_solve.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): screened_run_record CSV flattener

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Driver script — Sharp-mode smoke run on density solve

**Files:**
- Create: `monolayer/scripts/screened_monolayer_k323201.jl` (initial smoke version)

Before launching a full bandwidth sweep, we run one density solve at `Sharp` mode end-to-end. This validates the geometry, parser, normalization, and pipeline integration on a single ~5–30 min run.

- [ ] **Step 1: Create the script (smoke-only version)**

Create `monolayer/scripts/screened_monolayer_k323201.jl`:

```julia
using BoundaryIntegral
import BoundaryIntegral as BI

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve

include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1   = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2   = joinpath(REF_DIR, "graphene_00002.xsf")
const COQUI   = joinpath(REF_DIR, "_coqui_crpa_loc.out")

# Slab geometry
const LZ = 3.35
const Z_CENTER = LZ / 2
const EPS_IN = 2.4
const EPS_OUT = 1.0
const L = 90.0

# Solver tolerances (reused from the bilayer protocol)
const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12

# Density-channel bandwidth sweep (used by Task 7 main; Task 6 smoke uses only Sharp)
const BANDWIDTHS = [0.05, 0.1, 0.2, 0.4]

const OUT_CSV = joinpath(@__DIR__, "..", "data", "screened_monolayer_k323201.csv")

function _solve_kwargs()
    return (
        n_quad = N_QUAD,
        edge_refine_level = EDGE_REFINE_LEVEL,
        rhs_tol = RHS_TOL,
        lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL,
        gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER,
        max_depth = MAX_DEPTH,
        volume_tol = RHS_TOL,
    )
end

function _record_kwargs(mode_label, bandwidth)
    return (
        mode_label = mode_label, bandwidth = bandwidth,
        eps_in = EPS_IN, eps_out = EPS_OUT, Lz = LZ, L = L, z_center = Z_CENTER,
        source_tol = SOURCE_TOL, rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH,
    )
end

function smoke()
    println("Loading sources ...")
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER,
    )
    println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

    println("Running density Sharp solve ...")
    density_specs = density_target_specs(src)
    result = solve_screened_mode(
        src.vs1, density_specs,
        L, L, LZ, EPS_IN, EPS_OUT, BI.SharpScreening();
        _solve_kwargs()...,
    )
    println("  residual = ", result.residual,
            "  n_interface_points = ", result.n_interface_points)
    for pr in result.pair_results
        println("  ", pr.pair, "  u_total_ev = ", pr.u_total_ev,
                "  (u_int_ev = ", pr.u_int_ev, ", u_scatter_ev = ", pr.u_scatter_ev, ")")
    end
    coqui = parse_coqui_loc(COQUI)
    println("CoQui reference:")
    for ch in (:onsite, :nn)
        println("  ", ch, "  U_ijkl = ", coqui[ch].U_ijkl)
    end
    return result
end

# Run smoke when invoked directly (Task 6).
smoke()
```

- [ ] **Step 2: Run smoke (long compute, ~5–30 min)**

Same long-running pattern as Task 4 of the bare-channel plan: launch in background, poll until julia exits.

```bash
mkdir -p monolayer/logs
LOGFILE=monolayer/logs/screened_smoke_$(date +%Y%m%d-%H%M%S).log
echo "$LOGFILE" > /tmp/screened_smoke_logfile
julia --project=. monolayer/scripts/screened_monolayer_k323201.jl > "$LOGFILE" 2>&1 &
JULIA_PID=$!
echo "$JULIA_PID" > /tmp/screened_smoke_pid
echo "PID=$JULIA_PID  LOG=$LOGFILE"
```

Then poll (Bash with `timeout: 600000`):

```bash
JULIA_PID=$(cat /tmp/screened_smoke_pid)
LOGFILE=$(cat /tmp/screened_smoke_logfile)
while kill -0 "$JULIA_PID" 2>/dev/null; do
  sleep 20
  echo "  still running... tail:"
  tail -3 "$LOGFILE" | sed 's/^/    /'
done
wait "$JULIA_PID" 2>/dev/null
echo "exit code: $?"
cat "$LOGFILE"
```

If the 10-minute Bash timeout fires before julia finishes, immediately re-run the same loop (PID and LOGFILE persist in /tmp). Do not proceed until julia has exited.

Sanity criteria (record observed values; do not "fix" anything that disagrees with predictions):

- `Nphi1 ≈ Nphi2 ≈ 1.0` (within ~20%).
- `residual < 1e-3` and finite (GMRES converged).
- `onsite  u_total_ev` somewhere in `[7, 13]` eV (CoQui has 9.78 eV).
- `nn      u_total_ev` somewhere in `[3,  8]` eV (CoQui has 5.17 eV).
- `u_int_ev` should be close to our bare value at k_323201 (onsite 17.42, nn 8.84); the screening reduction is the `u_scatter_ev` term.

If `u_total_ev` is wildly out of range (e.g., negative, > 50, or NaN), STOP and report — likely a geometry (Z_CENTER) or normalization bug.

- [ ] **Step 3: Commit (script only, no CSV / no logs)**

```bash
git add monolayer/scripts/screened_monolayer_k323201.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): screened-monolayer Sharp-mode smoke driver

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Driver script — full sweep + CSV

**Files:**
- Modify: `monolayer/scripts/screened_monolayer_k323201.jl`

Replace the smoke-only invocation with a full sweep: density solve over Sharp + 4 soft bandwidths, Hund's solve at Sharp only.

- [ ] **Step 1: Replace `smoke()` and its invocation with `main()`**

Replace the entire `function smoke() ... end\nsmoke()` block at the bottom of `monolayer/scripts/screened_monolayer_k323201.jl` with:

```julia
function _density_modes()
    cases = NamedTuple{(:label, :bandwidth, :mode), Tuple{String, Union{Missing, Float64}, BI.AbstractScreeningMode}}[]
    push!(cases, (label = "Sharp", bandwidth = missing, mode = BI.SharpScreening()))
    for bw in BANDWIDTHS
        push!(cases, (label = "SoftMixInversePermittivity", bandwidth = bw, mode = BI.SoftMixInversePermittivity(bw)))
    end
    return cases
end

function _hund_modes()
    return [
        (label = "Sharp", bandwidth = missing, mode = BI.SharpScreening()),
    ]
end

function main()
    println("Loading sources ...")
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER,
    )
    println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

    coqui = parse_coqui_loc(COQUI)

    rows = NamedTuple[]

    density_specs = density_target_specs(src)
    for case in _density_modes()
        println("Density solve  mode=", case.label, "  bandwidth=", case.bandwidth)
        result = solve_screened_mode(
            src.vs1, density_specs,
            L, L, LZ, EPS_IN, EPS_OUT, case.mode;
            _solve_kwargs()...,
        )
        println("  residual = ", result.residual,
                "  n_interface_points = ", result.n_interface_points)
        for pr in result.pair_results
            ch_sym = pr.pair == "onsite" ? :onsite : :nn
            row0 = screened_run_record(pr, result; _record_kwargs(case.label, case.bandwidth)...)
            row = merge(row0, (
                coqui_v_ijkl = coqui[ch_sym].v_ijkl,
                coqui_U_ijkl = coqui[ch_sym].U_ijkl,
                diff_ev = coqui[ch_sym].U_ijkl - row0.u_total_ev,
            ))
            push!(rows, row)
            println("    ", pr.pair, "  ours=", row0.u_total_ev,
                    "  CoQui=", coqui[ch_sym].U_ijkl,
                    "  diff=", row.diff_ev)
        end
    end

    hund_specs = hund_target_specs(src)
    for case in _hund_modes()
        println("Hund's solve  mode=", case.label, "  bandwidth=", case.bandwidth)
        result = solve_screened_mode(
            src.vs_product, hund_specs,
            L, L, LZ, EPS_IN, EPS_OUT, case.mode;
            _solve_kwargs()...,
        )
        println("  residual = ", result.residual,
                "  n_interface_points = ", result.n_interface_points)
        for pr in result.pair_results
            row0 = screened_run_record(pr, result; _record_kwargs(case.label, case.bandwidth)...)
            row = merge(row0, (
                coqui_v_ijkl = coqui[:hund_sf].v_ijkl,
                coqui_U_ijkl = coqui[:hund_sf].U_ijkl,
                diff_ev = coqui[:hund_sf].U_ijkl - row0.u_total_ev,
            ))
            push!(rows, row)
            println("    ", pr.pair, "  ours=", row0.u_total_ev,
                    "  CoQui=", coqui[:hund_sf].U_ijkl,
                    "  diff=", row.diff_ev)
        end
    end

    mkpath(dirname(OUT_CSV))
    table = DataFrame(rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    summary_cols = [:channel, :mode, :bandwidth, :u_total_ev, :coqui_U_ijkl, :diff_ev, :sigma_residual]
    show(table[:, summary_cols]; allrows = true, allcols = true)
    println()
end

main()
```

- [ ] **Step 2: Run the full sweep (long compute, ~30 min – 3 hours)**

Same background-poll protocol:

```bash
mkdir -p monolayer/logs
LOGFILE=monolayer/logs/screened_sweep_$(date +%Y%m%d-%H%M%S).log
echo "$LOGFILE" > /tmp/screened_sweep_logfile
julia --project=. monolayer/scripts/screened_monolayer_k323201.jl > "$LOGFILE" 2>&1 &
JULIA_PID=$!
echo "$JULIA_PID" > /tmp/screened_sweep_pid
echo "PID=$JULIA_PID  LOG=$LOGFILE"
```

Poll (re-run as needed with 600000 ms timeout):

```bash
JULIA_PID=$(cat /tmp/screened_sweep_pid)
LOGFILE=$(cat /tmp/screened_sweep_logfile)
while kill -0 "$JULIA_PID" 2>/dev/null; do
  sleep 30
  echo "  still running... tail:"
  tail -3 "$LOGFILE" | sed 's/^/    /'
done
wait "$JULIA_PID" 2>/dev/null
echo "exit code: $?"
tail -50 "$LOGFILE"
```

Sanity criteria:

1. CSV has 11 rows (2 channels × 5 modes + 1 Hund's × 1 mode).
2. Every row has `sigma_residual < 1e-3` (GMRES converged).
3. `onsite u_total_ev` and `nn u_total_ev` are both in their respective expected bands across modes.
4. Hund's row's `u_total_ev` is order 0.05–0.15 eV (CoQui has 0.0937).

- [ ] **Step 3: Inspect CSV**

```bash
column -t -s, monolayer/data/screened_monolayer_k323201.csv | cut -c1-160 | head -15
```

- [ ] **Step 4: Commit script + CSV**

```bash
git add monolayer/scripts/screened_monolayer_k323201.jl monolayer/data/screened_monolayer_k323201.csv
git commit -m "$(cat <<'EOF'
feat(monolayer): full screened k_323201 sweep + CSV (5 density modes + Hund's Sharp)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Report

**Files:**
- Create: `monolayer/results/2026-05-13-screened-monolayer-k323201-report.md`

- [ ] **Step 1: Pull values from the CSV**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_k323201.csv", DataFrame)
show(df[:, [:channel, :mode, :bandwidth, :u_int_ev, :u_scatter_ev, :u_total_ev, :coqui_U_ijkl, :diff_ev, :sigma_residual]]; allrows=true, allcols=true)
println()
'
```

- [ ] **Step 2: Create the report**

Create `monolayer/results/2026-05-13-screened-monolayer-k323201-report.md`. Replace every `<fill>` with the actual CSV value; compute the verdict per the spec criteria (Strong: `|rel_err| ≤ 5%` for onsite and nn; Acceptable: `≤ 15%`; Disagreement otherwise).

```markdown
# Screened Monolayer Graphene at k_323201 — Report

Date: 2026-05-13

## Context

The 2026-05-13 bare-monolayer k-mesh sweep showed our direct Wannier integration and CoQui's bare V agree to ~10 mV across all 4 channels at k_323201, establishing this dataset as the joint validation point. This report extends the comparison to the screened (cRPA) `U_ijkl` table using a slab dielectric model.

Reference: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out`.

| channel  | v_ijkl (eV) | U_ijkl (eV) |
|----------|------------:|------------:|
| onsite   | 17.4029     | 9.7824      |
| nn       |  8.8319     | 5.1678      |
| hund_ph  |  0.1316     | 0.0937      |
| hund_sf  |  0.1316     | 0.0937      |

## Method

Slab dielectric model with two GMRES dielectric solves:

| parameter | value |
|-----------|------:|
| LZ        | 3.35 Å |
| eps_in    | 2.4    |
| eps_out   | 1.0    |
| Z_CENTER  | 1.675 Å (LZ/2, symmetric) |
| L (in-plane box) | 90 Å |
| source_tol | 1e-3 |
| rhs_tol   | 1e-3  |
| lhs_tol   | 1e-5  |
| GMRES atol/rtol | 1e-5 |

- Density-source solve (source = `|phi_1|^2`) evaluated at targets `|phi_1|^2` (→ onsite) and `|phi_2|^2` (→ nn), swept over `Sharp` + 4 soft-mix-inverse-permittivity bandwidths `[0.05, 0.1, 0.2, 0.4]`.
- Hund's-source solve (source = target = `phi_1·phi_2` signed product), `Sharp` mode only.
- eV conversion uses orbital-norm normalization (`Na, Nb = Nphi1, Nphi2`), consistent with the bare-channel Hund's-safe convention.

## Results

### Density channels (5 modes)

| mode               | bandwidth | onsite ours | onsite CoQui | onsite diff | nn ours | nn CoQui | nn diff |
|--------------------|----------:|------------:|-------------:|------------:|--------:|---------:|--------:|
| Sharp              |       —   | <fill>      | 9.7824       | <fill>      | <fill>  | 5.1678   | <fill>  |
| SoftMix bw=0.05    |     0.05  | <fill>      | 9.7824       | <fill>      | <fill>  | 5.1678   | <fill>  |
| SoftMix bw=0.10    |     0.10  | <fill>      | 9.7824       | <fill>      | <fill>  | 5.1678   | <fill>  |
| SoftMix bw=0.20    |     0.20  | <fill>      | 9.7824       | <fill>      | <fill>  | 5.1678   | <fill>  |
| SoftMix bw=0.40    |     0.40  | <fill>      | 9.7824       | <fill>      | <fill>  | 5.1678   | <fill>  |

### Hund's channel (Sharp only)

| mode  | hund_total_ev | CoQui U_ijkl | diff_ev |
|-------|--------------:|-------------:|--------:|
| Sharp | <fill>        | 0.0937       | <fill>  |

### Decomposition into direct vs scattered (Sharp mode)

| channel | u_int_ev (≈ bare) | u_scatter_ev (screening) | u_total_ev |
|---------|------------------:|-------------------------:|-----------:|
| onsite  | <fill>            | <fill>                   | <fill>     |
| nn      | <fill>            | <fill>                   | <fill>     |
| hund    | <fill>            | <fill>                   | <fill>     |

The `u_int_ev` should agree with the bare-channel CSV's `our_u_ev` at k_323201 (onsite 17.42, nn 8.84, hund 0.13). Significant disagreement here means the slab geometry or the source construction is wrong.

## Verdict

Spec thresholds (relative error vs CoQui U_ijkl):
- Strong agreement: `|rel_err| ≤ 5%` for onsite and nn
- Acceptable:       `|rel_err| ≤ 15%` for onsite and nn
- Disagreement:     otherwise

**Verdict: <fill>** — choose Strong / Acceptable / Disagreement based on the best-performing density mode. State which mode/bandwidth gives the cleanest match, and one sentence about whether the Hund's channel falls in the expected range.

## Reproducibility

- Driver: `monolayer/scripts/screened_monolayer_k323201.jl`
- CSV:    `monolayer/data/screened_monolayer_k323201.csv`
- Run log: `monolayer/logs/screened_sweep_*.log`

## Notes

- All GMRES residuals: `<fill: smallest>` to `<fill: largest>` (target < 1e-3).
- `n_interface_points` for the density vs Hund's solves: `<fill density>` / `<fill hund>`.
- Hund's channel: only `Sharp` mode in this pass. Soft modes deferred per the design.

## Next

Recommended follow-ups in order:

1. <fill: brief one-line judgment based on residuals> — either "current parameters are good enough; move to k-mesh sweep on screened side" or "sweep `eps_in` between 1.5 and 3.0 to find the empirical best-fit value" or "swap to a different slab-thickness value (`LZ` between 3.0 and 4.0)" depending on the verdict.
2. Replicate at k_252501 / k_161601 to see whether the same parameter set works across k-meshes (or whether the screened comparison is also k-mesh-sensitive in the way the bare comparison was for CoQui's side).
3. If Hund's `Sharp` falls outside `[0.05, 0.15]` eV, add the soft-mix sweep for the Hund's channel as a follow-up plan.
```

- [ ] **Step 3: Fill in every `<fill>` from the CSV and choose the verdict**

Open the file in editor, edit each cell.

- [ ] **Step 4: Commit**

```bash
git add monolayer/results/2026-05-13-screened-monolayer-k323201-report.md
git commit -m "$(cat <<'EOF'
docs(monolayer): screened k_323201 results + verdict

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Final verification

- [ ] **Step 1: Working tree clean, expected commits present**

```bash
git status
git log --oneline -12
```

Expected commits (newest first):

- `docs(monolayer): screened k_323201 results + verdict`
- `feat(monolayer): full screened k_323201 sweep + CSV (5 density modes + Hund's Sharp)`
- `feat(monolayer): screened-monolayer Sharp-mode smoke driver`
- `feat(monolayer): screened_run_record CSV flattener`
- `feat(monolayer): density and Hund's target_specs builders`
- `feat(monolayer): z-centered VS + product-source builder for screened solve`
- `feat(monolayer): CoQui loc.out parser for U_ijkl channels`
- `spec: screened monolayer k_323201 ...`

- [ ] **Step 2: Run full test suite**

```bash
julia --project=. test/runtests.jl
```

Expected: all bilayer_slab + monolayer tests pass, including the new `MonolayerScreenedSolve` testset.

- [ ] **Step 3: CSV sanity check**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_k323201.csv", DataFrame)
@assert nrow(df) == 11 "expected 11 rows, got $(nrow(df))"
@assert all(df.sigma_residual .< 1e-3) "some GMRES residual not converged"
println("OK: ", nrow(df), " rows, residuals all under 1e-3")
'
```

---

## Self-Review

**Spec coverage:**
- Spec § "Scope → In scope" — 2 GMRES solves, 4 channels, 5+1 modes, CoQui parsing, dated report: Tasks 2 (parser), 3 (sources), 4 (target_specs), 5 (record flattener), 6 (smoke), 7 (sweep), 8 (report).
- Spec § "File layout" — every new file is created by an explicit task.
- Spec § "Slab parameters" — encoded as `const` declarations in Task 6 Step 1.
- Spec § "Hund's source treatment" — covered by `monolayer_screened_sources` (Task 3) producing `vs_product` and by `hund_target_specs` (Task 4) using the orbital-norm convention.
- Spec § "Verdict criteria" — encoded in the report template (Task 8) and the Self-Review block of the report.

**Placeholder scan:**
- Report template has `<fill>` markers in the cells the implementer must populate from CSV. These are intentional and Task 8 Step 3 is the populate step. No `<fill>` outside the report.
- No `TBD` / `TODO` / "appropriate" anywhere in the plan.

**Type consistency:**
- `parse_coqui_loc` returns `Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}` — consumed in Task 7 main via `coqui[ch_sym].v_ijkl` and `.U_ijkl`. ✓
- `monolayer_screened_sources` returns a NamedTuple with `vs1, vs2, vs_product, Nphi1, Nphi2, shared_shift, z_center` — consumed in Task 4 (`density_target_specs(src)`, `hund_target_specs(src)`) and Task 6/7 (`src.vs1`, `src.vs_product`, `src.Nphi1`, etc). ✓
- `target_specs` returned by `density_target_specs` / `hund_target_specs` use the exact field names `pair, target_vs, Na, Nb` consumed by the bilayer solver (verified at planning time from `bilayer_slab/src/ScreenedOrbitalSolve.jl` line 305). ✓
- `screened_run_record` returns a NamedTuple with `channel, mode, bandwidth, u_int_raw, u_scatter_raw, u_total_raw, u_int_ev, u_scatter_ev, u_total_ev, Na, Nb, sigma_residual, n_interface_points, tkm_kmax, eps_in, eps_out, Lz, L, z_center, source_tol, rhs_tol, lhs_tol, gmres_atol, gmres_rtol, n_quad, edge_refine_level, max_order, max_depth` — merged with `(coqui_v_ijkl, coqui_U_ijkl, diff_ev)` in Task 7. The CSV column set is the union; consistent. ✓
- The `_solve_kwargs()` helper in Task 6 returns a NamedTuple consumed via `_solve_kwargs()...` splatting. The keys exactly match `solve_screened_mode`'s keyword arguments. ✓
- The `_record_kwargs(label, bw)` helper passes `mode_label, bandwidth` plus all 11 slab/tol parameters to `screened_run_record`; matches `screened_run_record`'s signature. ✓

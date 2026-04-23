# Monolayer Graphene Revision Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Archive the existing bilayer-slab workflow under `bilayer_slab/`, add a new parallel `monolayer/` workflow that computes the four bare Coulomb matrix elements (onsite, nearest-neighbour, two Hund's channels) from Malte Rösner's new monolayer Wannier XSF orbitals, and compare against his CoQui cRPA reference values.

**Architecture:** Shared top-level `Project.toml` environment. Existing bilayer code moves verbatim into `bilayer_slab/`; new monolayer code lives in `monolayer/` and uses the same `BoundaryIntegral` (`BI`) volume-integral primitives. The monolayer bare workflow reads signed Wannier orbitals directly from XSF (no squaring at load time), builds channel-specific `(f_src, f_tgt)` scalar-field pairs, integrates `1/|r-r'|` with `BI.TKM3D.ltkm3dc` on `VolumeSource` objects, and writes a CSV + markdown report. No dielectric solve, no GMRES. A mirror-pad convergence check provides a single diagnostic on xy tail convergence.

**Tech Stack:** Julia, `BoundaryIntegral.jl`, `CSV.jl`, `DataFrames.jl`, `Test.jl`. Reference XSF source path: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15/graphene_0000{1,2}.xsf`.

---

## File Structure

Paths are relative to `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/`.

**Moved (bilayer archive):**
- `src/*.jl` → `bilayer_slab/src/*.jl`
- `scripts/*.jl` → `bilayer_slab/scripts/*.jl`
- `data/*` → `bilayer_slab/data/*`
- `results/*` → `bilayer_slab/results/*`
- `test/test_screened_density_analysis.jl`, `test/test_screened_orbital_solve.jl`, `test/test_screened_hubbard_comparison.jl` → `bilayer_slab/test/`
- New `bilayer_slab/test/runtests.jl` that `include`s the three moved test files.

**New (monolayer workflow):**
- `monolayer/src/MonolayerOrbitalLoader.jl` — signed / squared XSF loader with shared centering shift + optional mirror-pad.
- `monolayer/src/MonolayerBareIntegrals.jl` — pair-source construction for the four channels + direct `1/|r-r'|` volume integration.
- `monolayer/scripts/bare_monolayer_graphene.jl` — run all four channels, write CSV, print summary.
- `monolayer/data/bare_monolayer_graphene.csv` — produced by the script (output artifact).
- `monolayer/results/2026-04-23-bare-monolayer-report.md` — written at the end of the plan.
- `monolayer/test/runtests.jl`, `monolayer/test/test_monolayer_orbital_loader.jl`, `monolayer/test/test_monolayer_bare_integrals.jl`.

**Shared (top level):**
- `test/runtests.jl` — rewritten as a thin dispatcher that includes both sub-suites.
- `docs/superpowers/specs/2026-04-23-monolayer-graphene-revision-design.md` — exists (spec).
- `results/2026-03-26-screened-graphene-report.md` is moved to `bilayer_slab/results/` and gets a historical header.

---

## Task 1: Reorganise repo — move bilayer workflow under `bilayer_slab/`

**Files:**
- Move: `src/`, `scripts/`, `data/`, `results/` → `bilayer_slab/src/`, `bilayer_slab/scripts/`, `bilayer_slab/data/`, `bilayer_slab/results/`
- Move: `test/test_screened_*.jl` → `bilayer_slab/test/`
- Create: `bilayer_slab/test/runtests.jl`
- Rewrite: `test/runtests.jl` as dispatcher

- [ ] **Step 1: Create the target directories**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
mkdir -p bilayer_slab/src bilayer_slab/scripts bilayer_slab/data bilayer_slab/results bilayer_slab/test
```

- [ ] **Step 2: Move the bilayer files with `git mv`**

```bash
git mv src/ScreenedDensityAnalysis.jl     bilayer_slab/src/ScreenedDensityAnalysis.jl
git mv src/ScreenedHubbardComparison.jl   bilayer_slab/src/ScreenedHubbardComparison.jl
git mv src/ScreenedOrbitalSolve.jl        bilayer_slab/src/ScreenedOrbitalSolve.jl
rmdir src

git mv scripts/plot_bilayer_hubbard_vs_reference.jl bilayer_slab/scripts/
git mv scripts/plot_screened_hubbard_vs_crpa.jl     bilayer_slab/scripts/
git mv scripts/screened_density_kz_decay_z_upsampled.jl bilayer_slab/scripts/
git mv scripts/screened_density_tkm.jl              bilayer_slab/scripts/
git mv scripts/screened_hubbard_graphene.jl         bilayer_slab/scripts/
rmdir scripts

git mv data/PhysRevLett.106.236805.pdf                bilayer_slab/data/
git mv data/screened_hubbard_graphene.csv             bilayer_slab/data/
git mv data/screened_hubbard_graphene_sharp.csv       bilayer_slab/data/
git mv data/screened_hubbard_graphene_soft.csv        bilayer_slab/data/
rmdir data

git mv results/2026-03-26-screened-graphene-report.md bilayer_slab/results/
rmdir results

git mv test/test_screened_density_analysis.jl    bilayer_slab/test/
git mv test/test_screened_orbital_solve.jl      bilayer_slab/test/
git mv test/test_screened_hubbard_comparison.jl bilayer_slab/test/
```

- [ ] **Step 3: Create `bilayer_slab/test/runtests.jl` matching the old contents**

```julia
using Test

include("test_screened_density_analysis.jl")
include("test_screened_orbital_solve.jl")
include("test_screened_hubbard_comparison.jl")
```

- [ ] **Step 4: Rewrite top-level `test/runtests.jl` as a dispatcher**

```julia
using Test

@testset "bilayer_slab" begin
    include(joinpath(@__DIR__, "..", "bilayer_slab", "test", "runtests.jl"))
end

@testset "monolayer" begin
    monolayer_runtests = joinpath(@__DIR__, "..", "monolayer", "test", "runtests.jl")
    if isfile(monolayer_runtests)
        include(monolayer_runtests)
    else
        @info "monolayer test suite not yet present; skipping"
    end
end
```

- [ ] **Step 5: Update path references inside moved bilayer files**

The moved files may contain path literals that were relative to the old `graphene/` root. For each moved `.jl` file, open it and repair any relative paths:

- `bilayer_slab/scripts/screened_hubbard_graphene.jl` uses `joinpath(@__DIR__, "..", "src", "ScreenedOrbitalSolve.jl")` — this remains correct after the move because `@__DIR__` for the moved script is now `bilayer_slab/scripts/`, and `../src/ScreenedOrbitalSolve.jl` correctly resolves to `bilayer_slab/src/ScreenedOrbitalSolve.jl`. **No change needed.**
- Same pattern: paths of the form `joinpath(@__DIR__, "..", "data", ...)` or `joinpath(@__DIR__, "..", "results", ...)` remain correct because the sibling directories moved together.
- `DEFAULT_ORBITAL_1` / `DEFAULT_ORBITAL_2` in `bilayer_slab/src/ScreenedOrbitalSolve.jl` use `joinpath(@__DIR__, "..", "..", "..", "density_data", ...)`. After the move, `@__DIR__` is `bilayer_slab/src`, so four levels up is `codes/`; but the original intent was `four_index_integral_solver/density_data/`. The original had three `..` levels from `graphene/src` → `codes/graphene` → `codes` → `four_index_integral_solver`. After the move we need **four** `..` levels: `bilayer_slab/src` → `graphene/bilayer_slab` → `graphene` → `codes` → `four_index_integral_solver`. Change the `..` count from 3 to 4:

```julia
const DEFAULT_ORBITAL_1 = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "density_data", "graphene_00001_5x5x1_shifted.xsf"))
const DEFAULT_ORBITAL_2 = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "density_data", "graphene_00002_5x5x1_shifted.xsf"))
```

Verify the new path resolves:

```bash
julia --project=. -e 'path = normpath(joinpath(pwd(), "bilayer_slab", "src", "..", "..", "..", "..", "density_data", "graphene_00001_5x5x1_shifted.xsf")); println(path); println(isfile(path))'
```

Expected: absolute path ending in `density_data/graphene_00001_5x5x1_shifted.xsf`, `true`.

- [ ] **Step 6: Verify the existing bilayer suite still passes**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
julia --project=. test/runtests.jl
```

Expected: all bilayer test sets pass; monolayer suite logs "not yet present; skipping".

- [ ] **Step 7: Add a historical header to the bilayer report**

Edit `bilayer_slab/results/2026-03-26-screened-graphene-report.md` — insert at the top, immediately after the `# Screened Graphene Results Report` heading:

```markdown
> **Historical note (2026-04-23):** This report documents the bilayer-slab dielectric approximation with scalar `eps_in/eps_out`. It has been superseded for the purpose of reference validation by the monolayer workflow under `graphene/monolayer/`, which compares directly against Malte Rösner's CoQui cRPA values. The numbers below remain valid within the slab model but are not a cRPA reproduction.
```

- [ ] **Step 8: Commit**

```bash
git add -A
git commit -m "chore: move bilayer-slab workflow under bilayer_slab/ and add test dispatcher"
```

---

## Task 2: Minimal loader — read signed XSF orbital and confirm grid metadata

**Files:**
- Create: `monolayer/src/MonolayerOrbitalLoader.jl`
- Create: `monolayer/test/test_monolayer_orbital_loader.jl`
- Create: `monolayer/test/runtests.jl`

**Reference XSF path:** `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15/graphene_00001.xsf` (and `..._00002.xsf`). Grid shape `150×150×192`, lattice `~12.24 × 10.60 × 14.92 Å`.

- [ ] **Step 1: Write the failing test**

Create `monolayer/test/test_monolayer_orbital_loader.jl`:

```julia
using Test

include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")

@testset "MonolayerOrbitalLoader.load_signed_xsf" begin
    @test isfile(XSF_1)
    datagrid = load_signed_xsf(XSF_1)
    @test size(datagrid.values) == (150, 150, 192)
    # Wannier orbitals are signed (not a probability density): expect both signs present.
    @test minimum(datagrid.values) < 0
    @test maximum(datagrid.values) > 0
end
```

Create `monolayer/test/runtests.jl`:

```julia
using Test

include("test_monolayer_orbital_loader.jl")
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
julia --project=. monolayer/test/runtests.jl
```

Expected: FAIL — `MonolayerOrbitalLoader` module does not exist.

- [ ] **Step 3: Implement the minimal loader**

Create `monolayer/src/MonolayerOrbitalLoader.jl`:

```julia
module MonolayerOrbitalLoader

using BoundaryIntegral
import BoundaryIntegral as BI

export load_signed_xsf, load_squared_xsf

"""
    load_signed_xsf(path)

Read a Wannier90 XSF file and return the raw signed datagrid (a NamedTuple
compatible with `BI.VolumeSource`). Values are not squared.
"""
function load_signed_xsf(path::AbstractString)
    _, datagrid = BI.read_xsf(path)
    return datagrid
end

"""
    load_squared_xsf(path)

Read a Wannier90 XSF file and return `|phi|^2` as a datagrid.
"""
function load_squared_xsf(path::AbstractString)
    _, raw = BI.read_xsf(path)
    values = copy(raw.values)
    values .*= values
    return merge(raw, (; values = values))
end

end # module
```

- [ ] **Step 4: Run test to verify it passes**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerOrbitalLoader.jl monolayer/test/runtests.jl monolayer/test/test_monolayer_orbital_loader.jl
git commit -m "feat(monolayer): add signed/squared XSF loader"
```

---

## Task 3: Shared centering shift applied to both orbitals

**Files:**
- Modify: `monolayer/src/MonolayerOrbitalLoader.jl`
- Modify: `monolayer/test/test_monolayer_orbital_loader.jl`

The bilayer workflow's invariant is: compute the centroid of orbital 1 once, apply the same shift to orbital 2. We preserve that here.

- [ ] **Step 1: Extend the failing test**

Append to `monolayer/test/test_monolayer_orbital_loader.jl`:

```julia
@testset "MonolayerOrbitalLoader.centered_monolayer_sources — shared shift" begin
    out = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = true)

    @test hasproperty(out, :vs1)
    @test hasproperty(out, :vs2)
    @test hasproperty(out, :shared_shift)
    @test length(out.shared_shift) == 3

    # The shift is defined to drive orbital 1's density-weighted centroid to the origin.
    weights_1 = out.vs1.weights .* out.vs1.density
    total_1 = sum(weights_1)
    centroid_1 = ntuple(d -> sum(out.vs1.positions[d, :] .* weights_1) / total_1, 3)
    @test all(abs.(centroid_1) .< 1e-8)

    # Orbital 2 gets the same shift, so its centroid sits at the A→B bond vector (~1.4232 Å in y).
    weights_2 = out.vs2.weights .* out.vs2.density
    total_2 = sum(weights_2)
    centroid_2 = ntuple(d -> sum(out.vs2.positions[d, :] .* weights_2) / total_2, 3)
    @test abs(centroid_2[1]) < 1e-2
    @test abs(centroid_2[2] - 1.4231684) < 5e-2  # within 50 mÅ of atom-to-atom offset
    @test abs(centroid_2[3]) < 1e-2
end
```

- [ ] **Step 2: Run the test and expect failure**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: FAIL — `centered_monolayer_sources` not defined.

- [ ] **Step 3: Implement the shared-shift loader**

Append to `monolayer/src/MonolayerOrbitalLoader.jl` before `end # module`:

```julia
export centered_monolayer_sources, shift_datagrid

function shift_datagrid(datagrid, shift::NTuple{3, <:Real})
    origin = ntuple(i -> Float64(datagrid.origin[i]) + Float64(shift[i]), 3)
    return merge(datagrid, (; origin = origin))
end

function _density_centroid(vs::BI.VolumeSource)
    weights = vs.weights .* vs.density
    total = sum(weights)
    iszero(total) && throw(ArgumentError("zero total density"))
    return ntuple(d -> sum(vs.positions[d, :] .* weights) / total, 3)
end

function _centering_shift(datagrid; tol::Real)
    vs = BI.VolumeSource(datagrid, tol = tol)
    centroid = _density_centroid(vs)
    return ntuple(i -> -Float64(centroid[i]), 3)
end

"""
    centered_monolayer_sources(; orbital_1, orbital_2, source_tol, square)

Load two monolayer Wannier XSF orbitals, compute a centering shift from orbital 1
(using squared density to locate the orbital core), apply the same shift to both,
and return `VolumeSource` objects. When `square = true`, both sources use |phi|^2;
when `square = false`, both use the signed phi directly.
"""
function centered_monolayer_sources(;
    orbital_1::AbstractString,
    orbital_2::AbstractString,
    source_tol::Real = 1e-3,
    square::Bool = true,
)
    centering_datagrid_1 = load_squared_xsf(orbital_1)
    shared_shift = _centering_shift(centering_datagrid_1; tol = source_tol)

    datagrid_1 = square ? load_squared_xsf(orbital_1) : load_signed_xsf(orbital_1)
    datagrid_2 = square ? load_squared_xsf(orbital_2) : load_signed_xsf(orbital_2)

    datagrid_1_shifted = shift_datagrid(datagrid_1, shared_shift)
    datagrid_2_shifted = shift_datagrid(datagrid_2, shared_shift)

    vs1 = BI.VolumeSource(datagrid_1_shifted, tol = source_tol)
    vs2 = BI.VolumeSource(datagrid_2_shifted, tol = source_tol)

    return (
        datagrid_1 = datagrid_1_shifted,
        datagrid_2 = datagrid_2_shifted,
        vs1 = vs1,
        vs2 = vs2,
        shared_shift = shared_shift,
    )
end
```

- [ ] **Step 4: Run the test and confirm pass**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerOrbitalLoader.jl monolayer/test/test_monolayer_orbital_loader.jl
git commit -m "feat(monolayer): shared centering shift loader"
```

---

## Task 4: Pair-source construction for four channels

**Files:**
- Create: `monolayer/src/MonolayerBareIntegrals.jl`
- Create: `monolayer/test/test_monolayer_bare_integrals.jl`
- Modify: `monolayer/test/runtests.jl`

Channel definitions (see spec):
- `:onsite` — `(f_src, f_tgt) = (|phi1|^2, |phi1|^2)` (use orbital 1 for both)
- `:nn` — `(|phi1|^2, |phi2|^2)`
- `:hund_sf` — `(phi1·phi2, phi2·phi1) = (phi1·phi2, phi1·phi2)` (signed product)
- `:hund_ph` — `(phi1·phi2, phi1·phi2)` (same as spin-flip for real orbitals — we still code both independently as two separate channel labels)

We intentionally implement all four on top of the loader's `VolumeSource` outputs. The product `phi1·phi2` is built from the signed datagrids, then passed through `BI.VolumeSource` just like a density.

- [ ] **Step 1: Write the failing test for pair-source shapes and orbital-1 self-pair**

Create `monolayer/test/test_monolayer_bare_integrals.jl`:

```julia
using Test
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
include(joinpath(@__DIR__, "..", "src", "MonolayerBareIntegrals.jl"))
using .MonolayerOrbitalLoader
using .MonolayerBareIntegrals

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")

@testset "channel_pair_sources — shape and norms" begin
    densities = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = true)
    signed    = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = false)

    for channel in (:onsite, :nn, :hund_sf, :hund_ph)
        pair = channel_pair_sources(densities, signed, channel)
        @test hasproperty(pair, :source)
        @test hasproperty(pair, :target)
        # Length agreement
        @test size(pair.source.positions, 1) == 3
        @test size(pair.target.positions, 1) == 3
    end

    # Density-density norms should be the orbital normalization (~1 for a Wannier function).
    onsite = channel_pair_sources(densities, signed, :onsite)
    n_on = sum(onsite.source.weights .* onsite.source.density)
    @test 0.8 < n_on < 1.2

    # Hund's signed-product pair integrates to ≈ 0 because phi1 and phi2 are orthogonal Wannier orbitals.
    hund = channel_pair_sources(densities, signed, :hund_sf)
    n_hund = sum(hund.source.weights .* hund.source.density)
    @test abs(n_hund) < 5e-2
end
```

Add to `monolayer/test/runtests.jl` (append):

```julia
include("test_monolayer_bare_integrals.jl")
```

- [ ] **Step 2: Run tests — expect failure**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: FAIL — `MonolayerBareIntegrals` missing.

- [ ] **Step 3: Implement pair-source construction**

Create `monolayer/src/MonolayerBareIntegrals.jl`:

```julia
module MonolayerBareIntegrals

using BoundaryIntegral
using LinearAlgebra
import BoundaryIntegral as BI

export channel_pair_sources, bare_channel_integral, compute_all_bare_channels, E2_4PIEPS0, to_eV

const E2_4PIEPS0 = 14.3996

to_eV(raw::Real, Na::Real, Nb::Real) = raw * 4π * E2_4PIEPS0 / (Na * Nb)

"""
Build `(source, target)` VolumeSource pair for a given bare-interaction channel.

- `densities` carries |phi1|^2 and |phi2|^2 (`square = true` loader output).
- `signed`    carries  phi1 and  phi2 (`square = false` loader output).
Both must have been produced with the **same** centering shift and source_tol.

Channel semantics:
- `:onsite`    — |phi1|^2 vs |phi1|^2
- `:nn`        — |phi1|^2 vs |phi2|^2
- `:hund_sf`   — (phi1*phi2) vs (phi1*phi2)
- `:hund_ph`   — (phi1*phi2) vs (phi1*phi2)  (same numerics as :hund_sf for real orbitals)
"""
function channel_pair_sources(densities, signed, channel::Symbol)
    if channel === :onsite
        return (source = densities.vs1, target = densities.vs1)
    elseif channel === :nn
        return (source = densities.vs1, target = densities.vs2)
    elseif channel === :hund_sf || channel === :hund_ph
        product_datagrid = _product_datagrid(signed.datagrid_1, signed.datagrid_2)
        vs_product = BI.VolumeSource(product_datagrid, tol = 1e-3)
        return (source = vs_product, target = vs_product)
    else
        throw(ArgumentError("unknown channel $channel"))
    end
end

function _product_datagrid(a, b)
    size(a.values) == size(b.values) || throw(ArgumentError("grid shapes differ"))
    a.origin == b.origin || throw(ArgumentError("grid origins differ"))
    values = a.values .* b.values
    return merge(a, (; values = values))
end

end # module
```

- [ ] **Step 4: Run the test — expect PASS**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerBareIntegrals.jl monolayer/test/test_monolayer_bare_integrals.jl monolayer/test/runtests.jl
git commit -m "feat(monolayer): pair-source construction for 4 bare channels"
```

---

## Task 5: Direct bare volume integral and eV conversion

**Files:**
- Modify: `monolayer/src/MonolayerBareIntegrals.jl`
- Modify: `monolayer/test/test_monolayer_bare_integrals.jl`

The bare integral is `∫∫ f_src(r) f_tgt(r') / |r-r'| dr dr'`. We evaluate it by:
1. Computing the potential from `f_src` at each `f_tgt` quadrature node via `BI.TKM3D.ltkm3dc`.
2. Taking the `f_tgt`-weighted dot product: `raw = dot(u, target.weights .* target.density)`.
3. Converting to eV: `u_ev = raw * 4π * e²/(4πε₀) / (Na * Nb)` where `Na = sum(source.weights.*source.density)` is the source normalization (and similarly `Nb`). For the density-density channels these are ~1; for the Hund's channel they are ~0, so we **do not normalize Hund's by product norms** — instead we use the orbital normalizations `||phi_i||² = 1` implicit in the XSF files. Concretely:
   - `:onsite`, `:nn`: divide by `(Na * Nb)` computed from the source/target VolumeSource densities.
   - `:hund_sf`, `:hund_ph`: divide by `(Nphi1 * Nphi2)` where `Nphi_i = ∫ |phi_i|^2 dr` (available from the `densities` VolumeSources), **not** by the product-integral which would be ~0.

This matches the orbital-average convention used by CoQui (the matrix element between normalized Wannier orbitals).

- [ ] **Step 1: Write failing test — density-density onsite lands near reference**

Append to `monolayer/test/test_monolayer_bare_integrals.jl`:

```julia
@testset "bare_channel_integral — onsite ≈ 17.43 eV, nn ≈ 8.84 eV" begin
    densities = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = true)
    signed    = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = false)

    onsite = bare_channel_integral(densities, signed, :onsite; volume_tol = 1e-3)
    @test abs(onsite.u_ev - 17.434191) / 17.434191 < 0.05  # 5% tolerance at this grid/source_tol

    nn = bare_channel_integral(densities, signed, :nn; volume_tol = 1e-3)
    @test abs(nn.u_ev - 8.839615) / 8.839615 < 0.05
end

@testset "bare_channel_integral — Hund's ≈ 0.131 eV and spin-flip == pair-hopping" begin
    densities = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = true)
    signed    = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = false)

    j_sf = bare_channel_integral(densities, signed, :hund_sf; volume_tol = 1e-3)
    j_ph = bare_channel_integral(densities, signed, :hund_ph; volume_tol = 1e-3)

    @test abs(j_sf.u_ev - 0.130804) / 0.130804 < 0.15  # Hund's is small — allow 15% at coarse tol
    @test abs(j_ph.u_ev - 0.130804) / 0.130804 < 0.15
    @test abs(j_sf.u_ev - j_ph.u_ev) < 1e-6            # independent-calc consistency
end
```

- [ ] **Step 2: Run the test — expect FAIL**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: FAIL — `bare_channel_integral` undefined.

- [ ] **Step 3: Implement `bare_channel_integral`**

Append to `monolayer/src/MonolayerBareIntegrals.jl` before `end # module`:

```julia
"""
    bare_channel_integral(densities, signed, channel; volume_tol = 1e-3, kmax = nothing)

Compute a single bare Coulomb matrix element in raw units and eV.

`densities` and `signed` are the outputs of `centered_monolayer_sources` with
`square = true` / `false`, sharing the same centering shift and `source_tol`.
"""
function bare_channel_integral(densities, signed, channel::Symbol; volume_tol::Real = 1e-3, kmax = nothing)
    pair = channel_pair_sources(densities, signed, channel)
    source = pair.source
    target = pair.target

    target_positions = target.positions
    target_weights = target.weights .* target.density

    charges = source.weights .* source.density
    resolved_kmax = isnothing(kmax) ? BI._estimate_tkm3dc_kmax(source) : Float64(kmax)
    values = BoundaryIntegral.TKM3D.ltkm3dc(
        Float64(volume_tol),
        source.positions;
        charges = charges,
        targets = target_positions,
        pgt = 1,
        kmax = resolved_kmax,
    )
    values.ier == 0 || error("TKM3D.ltkm3dc failed with ier=$(values.ier)")
    u_at_targets = real.(values.pottarg)
    u_raw = dot(u_at_targets, target_weights)

    # Normalize by orbital norms, not product norms (Hund's-safe convention).
    Nphi1 = sum(densities.vs1.weights .* densities.vs1.density)
    Nphi2 = sum(densities.vs2.weights .* densities.vs2.density)

    if channel === :onsite
        Na, Nb = Nphi1, Nphi1
    elseif channel === :nn
        Na, Nb = Nphi1, Nphi2
    elseif channel === :hund_sf || channel === :hund_ph
        Na, Nb = Nphi1, Nphi2
    else
        throw(ArgumentError("unknown channel $channel"))
    end

    u_ev = to_eV(u_raw, Na, Nb)
    return (
        channel = channel,
        u_raw = u_raw,
        u_ev = u_ev,
        Na = Na,
        Nb = Nb,
        n_source_points = size(source.positions, 2),
        n_target_points = size(target.positions, 2),
        tkm_kmax = resolved_kmax,
        volume_tol = Float64(volume_tol),
    )
end
```

- [ ] **Step 4: Run the test — expect PASS**

```bash
julia --project=. monolayer/test/runtests.jl
```

Expected: PASS. If the onsite/NN values are outside 5%, the loader or the normalization is wrong — do NOT loosen the tolerance; investigate.

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerBareIntegrals.jl monolayer/test/test_monolayer_bare_integrals.jl
git commit -m "feat(monolayer): direct bare 1/|r-r'| integral with eV conversion"
```

---

## Task 6: Synthetic Gaussian-pair validation

**Files:**
- Modify: `monolayer/test/test_monolayer_bare_integrals.jl`

Prove the `bare_channel_integral` primitive against a known analytic answer before we trust the XSF numbers.

Analytic reference: two normalized 3D Gaussians `g_α(r) = (α/π)^(3/2) exp(-α r²)` at distance `R` have Coulomb integral
`V(R) = erf(sqrt(α/2) R) / R` (atomic units, i.e., units where `e²/(4πε₀) = 1`). In our eV units with the `to_eV` conversion, `V_eV = erf(sqrt(α/2) R) / R * 4π * E2_4PIEPS0`. Wait — our convention bakes the `4π * E2_4PIEPS0` factor into `to_eV`, but that's designed to consume the extra `4π` that comes from writing `1/|r-r'|` vs `1/(4π|r-r'|)`. Since our `ltkm3dc` call uses bare `1/|r-r'|`, the analytic reference we should match in raw units is simply `erf(sqrt(α/2) R) / R`, and then we multiply by `E2_4PIEPS0` (not `4π * E2_4PIEPS0`) when checking against an eV value. **Implementation therefore checks against the raw value directly** to avoid any convention mismatch.

- [ ] **Step 1: Write the failing synthetic test**

Append to `monolayer/test/test_monolayer_bare_integrals.jl`:

```julia
using SpecialFunctions
import BoundaryIntegral as BI

@testset "bare integral against two-Gaussian analytic answer" begin
    α = 2.0
    R = 2.5  # separation in Å
    # Build two Gaussian densities on a cubic grid with enough extent.
    Lx = Ly = Lz = 8.0
    nx = ny = nz = 60
    xs = range(-Lx/2, Lx/2; length = nx)
    ys = range(-Ly/2, Ly/2; length = ny)
    zs = range(-Lz/2, Lz/2; length = nz)
    dens_a = Array{Float64}(undef, nx, ny, nz)
    dens_b = Array{Float64}(undef, nx, ny, nz)
    norm_prefactor = (α/π)^1.5
    for (iz, z) in enumerate(zs), (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        r2_a = x^2 + y^2 + z^2
        r2_b = x^2 + (y - R)^2 + z^2
        dens_a[ix, iy, iz] = norm_prefactor * exp(-α * r2_a)
        dens_b[ix, iy, iz] = norm_prefactor * exp(-α * r2_b)
    end
    origin = (-Lx/2, -Ly/2, -Lz/2)
    lattice = ((Lx, 0.0, 0.0), (0.0, Ly, 0.0), (0.0, 0.0, Lz))

    datagrid_a = (; origin = origin, lattice = lattice, values = dens_a)
    datagrid_b = (; origin = origin, lattice = lattice, values = dens_b)
    vs_a = BI.VolumeSource(datagrid_a, tol = 1e-4)
    vs_b = BI.VolumeSource(datagrid_b, tol = 1e-4)

    charges_a = vs_a.weights .* vs_a.density
    resolved_kmax = BI._estimate_tkm3dc_kmax(vs_a)
    vals = BI.TKM3D.ltkm3dc(1e-4, vs_a.positions;
                            charges = charges_a, targets = vs_b.positions,
                            pgt = 1, kmax = resolved_kmax)
    u_at_b = real.(vals.pottarg)
    u_raw = sum(u_at_b .* (vs_b.weights .* vs_b.density))

    analytic = erf(sqrt(α/2) * R) / R
    @test abs(u_raw - analytic) / analytic < 5e-3
end
```

Note: if the signature for constructing a `VolumeSource` from an ad-hoc NamedTuple datagrid differs in `BoundaryIntegral` (e.g., lattice vs primvec naming), adjust the field names to match `BI.read_xsf` outputs — inspect `dump(first(BI.read_xsf(XSF_1))[2])` once to confirm the expected shape, then adapt. The analytic test still compares against `erf(sqrt(α/2) R)/R`.

- [ ] **Step 2: Run — expect FAIL or PASS depending on field-name adaptation**

```bash
julia --project=. monolayer/test/runtests.jl
```

If the test errors on `BI.VolumeSource(datagrid_a, ...)` due to field naming, read `BI.read_xsf` output structure and adapt the NamedTuple keys; then re-run.

Expected once adapted: PASS (5e-3 agreement).

- [ ] **Step 3: Commit**

```bash
git add monolayer/test/test_monolayer_bare_integrals.jl
git commit -m "test(monolayer): synthetic Gaussian-pair analytic validation of bare integral"
```

---

## Task 7: `compute_all_bare_channels` convenience function

**Files:**
- Modify: `monolayer/src/MonolayerBareIntegrals.jl`
- Modify: `monolayer/test/test_monolayer_bare_integrals.jl`

- [ ] **Step 1: Failing test**

Append:

```julia
@testset "compute_all_bare_channels returns four labelled rows" begin
    rows = compute_all_bare_channels(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, volume_tol = 1e-3)
    @test length(rows) == 4
    labels = [r.channel for r in rows]
    @test Set(labels) == Set([:onsite, :nn, :hund_sf, :hund_ph])
    for r in rows
        @test r.u_ev > 0
    end
end
```

- [ ] **Step 2: Run — expect FAIL**

```bash
julia --project=. monolayer/test/runtests.jl
```

- [ ] **Step 3: Implement**

Append to `monolayer/src/MonolayerBareIntegrals.jl` before `end # module`:

```julia
const BARE_CHANNELS = (:onsite, :nn, :hund_sf, :hund_ph)

function compute_all_bare_channels(;
    orbital_1::AbstractString,
    orbital_2::AbstractString,
    source_tol::Real = 1e-3,
    volume_tol::Real = 1e-3,
)
    # Loader is imported at the call site; we re-invoke via the name in the calling module.
    densities = Main.MonolayerOrbitalLoader.centered_monolayer_sources(;
        orbital_1 = orbital_1, orbital_2 = orbital_2, source_tol = source_tol, square = true,
    )
    signed = Main.MonolayerOrbitalLoader.centered_monolayer_sources(;
        orbital_1 = orbital_1, orbital_2 = orbital_2, source_tol = source_tol, square = false,
    )
    rows = NamedTuple[]
    for channel in BARE_CHANNELS
        r = bare_channel_integral(densities, signed, channel; volume_tol = volume_tol)
        push!(rows, merge(r, (; source_tol = Float64(source_tol), shift_x = densities.shared_shift[1], shift_y = densities.shared_shift[2], shift_z = densities.shared_shift[3])))
    end
    return rows
end
```

If `Main.MonolayerOrbitalLoader` can't be reached (module scoping), replace with an explicit `using ..MonolayerOrbitalLoader` at the top of `MonolayerBareIntegrals.jl` (add the import) and call `centered_monolayer_sources` directly. The test will surface whichever path is broken.

- [ ] **Step 4: Run — expect PASS**

```bash
julia --project=. monolayer/test/runtests.jl
```

- [ ] **Step 5: Commit**

```bash
git add monolayer/src/MonolayerBareIntegrals.jl monolayer/test/test_monolayer_bare_integrals.jl
git commit -m "feat(monolayer): compute_all_bare_channels aggregate entry point"
```

---

## Task 8: Run script + CSV output

**Files:**
- Create: `monolayer/scripts/bare_monolayer_graphene.jl`
- Artifact: `monolayer/data/bare_monolayer_graphene.csv`

- [ ] **Step 1: Write the script**

Create `monolayer/scripts/bare_monolayer_graphene.jl`:

```julia
using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
include(joinpath(@__DIR__, "..", "src", "MonolayerBareIntegrals.jl"))
using .MonolayerOrbitalLoader
using .MonolayerBareIntegrals

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")
const SOURCE_TOL = 1e-3
const VOLUME_TOL = 1e-3

const REFERENCE_EV = Dict(
    :onsite  => 17.434191,
    :nn      => 8.839615,
    :hund_sf => 0.130804,
    :hund_ph => 0.130804,
)

const OUT_CSV = joinpath(@__DIR__, "..", "data", "bare_monolayer_graphene.csv")

function main()
    rows = compute_all_bare_channels(; orbital_1 = XSF_1, orbital_2 = XSF_2,
                                       source_tol = SOURCE_TOL, volume_tol = VOLUME_TOL)
    out_rows = NamedTuple[]
    for r in rows
        ref = REFERENCE_EV[r.channel]
        push!(out_rows, (
            channel = String(r.channel),
            u_raw = r.u_raw,
            u_ev = r.u_ev,
            u_ref_ev = ref,
            rel_err_pct = 100 * (r.u_ev - ref) / ref,
            Na = r.Na,
            Nb = r.Nb,
            n_source_points = r.n_source_points,
            n_target_points = r.n_target_points,
            tkm_kmax = r.tkm_kmax,
            source_tol = r.source_tol,
            volume_tol = r.volume_tol,
            shift_x = r.shift_x,
            shift_y = r.shift_y,
            shift_z = r.shift_z,
        ))
    end
    mkpath(dirname(OUT_CSV))
    table = DataFrame(out_rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    show(table; allrows = true, allcols = true)
    println()
end

main()
```

- [ ] **Step 2: Run the script**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
julia --project=. monolayer/scripts/bare_monolayer_graphene.jl
```

Expected: prints a 4-row table with `rel_err_pct` absolute values all below 5% (stretch: below 2%). `u_ev` close to `17.43`, `8.84`, `0.131`, `0.131`. CSV at `monolayer/data/bare_monolayer_graphene.csv`.

- [ ] **Step 3: Inspect the CSV**

```bash
cat monolayer/data/bare_monolayer_graphene.csv
```

- [ ] **Step 4: Commit**

```bash
git add monolayer/scripts/bare_monolayer_graphene.jl monolayer/data/bare_monolayer_graphene.csv
git commit -m "feat(monolayer): run bare-channel script and land first result CSV"
```

---

## Task 9: Mirror-pad convergence diagnostic

**Files:**
- Modify: `monolayer/src/MonolayerOrbitalLoader.jl`
- Modify: `monolayer/scripts/bare_monolayer_graphene.jl`
- Modify: `monolayer/data/bare_monolayer_graphene.csv`

"Mirror-pad level k" means: take the signed 3-D datagrid, reflect it `k` times about each of the x and y boundaries, stack the reflections into a `(2k+1)×(2k+1)×1` tiled grid covering `k` extra image cells on each side in plane. Level `0` is the raw XSF grid. We keep z fixed (15 Å vacuum is already large).

The purpose is a sanity check: if mirror-pad level 1 shifts any channel by more than ~2% relative to level 0, the xy tail is not converged and we should switch to the FFT/madelung fallback (out of scope here — flag it in the report).

- [ ] **Step 1: Add `mirror_pad_xy` to the loader**

Append to `monolayer/src/MonolayerOrbitalLoader.jl` before `end # module`:

```julia
export mirror_pad_xy

"""
    mirror_pad_xy(datagrid, level::Integer)

Return a new datagrid tiled `(2 level + 1)` times in x and y using mirror
reflections. Level 0 is a no-op. z is untouched.
"""
function mirror_pad_xy(datagrid, level::Integer)
    level >= 0 || throw(ArgumentError("level must be ≥ 0"))
    level == 0 && return datagrid

    vals = datagrid.values
    nx, ny, nz = size(vals)
    tile = 2 * level + 1

    # Build a (tile*nx, tile*ny, nz) array by reflecting vals along x and y for each tile index.
    padded = Array{eltype(vals)}(undef, tile * nx, tile * ny, nz)
    for jy in 0:tile-1
        flip_y = isodd(jy - level) ? true : false  # center tile (jy == level) is not flipped
        src_y = flip_y ? reverse(vals, dims = 2) : vals
        for jx in 0:tile-1
            flip_x = isodd(jx - level) ? true : false
            src = flip_x ? reverse(src_y, dims = 1) : src_y
            x_range = (jx * nx + 1):((jx + 1) * nx)
            y_range = (jy * ny + 1):((jy + 1) * ny)
            padded[x_range, y_range, :] .= src
        end
    end

    lat = datagrid.lattice
    a1 = lat[1]; a2 = lat[2]; a3 = lat[3]
    new_lat = (
        (tile * a1[1], tile * a1[2], tile * a1[3]),
        (tile * a2[1], tile * a2[2], tile * a2[3]),
        a3,
    )
    # Shift origin so the center tile aligns with the original datagrid.
    o = datagrid.origin
    new_origin = (
        o[1] - level * a1[1] - level * a2[1],
        o[2] - level * a1[2] - level * a2[2],
        o[3],
    )

    return merge(datagrid, (; values = padded, lattice = new_lat, origin = new_origin))
end
```

- [ ] **Step 2: Expose a `mirror_pad_level` parameter in `bare_channel_integral` (pass-through via new loader call)**

Instead of deeply threading a parameter through `centered_monolayer_sources`, the script-level logic will:
1. Load orbital 1 and 2 datagrids (both signed and squared).
2. Apply `mirror_pad_xy` at the desired level to each.
3. Compute the centering shift from the padded squared orbital 1.
4. Build `VolumeSource` objects from padded+shifted grids.
5. Hand these to `bare_channel_integral` through a tiny helper that mirrors the shape of `centered_monolayer_sources`' return value but with the padded inputs.

Add to `monolayer/src/MonolayerOrbitalLoader.jl` before `end # module`:

```julia
export centered_monolayer_sources_padded

"""
    centered_monolayer_sources_padded(; orbital_1, orbital_2, source_tol, square, mirror_pad_level)

Same as `centered_monolayer_sources`, but each orbital is first mirror-padded
`mirror_pad_level` times in x and y before centering and VolumeSource construction.
"""
function centered_monolayer_sources_padded(;
    orbital_1::AbstractString,
    orbital_2::AbstractString,
    source_tol::Real = 1e-3,
    square::Bool = true,
    mirror_pad_level::Integer = 0,
)
    base_1 = square ? load_squared_xsf(orbital_1) : load_signed_xsf(orbital_1)
    base_2 = square ? load_squared_xsf(orbital_2) : load_signed_xsf(orbital_2)
    padded_1 = mirror_pad_xy(base_1, mirror_pad_level)
    padded_2 = mirror_pad_xy(base_2, mirror_pad_level)

    centering_1 = mirror_pad_xy(load_squared_xsf(orbital_1), mirror_pad_level)
    shared_shift = _centering_shift(centering_1; tol = source_tol)

    padded_1_shifted = shift_datagrid(padded_1, shared_shift)
    padded_2_shifted = shift_datagrid(padded_2, shared_shift)

    vs1 = BI.VolumeSource(padded_1_shifted, tol = source_tol)
    vs2 = BI.VolumeSource(padded_2_shifted, tol = source_tol)

    return (
        datagrid_1 = padded_1_shifted,
        datagrid_2 = padded_2_shifted,
        vs1 = vs1,
        vs2 = vs2,
        shared_shift = shared_shift,
    )
end
```

- [ ] **Step 3: Update `compute_all_bare_channels` to accept `mirror_pad_level`**

Replace the body of `compute_all_bare_channels` in `monolayer/src/MonolayerBareIntegrals.jl` with:

```julia
function compute_all_bare_channels(;
    orbital_1::AbstractString,
    orbital_2::AbstractString,
    source_tol::Real = 1e-3,
    volume_tol::Real = 1e-3,
    mirror_pad_level::Integer = 0,
)
    densities = Main.MonolayerOrbitalLoader.centered_monolayer_sources_padded(;
        orbital_1 = orbital_1, orbital_2 = orbital_2, source_tol = source_tol, square = true,
        mirror_pad_level = mirror_pad_level,
    )
    signed = Main.MonolayerOrbitalLoader.centered_monolayer_sources_padded(;
        orbital_1 = orbital_1, orbital_2 = orbital_2, source_tol = source_tol, square = false,
        mirror_pad_level = mirror_pad_level,
    )
    rows = NamedTuple[]
    for channel in BARE_CHANNELS
        r = bare_channel_integral(densities, signed, channel; volume_tol = volume_tol)
        push!(rows, merge(r, (;
            source_tol = Float64(source_tol),
            mirror_pad_level = Int(mirror_pad_level),
            shift_x = densities.shared_shift[1],
            shift_y = densities.shared_shift[2],
            shift_z = densities.shared_shift[3],
        )))
    end
    return rows
end
```

- [ ] **Step 4: Update the run script to sweep `mirror_pad_level ∈ {0, 1}`**

Edit `monolayer/scripts/bare_monolayer_graphene.jl`. Replace the body of `main` with:

```julia
function main()
    all_rows = NamedTuple[]
    for level in (0, 1)
        println("--- mirror_pad_level = $level ---")
        rows = compute_all_bare_channels(; orbital_1 = XSF_1, orbital_2 = XSF_2,
                                           source_tol = SOURCE_TOL, volume_tol = VOLUME_TOL,
                                           mirror_pad_level = level)
        for r in rows
            ref = REFERENCE_EV[r.channel]
            push!(all_rows, (
                channel = String(r.channel),
                mirror_pad_level = r.mirror_pad_level,
                u_raw = r.u_raw,
                u_ev = r.u_ev,
                u_ref_ev = ref,
                rel_err_pct = 100 * (r.u_ev - ref) / ref,
                Na = r.Na,
                Nb = r.Nb,
                n_source_points = r.n_source_points,
                n_target_points = r.n_target_points,
                tkm_kmax = r.tkm_kmax,
                source_tol = r.source_tol,
                volume_tol = r.volume_tol,
                shift_x = r.shift_x,
                shift_y = r.shift_y,
                shift_z = r.shift_z,
            ))
        end
    end
    mkpath(dirname(OUT_CSV))
    table = DataFrame(all_rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    show(table; allrows = true, allcols = true)
    println()
end
```

- [ ] **Step 5: Re-run the script**

```bash
julia --project=. monolayer/scripts/bare_monolayer_graphene.jl
```

Expected: 8 rows (4 channels × 2 pad levels). Compute for each channel `|u_ev(level=1) - u_ev(level=0)| / u_ev(level=0)` — should be sub-percent for `:onsite`, `:nn`; possibly a few percent for the small Hund's integrals. If any density-channel drift exceeds 2%, note it in the report and flag the FFT+madelung fallback as needed.

- [ ] **Step 6: Commit**

```bash
git add monolayer/src/MonolayerOrbitalLoader.jl monolayer/src/MonolayerBareIntegrals.jl monolayer/scripts/bare_monolayer_graphene.jl monolayer/data/bare_monolayer_graphene.csv
git commit -m "feat(monolayer): mirror-pad convergence diagnostic at level {0, 1}"
```

---

## Task 10: Report

**Files:**
- Create: `monolayer/results/2026-04-23-bare-monolayer-report.md`

- [ ] **Step 1: Write the report**

Template to fill in with the actual numbers from the CSV:

```markdown
# Bare Monolayer Graphene Interactions — Report

Date: 2026-04-23

## Reference (CoQui cRPA; Malte Rösner)

Source: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15/_coqui_thc_crpa.out`

Bare interactions (orbital-average):

- intra-orbital (onsite)  = 17.434191 eV
- inter-orbital (NN)      = 8.839615 eV
- Hund's (spin-flip)      = 0.130804 eV
- Hund's (pair-hopping)   = 0.130804 eV

## This Run

Orbitals: `graphene_00001.xsf`, `graphene_00002.xsf` (150×150×192 grid over ~12.24 × 10.60 × 14.92 Å).

Parameters: `source_tol = 1e-3`, `volume_tol = 1e-3`, direct `1/|r-r'|` via `BI.TKM3D.ltkm3dc`, no dielectric.

| Channel  | `u_ev` | `u_ref_ev` | `rel_err_pct` | pad level |
|----------|-------:|-----------:|--------------:|----------:|
| onsite   | <fill> | 17.434191  | <fill>        | 0         |
| nn       | <fill> | 8.839615   | <fill>        | 0         |
| hund_sf  | <fill> | 0.130804   | <fill>        | 0         |
| hund_ph  | <fill> | 0.130804   | <fill>        | 0         |

## Mirror-Pad Convergence

| Channel | pad 0 `u_ev` | pad 1 `u_ev` | drift % |
|---------|-------------:|-------------:|--------:|
| onsite  | <fill>       | <fill>       | <fill>  |
| nn      | <fill>       | <fill>       | <fill>  |
| hund_sf | <fill>       | <fill>       | <fill>  |
| hund_ph | <fill>       | <fill>       | <fill>  |

## Verdict

- density channels (onsite, nn) within <fill>% of the cRPA bare reference
- Hund's channels within <fill>% of the cRPA bare reference
- Spin-flip and pair-hopping channels agree numerically (difference < 1e-6 eV), confirming the real-orbital identity
- Mirror-pad drift is <fill>%, indicating <fill: xy tail is converged / xy tail requires FFT+madelung follow-up>

## Next

- await Malte's effective-thickness / ε_eff(q) extraction script for the screened channel validation
- await his bilayer monolayer example before any bilayer revision
```

- [ ] **Step 2: Fill in the `<fill>` values from `monolayer/data/bare_monolayer_graphene.csv`**

Open the CSV, copy the numbers into the appropriate cells.

- [ ] **Step 3: Commit**

```bash
git add monolayer/results/2026-04-23-bare-monolayer-report.md
git commit -m "docs(monolayer): add bare-channel results report"
```

---

## Task 11: Final end-to-end verification

- [ ] **Step 1: Run the full test suite**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
julia --project=. test/runtests.jl
```

Expected: all bilayer_slab + monolayer tests pass.

- [ ] **Step 2: Re-run the monolayer script clean**

```bash
julia --project=. monolayer/scripts/bare_monolayer_graphene.jl
```

Expected: 8-row CSV output, console summary matches committed `bare_monolayer_graphene.csv`.

- [ ] **Step 3: `git status` clean; log summary**

```bash
git status
git log --oneline -15
```

Expected: working tree clean; the new commits listed in order.

---

## Self-Review

Spec coverage:

- Repo reorganization under `bilayer_slab/` + dispatcher test: Task 1.
- New `monolayer/` tree: Tasks 2–10.
- Four bare channels (`onsite`, `nn`, `hund_sf`, `hund_ph`): Tasks 4, 5, 7, 8.
- Signed-orbital pipeline for Hund's: Task 2 (loader), Task 4 (product datagrid), Task 5 (integrator).
- Mirror-pad convergence check: Task 9.
- Synthetic Gaussian validation: Task 6.
- CSV output + report: Tasks 8 and 10.
- Out-of-scope items (screened monolayer, bilayer example, FFT+madelung fallback, extra shells) are left untouched, as the spec requires.

Known loose ends made explicit inside the plan:

- Task 6 explicitly flags that the synthetic test may need field-name adaptation for `BI.VolumeSource`. The plan instructs the implementer to `dump` one real XSF-loaded datagrid and align names.
- Task 7 warns that `Main.MonolayerOrbitalLoader` scoping may break; a fallback (`using ..MonolayerOrbitalLoader`) is provided.
- Task 9 frames a > 2% drift as a trigger for a separate FFT+madelung plan rather than scope creep here.

module MonolayerScreenedSolve

using BoundaryIntegral
import BoundaryIntegral as BI

include(joinpath(@__DIR__, "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

export parse_coqui_loc, monolayer_screened_sources, density_target_specs, hund_target_specs

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

end # module

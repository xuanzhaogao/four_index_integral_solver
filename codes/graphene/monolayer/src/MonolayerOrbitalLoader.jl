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

export mirror_pad_xy

"""
    mirror_pad_xy(datagrid, level::Integer)

Return a new datagrid tiled `(2 level + 1)` times in x and y using mirror
reflections. Level 0 is a no-op. z is untouched.

The BoundaryIntegral XSF datagrid follows the Wannier90 half-open convention,
with fields `nx, ny, nz, origin, A, B, C, values`. The stored `A`, `B`, `C`
span `(n-1)/n` of the "true" periodic cell; `true_cell_vectors` scales them
back. After mirror padding we keep that convention intact by computing the
true cell, extending it, and re-applying the `(n_new - 1)/n_new` factor.
"""
function mirror_pad_xy(datagrid, level::Integer)
    level >= 0 || throw(ArgumentError("level must be ≥ 0"))
    level == 0 && return datagrid

    vals = datagrid.values
    nx, ny, nz = size(vals)
    tile = 2 * level + 1

    padded = Array{eltype(vals)}(undef, tile * nx, tile * ny, nz)
    for jy in 0:tile-1
        flip_y = isodd(jy - level)          # center tile (jy == level) is not flipped
        src_y = flip_y ? reverse(vals, dims = 2) : vals
        for jx in 0:tile-1
            flip_x = isodd(jx - level)
            src = flip_x ? reverse(src_y, dims = 1) : src_y
            x_range = (jx * nx + 1):((jx + 1) * nx)
            y_range = (jy * ny + 1):((jy + 1) * ny)
            padded[x_range, y_range, :] .= src
        end
    end

    # Rescale A, B per the Wannier90 half-open convention so true_cell_vectors
    # produces tile*At, tile*Bt after padding.
    A = datagrid.A
    B = datagrid.B
    At = ntuple(i -> A[i] * nx / (nx - 1), 3)
    Bt = ntuple(i -> B[i] * ny / (ny - 1), 3)
    new_nx = tile * nx
    new_ny = tile * ny
    new_At = ntuple(i -> tile * At[i], 3)
    new_Bt = ntuple(i -> tile * Bt[i], 3)
    new_A = ntuple(i -> new_At[i] * (new_nx - 1) / new_nx, 3)
    new_B = ntuple(i -> new_Bt[i] * (new_ny - 1) / new_ny, 3)

    # Shift the origin so the center tile occupies its original position.
    o = datagrid.origin
    new_origin = ntuple(i -> o[i] - level * At[i] - level * Bt[i], 3)

    return merge(datagrid, (;
        values = padded,
        nx = new_nx,
        ny = new_ny,
        A = new_A,
        B = new_B,
        origin = new_origin,
    ))
end

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

end # module

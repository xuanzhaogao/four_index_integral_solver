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

end # module

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

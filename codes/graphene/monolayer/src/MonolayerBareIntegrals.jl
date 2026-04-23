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

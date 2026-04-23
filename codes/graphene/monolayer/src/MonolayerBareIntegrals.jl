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

end # module

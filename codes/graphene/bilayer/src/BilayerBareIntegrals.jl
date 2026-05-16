module BilayerBareIntegrals

# Bilayer bare-Coulomb integrals computed by reusing the monolayer pipeline
# (MonolayerOrbitalLoader + MonolayerBareIntegrals) on every distinct orbital
# pair of a 4-orbital bilayer Wannier basis.
#
# Channels per pair (i, j):
#   - i == j → only :onsite (the other channels would just duplicate it).
#   - i <  j → :nn (|phi_i|^2 vs |phi_j|^2), :hund_sf, :hund_ph.

using Statistics

export compute_all_bilayer_bare_channels, BILAYER_REF_AVG_EV

# CoQui orbital-averaged bare reference for k_161601_nb_288_c_15:
# /mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/bilayer/k_161601_nb_288_c_15/_coqui_thc_crpa.out
const BILAYER_REF_AVG_EV = (
    intra   = 17.669390,
    inter   = 5.785749,
    hund_sf = 0.041867,
    hund_ph = 0.041867,
)

"""
    compute_all_bilayer_bare_channels(; orbitals, source_tol, volume_tol, mirror_pad_level)

Run the monolayer per-pair bare-channel solver on every distinct orbital pair of
a 4-orbital bilayer XSF set. `orbitals` is a `NTuple{4, String}` of paths.

Returns a `Vector{NamedTuple}` with one row per (i, j, channel), annotated with
the orbital indices, the channel name, raw & eV values, normalizations,
TKM kmax, and the shared centering shift applied per pair.

Resolution of the monolayer driver is dynamic via `Main.MonolayerBareIntegrals`,
matching the lookup pattern in `MonolayerBareIntegrals.compute_all_bare_channels`
itself (the loader / integrator modules are expected to be `include`d into
`Main` by the calling script).
"""
function compute_all_bilayer_bare_channels(;
    orbitals::NTuple{4, <:AbstractString},
    source_tol::Real = 1e-3,
    volume_tol::Real = 1e-3,
    mirror_pad_level::Integer = 0,
)
    monolayer_solver = Main.MonolayerBareIntegrals.compute_all_bare_channels
    rows = NamedTuple[]
    for i in 1:4, j in i:4
        keep = (i == j) ? (:onsite,) : (:nn, :hund_sf, :hund_ph)
        sub = monolayer_solver(;
            orbital_1 = orbitals[i],
            orbital_2 = orbitals[j],
            source_tol = source_tol,
            volume_tol = volume_tol,
            mirror_pad_level = mirror_pad_level,
        )
        for r in sub
            r.channel in keep || continue
            push!(rows, (
                i = i,
                j = j,
                channel = String(r.channel),
                u_raw = r.u_raw,
                u_ev = r.u_ev,
                Na = r.Na,
                Nb = r.Nb,
                n_source_points = r.n_source_points,
                n_target_points = r.n_target_points,
                tkm_kmax = r.tkm_kmax,
                source_tol = r.source_tol,
                volume_tol = r.volume_tol,
                mirror_pad_level = r.mirror_pad_level,
                shift_x = r.shift_x,
                shift_y = r.shift_y,
                shift_z = r.shift_z,
            ))
        end
    end
    return rows
end

end # module

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

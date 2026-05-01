using CSV
using DataFrames
using Printf

include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
include(joinpath(@__DIR__, "..", "src", "MonolayerBareIntegrals.jl"))
using .MonolayerOrbitalLoader
using .MonolayerBareIntegrals

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")

const TOL_GRID = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
const MIRROR_PAD_LEVEL = 0

const REFERENCE_EV = Dict(
    :onsite  => 17.434191,
    :nn      => 8.839615,
    :hund_sf => 0.130804,
    :hund_ph => 0.130804,
)

const PRL_REF_EV = Dict(
    :onsite  => 17.0,
    :nn      => 8.5,
)

const OUT_CSV = joinpath(@__DIR__, "..", "data", "bare_monolayer_tol_sweep.csv")

function main()
    println("Threads: ", Threads.nthreads())
    println("Tolerance grid: ", TOL_GRID)
    all_rows = NamedTuple[]
    for tol in TOL_GRID
        @printf("\n=== source_tol = volume_tol = %.1e ===\n", tol)
        flush(stdout)
        t0 = time()
        rows = compute_all_bare_channels(; orbital_1 = XSF_1, orbital_2 = XSF_2,
                                          source_tol = tol, volume_tol = tol,
                                          mirror_pad_level = MIRROR_PAD_LEVEL)
        dt = time() - t0
        @printf("  elapsed: %.1f s\n", dt)
        for r in rows
            ref = REFERENCE_EV[r.channel]
            push!(all_rows, (
                channel = String(r.channel),
                tol = Float64(tol),
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
                elapsed_s = dt,
            ))
            @printf("  %-8s u_ev = %10.6f  ref = %10.6f  err = %+6.2f%%  n_src = %d\n",
                    String(r.channel), r.u_ev, ref, 100*(r.u_ev - ref)/ref, r.n_source_points)
            flush(stdout)
        end
        # Incremental write so partial results survive a crash / kill.
        mkpath(dirname(OUT_CSV))
        CSV.write(OUT_CSV, DataFrame(all_rows))
    end
    println("\nWrote $(length(all_rows)) rows to $(OUT_CSV)")
    show(DataFrame(all_rows); allrows = true, allcols = true)
    println()
end

main()

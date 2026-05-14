using BoundaryIntegral
import BoundaryIntegral as BI

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve

include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1   = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2   = joinpath(REF_DIR, "graphene_00002.xsf")
const COQUI   = joinpath(REF_DIR, "_coqui_crpa_loc.out")

# Calibrated slab parameters from 2026-05-14 eps_in sweep
const LZ       = 3.35
const Z_CENTER = LZ / 2
const EPS_IN   = 3.5
const EPS_OUT  = 1.0
const L        = 90.0

const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12

const OUT_CSV = joinpath(@__DIR__, "..", "data", "screened_monolayer_hund_calibrated.csv")

function main()
    println("Parsing CoQui reference ...")
    coqui = parse_coqui_loc(COQUI)
    println("  hund_sf U_ijkl = ", coqui[:hund_sf].U_ijkl)

    println("Building sources at z_center = ", Z_CENTER, " ...")
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER,
    )
    println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

    hund_specs = hund_target_specs(src)
    println("--- Hund's solve  eps_in = ", EPS_IN, "  LZ = ", LZ, " ---")
    result = solve_screened_mode(
        src.vs_product, hund_specs,
        L, L, LZ, EPS_IN, EPS_OUT, BI.SharpScreening();
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
        rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH,
        volume_tol = RHS_TOL,
    )
    println("  residual = ", result.residual,
            "  n_interface_points = ", result.n_interface_points)

    rows = NamedTuple[]
    for pr in result.pair_results
        row0 = screened_run_record(pr, result;
            mode_label = "Sharp", bandwidth = missing,
            eps_in = EPS_IN, eps_out = EPS_OUT, Lz = LZ, L = L, z_center = Z_CENTER,
            source_tol = SOURCE_TOL, rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
            gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
            n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
            max_order = MAX_ORDER, max_depth = MAX_DEPTH,
        )
        U_ref = coqui[:hund_sf].U_ijkl
        row = merge(row0, (
            coqui_v_ijkl = coqui[:hund_sf].v_ijkl,
            coqui_U_ijkl = U_ref,
            diff_ev = U_ref - row0.u_total_ev,
            rel_err_pct = 100 * (row0.u_total_ev - U_ref) / U_ref,
        ))
        push!(rows, row)
        println("    ", pr.pair, "  ours=", row0.u_total_ev,
                "  CoQui=", U_ref,
                "  rel_err_pct=", row.rel_err_pct)
    end

    mkpath(dirname(OUT_CSV))
    table = DataFrame(rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    summary_cols = [:channel, :eps_in, :Lz, :u_int_ev, :u_scatter_ev, :u_total_ev, :coqui_U_ijkl, :rel_err_pct, :sigma_residual]
    show(table[:, summary_cols]; allrows = true, allcols = true)
    println()
end

main()

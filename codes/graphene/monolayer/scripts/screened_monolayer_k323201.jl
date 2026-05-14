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

# Slab geometry
const LZ = 3.35
const Z_CENTER = LZ / 2
const EPS_IN = 2.4
const EPS_OUT = 1.0
const L = 90.0

# Solver tolerances (reused from the bilayer protocol)
const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12

# Density-channel bandwidth sweep (used by Task 7 main; Task 6 smoke uses only Sharp)
const BANDWIDTHS = [0.05, 0.1, 0.2, 0.4]

const OUT_CSV = joinpath(@__DIR__, "..", "data", "screened_monolayer_k323201.csv")

function _solve_kwargs()
    return (
        n_quad = N_QUAD,
        edge_refine_level = EDGE_REFINE_LEVEL,
        rhs_tol = RHS_TOL,
        lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL,
        gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER,
        max_depth = MAX_DEPTH,
        volume_tol = RHS_TOL,
    )
end

function _record_kwargs(mode_label, bandwidth)
    return (
        mode_label = mode_label, bandwidth = bandwidth,
        eps_in = EPS_IN, eps_out = EPS_OUT, Lz = LZ, L = L, z_center = Z_CENTER,
        source_tol = SOURCE_TOL, rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH,
    )
end

function _density_modes()
    cases = NamedTuple{(:label, :bandwidth, :mode), Tuple{String, Union{Missing, Float64}, BI.AbstractScreeningMode}}[]
    push!(cases, (label = "Sharp", bandwidth = missing, mode = BI.SharpScreening()))
    for bw in BANDWIDTHS
        push!(cases, (label = "SoftMixInversePermittivity", bandwidth = bw, mode = BI.SoftMixInversePermittivity(bw)))
    end
    return cases
end

function _hund_modes()
    return [
        (label = "Sharp", bandwidth = missing, mode = BI.SharpScreening()),
    ]
end

function main()
    println("Loading sources ...")
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER,
    )
    println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

    coqui = parse_coqui_loc(COQUI)

    rows = NamedTuple[]

    density_specs = density_target_specs(src)
    for case in _density_modes()
        println("Density solve  mode=", case.label, "  bandwidth=", case.bandwidth)
        result = solve_screened_mode(
            src.vs1, density_specs,
            L, L, LZ, EPS_IN, EPS_OUT, case.mode;
            _solve_kwargs()...,
        )
        println("  residual = ", result.residual,
                "  n_interface_points = ", result.n_interface_points)
        for pr in result.pair_results
            ch_sym = pr.pair == "onsite" ? :onsite : :nn
            row0 = screened_run_record(pr, result; _record_kwargs(case.label, case.bandwidth)...)
            row = merge(row0, (
                coqui_v_ijkl = coqui[ch_sym].v_ijkl,
                coqui_U_ijkl = coqui[ch_sym].U_ijkl,
                diff_ev = coqui[ch_sym].U_ijkl - row0.u_total_ev,
            ))
            push!(rows, row)
            println("    ", pr.pair, "  ours=", row0.u_total_ev,
                    "  CoQui=", coqui[ch_sym].U_ijkl,
                    "  diff=", row.diff_ev)
        end
    end

    hund_specs = hund_target_specs(src)
    for case in _hund_modes()
        println("Hund's solve  mode=", case.label, "  bandwidth=", case.bandwidth)
        result = solve_screened_mode(
            src.vs_product, hund_specs,
            L, L, LZ, EPS_IN, EPS_OUT, case.mode;
            _solve_kwargs()...,
        )
        println("  residual = ", result.residual,
                "  n_interface_points = ", result.n_interface_points)
        for pr in result.pair_results
            row0 = screened_run_record(pr, result; _record_kwargs(case.label, case.bandwidth)...)
            row = merge(row0, (
                coqui_v_ijkl = coqui[:hund_sf].v_ijkl,
                coqui_U_ijkl = coqui[:hund_sf].U_ijkl,
                diff_ev = coqui[:hund_sf].U_ijkl - row0.u_total_ev,
            ))
            push!(rows, row)
            println("    ", pr.pair, "  ours=", row0.u_total_ev,
                    "  CoQui=", coqui[:hund_sf].U_ijkl,
                    "  diff=", row.diff_ev)
        end
    end

    mkpath(dirname(OUT_CSV))
    table = DataFrame(rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    summary_cols = [:channel, :mode, :bandwidth, :u_total_ev, :coqui_U_ijkl, :diff_ev, :sigma_residual]
    show(table[:, summary_cols]; allrows = true, allcols = true)
    println()
end

main()

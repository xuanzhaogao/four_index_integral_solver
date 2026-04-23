using BoundaryIntegral
import BoundaryIntegral as BI

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve

const L = 90.0
const LZ = 6.7
const Z_CENTER = LZ / 4
const EPS_IN = 2.4
const EPS_OUT = 1.0
const BANDWIDTHS = [0.05, 0.1, 0.2, 0.4]
const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12
const OUTPUT_FILE = joinpath(@__DIR__, "..", "data", "screened_hubbard_graphene.csv")

function mode_cases(args = ARGS)
    requested = isempty(args) ? nothing : args
    cases = NamedTuple{(:label, :bandwidth, :mode), Tuple{String, Union{Missing, Float64}, BI.AbstractScreeningMode}}[]
    if isnothing(requested) || any(lowercase(arg) == "sharp" for arg in requested)
        push!(cases, (label = "Sharp", bandwidth = missing, mode = BI.SharpScreening()))
    end

    bandwidths = if isnothing(requested)
        BANDWIDTHS
    else
        parsed = Float64[]
        for arg in requested
            lowercase(arg) == "sharp" && continue
            push!(parsed, parse(Float64, arg))
        end
        parsed
    end
    for bandwidth in bandwidths
        push!(cases, (label = "SoftMixInversePermittivity", bandwidth = bandwidth, mode = BI.SoftMixInversePermittivity(bandwidth)))
    end
    return cases
end

function main(args = ARGS)
    sources = centered_graphene_sources(tol = SOURCE_TOL, z_center = Z_CENTER)
    specs = pair_specs(sources.vs1, sources.vs2)

    rows = NamedTuple[]
    for case in mode_cases(args)
        println("Solving mode $(case.label) bandwidth=$(case.bandwidth)")
        result = solve_screened_mode(
            sources.vs1,
            specs,
            L,
            L,
            LZ,
            EPS_IN,
            EPS_OUT,
            case.mode;
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

        for pair_result in result.pair_results
            push!(rows, (
                mode = case.label,
                bandwidth = case.bandwidth,
                pair = pair_result.pair,
                u_int_raw = pair_result.u_int_raw,
                u_scatter_raw = pair_result.u_scatter_raw,
                u_total_raw = pair_result.u_total_raw,
                u_int_ev = pair_result.u_int_ev,
                u_scatter_ev = pair_result.u_scatter_ev,
                u_total_ev = pair_result.u_total_ev,
                source_norm = pair_result.Na,
                target_norm = pair_result.Nb,
                n_source_points = length(result.screened_vs.density),
                n_target_points = length(pair_result.target_vs.density),
                n_interface_points = result.n_interface_points,
                sigma_residual = result.residual,
                tkm_kmax = result.tkm_kmax,
                source_tol = SOURCE_TOL,
                n_quad = N_QUAD,
                edge_refine_level = EDGE_REFINE_LEVEL,
                rhs_tol = RHS_TOL,
                lhs_tol = LHS_TOL,
                gmres_atol = GMRES_ATOL,
                gmres_rtol = GMRES_RTOL,
                slab_thickness = LZ,
                orbital_z_center = Z_CENTER,
                shift_x = sources.shared_shift[1],
                shift_y = sources.shared_shift[2],
                shift_z = sources.shared_shift[3],
            ))
        end
    end

    table = DataFrame(rows)
    CSV.write(OUTPUT_FILE, table)
    println("Wrote $(nrow(table)) rows to $(OUTPUT_FILE)")
    show(table; allrows = true, allcols = true)
    println()
end

main()

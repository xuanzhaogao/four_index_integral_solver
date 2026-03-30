using BoundaryIntegral
const BI = BoundaryIntegral

using LinearAlgebra
using Krylov
using CSV
using DataFrames
using Printf
using Logging

const CONFIGS = [
    (2, 2),
    (2, 3),
    (2, 4),
    (4, 2),
    (4, 3),
    (4, 4),
    (6, 3),
    (6, 4),
]

const SUBSTRATE_CENTER = (0.0, 0.0, -15.0)
const SUBSTRATE_SIZE = (30.0, 30.0, 30.0)
const FILM_CENTER = (0.0, 0.0, 1.4)
const FILM_SIZE = (20.0, 20.0, 2.8)

const EPSES = [10.0, 2.4]
const EPS_OUT = 1.0
const EPS_SRC = 2.4

const A_CC = 1.42
const A1 = 2.465
const TARGETS = [
    (0.0, 0.0, 1.4),
    (A_CC, 0.0, 1.4),
    (A1, 0.0, 1.4),
    (A1 + A_CC, 0.0, 1.4),
]
const TARGET_LABELS = ("U_00", "U_01", "U_02", "U_03")

const PS = BI.PointSource((0.0, 0.0, 1.4), 1.0)
const FMM_TOL = 1.0e-6
const MAX_ORDER = 128
const GMRES_RTOL = 1.0e-8

const OUTPUT_DIR = "data"
const OUTPUT_PATH = joinpath(OUTPUT_DIR, "convergence.csv")

function empty_results()
    DataFrame(
        n_quad = Int[],
        edge_refine = Int[],
        n_pts = Int[],
        n_iter = Int[],
        residual = Float64[],
        U_00 = Float64[],
        U_01 = Float64[],
        U_02 = Float64[],
        U_03 = Float64[],
    )
end

function load_results(path::AbstractString)
    isfile(path) ? DataFrame(CSV.File(path)) : empty_results()
end

function completed_configs(results::DataFrame)
    Set((Int(row.n_quad), Int(row.edge_refine)) for row in eachrow(results))
end

function sorted_results(results::DataFrame)
    nrow(results) == 0 ? results : sort(results, [:n_quad, :edge_refine])
end

function persist_results(results::DataFrame, path::AbstractString)
    mkpath(dirname(path))
    CSV.write(path, sorted_results(results))
    nothing
end

function targets_matrix_3xn(targets)
    mat = Matrix{Float64}(undef, 3, length(targets))
    for (j, trg) in enumerate(targets)
        mat[:, j] .= trg
    end
    mat
end

const BOXES = [
    (center = SUBSTRATE_CENTER, Lx = SUBSTRATE_SIZE[1], Ly = SUBSTRATE_SIZE[2], Lz = SUBSTRATE_SIZE[3]),
    (center = FILM_CENTER, Lx = FILM_SIZE[1], Ly = FILM_SIZE[2], Lz = FILM_SIZE[3]),
]

function build_interface(n_quad::Int, edge_refine_level::Int)
    l_ec = FILM_SIZE[3] / 2.0^edge_refine_level * 1.01
    interface = BI.multi_dielectric_box3d(n_quad, l_ec, BOXES, EPSES, EPS_OUT)
    return interface
end

function evaluate_scattered_potential(interface, sigma::AbstractVector)
    ntargets = length(TARGETS)
    targets_3xn = targets_matrix_3xn(TARGETS)

    pot = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(
        interface,
        targets_3xn,
        FMM_TOL,
        FMM_TOL,
        5.0;
        include_edges_src = false,
    )
    return vec(pot * sigma)
end

function run_configuration(n_quad::Int, edge_refine_level::Int)
    @printf("\nRunning configuration n_quad=%d, edge_refine=%d\n", n_quad, edge_refine_level)

    interface = build_interface(n_quad, edge_refine_level)
    n_pts = BI.num_points(interface)
    @printf("  interface points: %d\n", n_pts)

    lhs = BI.lhs_dielectric_box3d_fmm3d_corrected(interface, FMM_TOL, FMM_TOL, MAX_ORDER)
    rhs = BI.rhs_dielectric_box3d(interface, PS, EPS_SRC)
    sigma, stats = Krylov.gmres(lhs, rhs; rtol = GMRES_RTOL, atol = GMRES_RTOL, verbose = 1)

    residual = norm(lhs * sigma - rhs) / max(norm(rhs), eps(Float64))
    n_iter = stats.niter
    @printf("  gmres iterations: %d, residual: %.6e\n", n_iter, residual)

    values = evaluate_scattered_potential(interface, vec(sigma))
    for (label, value) in zip(TARGET_LABELS, values)
        @printf("  %s = %.12e\n", label, value)
    end

    (
        n_quad = n_quad,
        edge_refine = edge_refine_level,
        n_pts = n_pts,
        n_iter = n_iter,
        residual = residual,
        U_00 = Float64(values[1]),
        U_01 = Float64(values[2]),
        U_02 = Float64(values[3]),
        U_03 = Float64(values[4]),
    )
end

function print_summary(results::DataFrame)
    if nrow(results) == 0
        println("\nNo completed configurations to summarize.")
        return nothing
    end

    ordered = sorted_results(results)
    finest_pair = last(CONFIGS)
    ref_rows = ordered[(ordered.n_quad .== finest_pair[1]) .& (ordered.edge_refine .== finest_pair[2]), :]
    ref_row = nrow(ref_rows) > 0 ? ref_rows[1, :] : ordered[end, :]

    println("\nSummary relative to finest available result:")
    @printf("  reference = (n_quad=%d, edge_refine=%d)\n", Int(ref_row.n_quad), Int(ref_row.edge_refine))
    @printf(
        "%6s %6s %8s %8s %12s %15s %11s %15s %11s %15s %11s %15s %11s\n",
        "nq",
        "eref",
        "n_pts",
        "n_iter",
        "residual",
        "U_00",
        "err_00",
        "U_01",
        "err_01",
        "U_02",
        "err_02",
        "U_03",
        "err_03",
    )

    for row in eachrow(ordered)
        err_00 = abs(row.U_00 - ref_row.U_00)
        err_01 = abs(row.U_01 - ref_row.U_01)
        err_02 = abs(row.U_02 - ref_row.U_02)
        err_03 = abs(row.U_03 - ref_row.U_03)

        @printf(
            "%6d %6d %8d %8d %12.3e %15.8e %11.3e %15.8e %11.3e %15.8e %11.3e %15.8e %11.3e\n",
            Int(row.n_quad),
            Int(row.edge_refine),
            Int(row.n_pts),
            Int(row.n_iter),
            Float64(row.residual),
            Float64(row.U_00),
            err_00,
            Float64(row.U_01),
            err_01,
            Float64(row.U_02),
            err_02,
            Float64(row.U_03),
            err_03,
        )
    end

    nothing
end

function main()
    mkpath(OUTPUT_DIR)
    results = load_results(OUTPUT_PATH)
    done = completed_configs(results)
    failures = Tuple{Int, Int}[]

    if nrow(results) > 0
        println("Loaded existing results from $OUTPUT_PATH")
    end

    for (n_quad, edge_refine_level) in CONFIGS
        config = (n_quad, edge_refine_level)

        if config in done
            @printf("\nSkipping existing configuration n_quad=%d, edge_refine=%d\n", n_quad, edge_refine_level)
            continue
        end

        try
            row = run_configuration(n_quad, edge_refine_level)
            push!(results, row)
            persist_results(results, OUTPUT_PATH)
            push!(done, config)
            println("  saved partial results to $OUTPUT_PATH")
        catch err
            push!(failures, config)
            persist_results(results, OUTPUT_PATH)
            @error "Configuration failed" n_quad edge_refine_level exception = (err, catch_backtrace())
        end
    end

    print_summary(results)

    if !isempty(failures)
        println("\nFailed configurations:")
        for (n_quad, edge_refine_level) in failures
            @printf("  (%d, %d)\n", n_quad, edge_refine_level)
        end
    end

    println("\nResults written to $OUTPUT_PATH")
    nothing
end

main()

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Krylov
using CSV
using DataFrames
using Printf
using Logging

const config_ref = (4, 5)
const ps = [3, 4]
const rs = [1, 2, 3, 4]

const SUBSTRATE_CENTER = (0.0, 0.0, -5.0)
const SUBSTRATE_SIZE = (10.0, 10.0, 10.0)
const FILM_CENTER = (0.0, 0.0, 0.5)
const FILM_SIZE = (5.0, 5.0, 1.0)

const EPSES = [10.0, 2.4]
const EPS_OUT = 1.0
const EPS_SRC = 2.4

const PS = BI.PointSource((0.2, 0.3, 0.4), 1.0)
const FMM_TOL = 1.0e-6
const MAX_ORDER = 128
const GMRES_RTOL = 1.0e-4

const OUTPUT_DIR = "data"
const OUTPUT_PATH = joinpath(OUTPUT_DIR, "singlecharge_convergence.csv")

const TARGETS = [(randn(), randn(), 0.5 + 0.5 * randn()) for _ in 1:100]

const BOXES = [
    (center = SUBSTRATE_CENTER, Lx = SUBSTRATE_SIZE[1], Ly = SUBSTRATE_SIZE[2], Lz = SUBSTRATE_SIZE[3]),
    (center = FILM_CENTER, Lx = FILM_SIZE[1], Ly = FILM_SIZE[2], Lz = FILM_SIZE[3]),
]

function targets_matrix_3xn(targets)
    mat = Matrix{Float64}(undef, 3, length(targets))
    for (j, trg) in enumerate(targets)
        mat[:, j] .= trg
    end
    mat
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

    l_ec = 1.0 / 2^edge_refine_level
    interface = BI.multi_dielectric_box3d_rhs_adaptive(n_quad, l_ec, BOXES, EPSES, PS, EPS_SRC, FMM_TOL)
    n_pts = BI.num_points(interface)
    @printf("  interface points: %d\n", n_pts)

    lhs = BI.lhs_dielectric_box3d_fmm3d_corrected(interface, FMM_TOL, FMM_TOL, MAX_ORDER, include_edges_src = false, include_edges_trg = false)
    rhs = BI.rhs_dielectric_box3d(interface, PS, EPS_SRC)
    sigma, stats = Krylov.gmres(lhs, rhs; rtol = GMRES_RTOL, atol = GMRES_RTOL, verbose = 1)

    residual = norm(lhs * sigma - rhs) / max(norm(rhs), eps(Float64))
    n_iter = stats.niter
    @printf("  gmres iterations: %d, residual: %.6e\n", n_iter, residual)

    pot = evaluate_scattered_potential(interface, vec(sigma))

    return n_quad, edge_refine_level, n_pts, n_iter, residual, pot
end

function main()
    df = CSV.write(OUTPUT_PATH, DataFrame(
        n_quad = Int[],
        edge_refine = Int[],
        n_pts = Int[],
        n_iter = Int[],
        residual = Float64[],
        l2_error = Float64[],
    ))

    # run the reference configuration
    _, _, _, _, _, ref_pot = run_configuration(config_ref...)
    @show ref_pot

    for n_quad in ps
        for edge_refine in rs
            _, _, n_pts, n_iter, residual, pot = run_configuration(n_quad, edge_refine)
            @show n_quad, edge_refine, n_pts, n_iter, residual
            l2_err = norm(pot .- ref_pot) / norm(ref_pot)
            @printf("  L2 error compared to reference: %.6e\n", l2_err)
            CSV.write(df, DataFrame(
                n_quad = n_quad,
                edge_refine = edge_refine,
                n_pts = n_pts,
                n_iter = n_iter,
                residual = residual,
                l2_error = l2_err,
            ); append = true)
        end
    end

    nothing
end

main()
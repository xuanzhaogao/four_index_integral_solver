# timing_density_solve.jl
#
# Per-step timing harness for the *actual* monolayer screened density solve
# (the onsite + nn channels), mirroring screened_monolayer_k323201.jl exactly:
#   - real Wannier orbital density (squared phi_1) as the continuous source
#   - real slab geometry L = 90 A, Lz = 3.35 A, eps_in = 2.4, Sharp screening
#
# It inlines the body of ScreenedOrbitalSolve.solve_screened_mode so each
# stage can be timed individually:
#   1. screened_volume_source        (RHS density prep)
#   2. interface build               (single_dielectric_box3d_rhs_adaptive)
#   3. RHS assembly                  (Rhs_dielectric_box3d_fmm3d)
#   4. LHS assembly  <-- near-correction upsampling lives here
#                                    (Lhs_dielectric_box3d_fmm3d_corrected)
#   5. GMRES solve
#   6. post-processing per target    (volume potential + scattered potential)
#
# A small warm-up solve is run first to remove JIT-compile time, so the
# reported numbers are warm (steady-state) costs. Source loading is timed
# separately (one-time, includes XSF read + compile).

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "..", "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve
include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

# ---------------------------------------------------------------------------
# Parameters -- identical to screened_monolayer_k323201.jl
# ---------------------------------------------------------------------------
const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1   = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2   = joinpath(REF_DIR, "graphene_00002.xsf")

const LZ        = 3.35
const Z_CENTER  = LZ / 2
const EPS_IN    = 2.4
const EPS_OUT   = 1.0
const L         = 90.0

const SOURCE_TOL = 1e-3
const N_QUAD     = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL    = 1e-3
const LHS_TOL    = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER  = 64
const MAX_DEPTH  = 12
const VOLUME_TOL = RHS_TOL
const MODE       = BI.SharpScreening()

# ---------------------------------------------------------------------------
# Instrumented copy of solve_screened_mode's body.
# Returns (times::Vector{Pair{String,Float64}}, n_pts, residual, pair_lines).
# ---------------------------------------------------------------------------
function timed_solve(source_vs, target_specs, Lx, Ly, Lz, eps_in, eps_out, mode;
                     n_quad, edge_refine_level, rhs_tol, lhs_tol,
                     gmres_atol, gmres_rtol, max_order, max_depth, volume_tol,
                     verbose::Bool)
    times = Pair{String,Float64}[]
    l_ec  = Float64(Lz) / 2.0^edge_refine_level * 1.01

    t = @elapsed screened_vs = BI.screened_volume_source(
        Float64(Lx), Float64(Ly), Float64(Lz), source_vs,
        Float64(eps_in), Float64(eps_out), mode; tol = Float64(rhs_tol))
    push!(times, "1. screened_volume_source" => t)

    t = @elapsed resolved_kmax = BI._estimate_tkm3dc_kmax(screened_vs)
    push!(times, "2. estimate tkm kmax" => t)

    t = @elapsed interface = BI.single_dielectric_box3d_rhs_adaptive(
        Float64(Lx), Float64(Ly), Float64(Lz), Int(n_quad), screened_vs,
        1.0, l_ec, Float64(rhs_tol), Float64(eps_in), Float64(eps_out), Float64;
        max_depth = Int(max_depth), tkm_kmax = resolved_kmax)
    push!(times, "3. interface build (adaptive box)" => t)

    t = @elapsed rhs = BI.Rhs_dielectric_box3d_fmm3d(interface, screened_vs, 1.0, Float64(rhs_tol))
    push!(times, "4. RHS assembly (FMM)" => t)

    t = @elapsed lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(
        interface, Float64(lhs_tol), Float64(lhs_tol), Int(max_order);
        correct_edges = false)
    push!(times, "5. LHS assembly (near-correction, correct_edges=false)" => t)

    local sigma, status
    t = @elapsed begin
        sigma, status = Krylov.gmres(lhs, rhs; atol = Float64(gmres_atol),
                                      rtol = Float64(gmres_rtol), verbose = 0)
    end
    push!(times, "6. GMRES solve" => t)

    residual = norm(lhs * sigma - rhs) / max(norm(rhs), eps(Float64))

    pair_lines = String[]
    for spec in target_specs
        tpos = spec.target_vs.positions
        tw   = spec.target_vs.weights .* spec.target_vs.density
        tmat = Matrix{Float64}(tpos)
        tv = @elapsed u_int = ScreenedOrbitalSolve.evaluate_volume_potential(
            screened_vs, tmat; tol = volume_tol, kmax = resolved_kmax)
        # call the FMM-corrected operator directly: the ScreenedOrbitalSolve
        # wrapper forwards include_edges_src, which the current BI branch dropped.
        ts = @elapsed begin
            scatter_op = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(
                interface, tmat, Float64(lhs_tol), Float64(lhs_tol), 5.0)
            u_scatter = scatter_op * sigma
        end
        push!(times, "7. post: $(spec.pair) volume potential (TKM3D)" => tv)
        push!(times, "8. post: $(spec.pair) scattered potential (FMM)" => ts)
        u_raw = dot(u_int, tw) + dot(u_scatter, tw)
        push!(pair_lines, @sprintf("    %-8s u_total_raw = %.6e   (Nt=%d)", spec.pair, u_raw, length(tw)))
    end

    return (times = times, n_pts = BI.num_points(interface),
            residual = residual, pair_lines = pair_lines,
            n_src = length(screened_vs.density))
end

# ---------------------------------------------------------------------------
function main()
    @info "Loading real Wannier orbital sources (113 MB XSF read; cold, one-time)"
    t_load = @elapsed src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER)
    @printf("  source load: %.2f s   (Nphi1=%.4f, Nphi2=%.4f, src pts=%d)\n",
            t_load, src.Nphi1, src.Nphi2, length(src.vs1.density))

    specs = density_target_specs(src)   # onsite -> vs1, nn -> vs2

    # ---- warm-up at tiny size to compile every stage ----
    @info "Warm-up solve (tiny box) to trigger JIT compilation ..."
    gsrc = BI.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 8, 1e-3)
    warm_specs = [(pair = "warm", target_vs = gsrc, Na = 1.0, Nb = 1.0)]
    twarm = @elapsed timed_solve(gsrc, warm_specs, 5.0, 5.0, 1.0, EPS_IN, EPS_OUT, MODE;
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL, rhs_tol = RHS_TOL,
        lhs_tol = LHS_TOL, gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH, volume_tol = VOLUME_TOL, verbose = false)
    @printf("  warm-up wall time: %.2f s\n", twarm)

    # ---- real, warm, timed run ----
    @info "Timed density solve  L=$L  Lz=$LZ  eps_in=$EPS_IN  Sharp screening"
    t_total = @elapsed res = timed_solve(src.vs1, specs, L, L, LZ, EPS_IN, EPS_OUT, MODE;
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL, rhs_tol = RHS_TOL,
        lhs_tol = LHS_TOL, gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH, volume_tol = VOLUME_TOL, verbose = true)

    println()
    println("="^64)
    @printf("  interface points : %d\n", res.n_pts)
    @printf("  screened src pts  : %d\n", res.n_src)
    @printf("  GMRES residual    : %.3e\n", res.residual)
    println("-"^64)
    solve_sum = 0.0
    for (label, secs) in res.times
        @printf("  %-44s %8.2f s\n", label, secs)
        solve_sum += secs
    end
    println("-"^64)
    @printf("  %-44s %8.2f s\n", "sum of timed steps", solve_sum)
    @printf("  %-44s %8.2f s\n", "total solve wall (warm)", t_total)
    @printf("  %-44s %8.2f s\n", "one-time source load (cold)", t_load)
    println("="^64)
    for ln in res.pair_lines
        println(ln)
    end
end

main()

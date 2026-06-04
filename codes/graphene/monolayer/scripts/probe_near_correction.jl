# probe_near_correction.jl
#
# Split the LHS near-correction cost (the 143 s step) into its two parts and
# inspect whether the SourceCache optimization is doing its job:
#
#   A. build_neighbor_list      -- the DECISION step (check_quad_order3d
#                                  self-convergence loop, NOT optimized)
#   B. laplace3d_DT_corrections -- the ASSEMBLY step (SourceCache: builds the
#                                  moments-to-nodal tensor once per source panel)
#
# Also reports the distribution of per-pair upsampling orders and quantifies
# the "n_up_max per source panel" inflation: the assembly integrates every
# neighbor of a source at that source's WORST (largest) order, not the
# per-pair order.

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "..", "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve
include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")
const LZ = 3.35; const Z_CENTER = LZ / 2
const EPS_IN = 2.4; const EPS_OUT = 1.0; const L = 90.0
const SOURCE_TOL = 1e-3; const N_QUAD = 6; const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3; const LHS_TOL = 1e-5; const MAX_ORDER = 64; const MAX_DEPTH = 12

build_interface(src_vs, Lx, Ly, Lz) = begin
    l_ec = Float64(Lz) / 2.0^EDGE_REFINE_LEVEL * 1.01
    kmax = BI._estimate_tkm3dc_kmax(src_vs)
    sv = BI.screened_volume_source(Float64(Lx), Float64(Ly), Float64(Lz), src_vs,
        EPS_IN, EPS_OUT, BI.SharpScreening(); tol = RHS_TOL)
    BI.single_dielectric_box3d_rhs_adaptive(Float64(Lx), Float64(Ly), Float64(Lz),
        N_QUAD, sv, 1.0, l_ec, RHS_TOL, EPS_IN, EPS_OUT, Float64;
        max_depth = MAX_DEPTH, tkm_kmax = kmax)
end

function probe(interface; label::String)
    adaptive_cfg = BI.AdaptiveConfig(LHS_TOL, sqrt(eps(Float64)), 0, 20)

    t_nl = @elapsed nl = BI.build_neighbor_list(interface, MAX_ORDER, LHS_TOL;
        range_factor = 5.0, correct_edges = true, adaptive_cfg = adaptive_cfg)
    upsample = nl.upsample; adaptive = nl.adaptive

    t_asm = @elapsed corr = BI.laplace3d_DT_corrections(interface, upsample, adaptive)

    # --- per-pair order histogram ---
    orders = sort(collect(values(upsample)))
    hist = Dict{Int,Int}()
    for o in orders; hist[o] = get(hist, o, 0) + 1; end

    # --- n_up_max per source panel + work inflation ---
    src_to = Dict{Int, Vector{Tuple{Int,Int}}}()
    for ((i, j), nup) in upsample
        push!(get!(() -> Tuple{Int,Int}[], src_to, i), (j, nup))
    end
    npts(p) = length(interface.panels[p].points)
    work_pairwise = 0.0   # Σ_pairs np_trg * n_up_pair^2   (ideal, per-pair order)
    work_actual   = 0.0   # Σ_pairs np_trg * n_up_max(src)^2 (what SourceCache does)
    for (i, lst) in src_to
        nmax = maximum(n for (_, n) in lst)
        for (j, nup) in lst
            work_pairwise += npts(j) * nup^2
            work_actual   += npts(j) * nmax^2
        end
    end

    println("\n", "="^60, "  ", label)
    @printf("  interface points       : %d\n", BI.num_points(interface))
    @printf("  panels                  : %d\n", length(interface.panels))
    @printf("  upsample pairs          : %d\n", length(upsample))
    @printf("  adaptive pairs          : %d\n", length(adaptive))
    println("-"^60)
    @printf("  A. build_neighbor_list  (decision) : %8.2f s\n", t_nl)
    @printf("  B. laplace3d_DT_corrections (assy) : %8.2f s\n", t_asm)
    @printf("     near-correction total           : %8.2f s\n", t_nl + t_asm)
    println("-"^60)
    println("  per-pair upsample order histogram:")
    for o in sort(collect(keys(hist)))
        @printf("     n_up = %3d : %6d pairs\n", o, hist[o])
    end
    @printf("  source panels with upsample nbrs   : %d\n", length(src_to))
    @printf("  kernel-eval work, per-pair order   : %.3e\n", work_pairwise)
    @printf("  kernel-eval work, n_up_max (actual): %.3e\n", work_actual)
    @printf("  inflation from n_up_max grouping   : %.2fx\n", work_actual / work_pairwise)
    println("="^60)
end

function main()
    @info "loading sources"
    src = monolayer_screened_sources(orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER)

    @info "warm-up (tiny interface) to compile probe path"
    g = BI.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 8, 1e-3)
    probe(build_interface(g, 5.0, 5.0, 1.0); label = "WARMUP (discard)")

    @info "building real interface (L=$L, Lz=$LZ)"
    t_if = @elapsed interface = build_interface(src.vs1, L, L, LZ)
    @printf("  interface build: %.2f s\n", t_if)

    probe(interface; label = "REAL  L=90 Lz=3.35")
end

main()

# exp_route_both_edge_adaptive.jl
#
# Experiment: route the "both panels is_edge" near-correction pairs from the
# UPSAMPLE path (uniform order-64) to the ADAPTIVE quadtree path, and measure
# (a) near-correction assembly time and (b) solution accuracy.
#
# Fast surface: a POINT source at the slab center drives single_dielectric_
# box3d_rhs_adaptive, so the box + geometric edge bands build in seconds
# (no 355k-point volume source). The edge-band geometry -- which generates
# the both-edge pairs -- is reproduced.
#
# Accuracy check: Gauss-law identity for a point charge of strength 1 embedded
# in the box dielectric:   dot(all_weights, sigma) + 1/eps_in  ==  1.
# Reported for baseline, modified routing, and a tight-tol "gold" (adaptive
# both-edge at atol=1e-9) so we can tell which routing is actually more
# accurate, not just whether the two agree.

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov, LinearAlgebra, Printf

const LZ = 1.0; const EPS_IN = 2.4; const EPS_OUT = 1.0; const L = 10.0
const N_QUAD = 6; const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3; const LHS_TOL = 1e-5; const FMM_TOL = 1e-5
const MAX_ORDER = 64; const MAX_DEPTH = 12
const GTOL = 1e-7

const l_ec = LZ / 2.0^EDGE_REFINE_LEVEL * 1.01

# dielectric LHS from arbitrary correction operator (mirrors
# lhs_dielectric_box3d_fmm3d_corrected, but with injectable corrections)
function make_lhs(interface, D_base, corr)
    n = BI.num_points(interface)
    panels = interface.panels
    function apply(charges::AbstractVector{Float64})
        y = D_base * charges
        y .+= corr * charges
        off = 0
        @inbounds for i in 1:length(panels)
            ei = interface.eps_in[i]; eo = interface.eps_out[i]
            t = 0.5 * (eo + ei) / (eo - ei)
            np = BI.num_points(panels[i])
            for j in 1:np
                y[off + j] += t * charges[off + j]
            end
            off += np
        end
        return y
    end
    return BI.LinearMap{Float64}(apply, n, n)
end

gauss_err(ws, x) = abs(dot(ws, x) + 1.0 / EPS_IN - 1.0)

function main()
    ps = BI.PointSource((0.0, 0.0, 0.0), 1.0)

    @info "warmup (tiny box)"
    let psw = BI.PointSource((0.0,0.0,0.0), 1.0)
        ifw = BI.single_dielectric_box3d_rhs_adaptive(5.0,5.0,1.0, N_QUAD, psw, EPS_IN,
            1.0/2^EDGE_REFINE_LEVEL*1.01, RHS_TOL, EPS_IN, EPS_OUT, Float64; max_depth=4)
        cfgw = BI.AdaptiveConfig(LHS_TOL, sqrt(eps(Float64)), 0, 20)
        nlw = BI.build_neighbor_list(ifw, MAX_ORDER, LHS_TOL; range_factor=5.0, correct_edges=true, adaptive_cfg=cfgw)
        BI.laplace3d_DT_corrections(ifw, nlw.upsample, nlw.adaptive)
        BI.laplace3d_DT_fmm3d(ifw, FMM_TOL)
    end

    @info "building point-source interface  L=$L Lz=$LZ"
    t_if = @elapsed interface = BI.single_dielectric_box3d_rhs_adaptive(
        L, L, LZ, N_QUAD, ps, EPS_IN, l_ec, RHS_TOL, EPS_IN, EPS_OUT, Float64; max_depth=MAX_DEPTH)
    @printf("  interface build: %.2f s   panels=%d  points=%d\n",
            t_if, length(interface.panels), BI.num_points(interface))

    rhs = BI.rhs_dielectric_box3d(interface, ps, EPS_IN)
    ws  = BI.all_weights(interface)
    panels = interface.panels

    t_dbase = @elapsed D_base = BI.laplace3d_DT_fmm3d(interface, FMM_TOL)
    @printf("  far FMM base build: %.2f s\n", t_dbase)

    cfg = BI.AdaptiveConfig(LHS_TOL, sqrt(eps(Float64)), 0, 20)
    cfg_gold = BI.AdaptiveConfig(1e-9, sqrt(eps(Float64)), 0, 24)

    t_nl = @elapsed nl = BI.build_neighbor_list(interface, MAX_ORDER, LHS_TOL;
        range_factor=5.0, correct_edges=true, adaptive_cfg=cfg)
    upsample = nl.upsample; adaptive = nl.adaptive
    @printf("  build_neighbor_list: %.2f s   upsample=%d adaptive=%d\n",
            t_nl, length(upsample), length(adaptive))

    # ---- modified routing: both-edge upsample pairs -> adaptive ----
    moved = 0
    ups_mod = Dict{Tuple{Int,Int},Int}()
    adp_mod = copy(adaptive)
    ups_gold = Dict{Tuple{Int,Int},Int}()
    adp_gold = Dict{Tuple{Int,Int},BI.AdaptiveConfig}()
    for (k,v) in adaptive; adp_gold[k] = cfg_gold; end
    for ((i,j),nup) in upsample
        if panels[i].is_edge && panels[j].is_edge
            adp_mod[(i,j)]  = cfg
            adp_gold[(i,j)] = cfg_gold
            moved += 1
        else
            ups_mod[(i,j)]  = nup
            ups_gold[(i,j)] = nup
        end
    end
    @printf("  both-edge pairs moved to adaptive: %d\n", moved)
    @printf("  -> baseline:  upsample=%d adaptive=%d\n", length(upsample), length(adaptive))
    @printf("  -> modified:  upsample=%d adaptive=%d\n", length(ups_mod), length(adp_mod))

    # ---- assemble corrections (this is the cost we care about) ----
    t_base = @elapsed corr_base = BI.laplace3d_DT_corrections(interface, upsample, adaptive)
    t_mod  = @elapsed corr_mod  = BI.laplace3d_DT_corrections(interface, ups_mod, adp_mod)
    t_gold = @elapsed corr_gold = BI.laplace3d_DT_corrections(interface, ups_gold, adp_gold)

    # ---- solve all three ----
    lhs_base = make_lhs(interface, D_base, corr_base)
    lhs_mod  = make_lhs(interface, D_base, corr_mod)
    lhs_gold = make_lhs(interface, D_base, corr_gold)
    x_base,_ = Krylov.gmres(lhs_base, rhs; atol=GTOL, rtol=GTOL)
    x_mod, _ = Krylov.gmres(lhs_mod,  rhs; atol=GTOL, rtol=GTOL)
    x_gold,_ = Krylov.gmres(lhs_gold, rhs; atol=GTOL, rtol=GTOL)

    println("\n", "="^60)
    @printf("  near-correction assembly time:\n")
    @printf("     baseline (upsample order-64) : %8.2f s\n", t_base)
    @printf("     modified (both-edge->adaptive): %8.2f s   (%.2fx)\n", t_mod, t_base/t_mod)
    @printf("     gold     (both-edge->adp 1e-9): %8.2f s\n", t_gold)
    println("-"^60)
    @printf("  Gauss-law error |dot(ws,x)+1/eps-1|:\n")
    @printf("     baseline : %.3e\n", gauss_err(ws, x_base))
    @printf("     modified : %.3e\n", gauss_err(ws, x_mod))
    @printf("     gold     : %.3e\n", gauss_err(ws, x_gold))
    println("-"^60)
    @printf("  ||x_mod  - x_gold|| / ||x_gold|| : %.3e\n", norm(x_mod  - x_gold)/norm(x_gold))
    @printf("  ||x_base - x_gold|| / ||x_gold|| : %.3e\n", norm(x_base - x_gold)/norm(x_gold))
    @printf("  ||x_mod  - x_base|| / ||x_base|| : %.3e\n", norm(x_mod  - x_base)/norm(x_base))
    println("="^60)
end

main()

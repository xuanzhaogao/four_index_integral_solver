# probe_pair_isedge.jl
# For each near-correction upsample pair, check whether the panels are flagged
# is_edge (refined edge-band panels). Tally both/one/neither, all pairs and the
# n_up=64 capped subset.

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "..", "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve
include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf"); const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")
const LZ = 3.35; const Z_CENTER = LZ/2; const EPS_IN = 2.4; const EPS_OUT = 1.0; const L = 90.0
const SOURCE_TOL = 1e-3; const N_QUAD = 6; const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3; const LHS_TOL = 1e-5; const MAX_ORDER = 64; const MAX_DEPTH = 12

build_interface(src_vs, Lx, Ly, Lz) = begin
    l_ec = Float64(Lz)/2.0^EDGE_REFINE_LEVEL*1.01
    kmax = BI._estimate_tkm3dc_kmax(src_vs)
    sv = BI.screened_volume_source(Float64(Lx),Float64(Ly),Float64(Lz),src_vs,EPS_IN,EPS_OUT,BI.SharpScreening();tol=RHS_TOL)
    BI.single_dielectric_box3d_rhs_adaptive(Float64(Lx),Float64(Ly),Float64(Lz),N_QUAD,sv,1.0,l_ec,RHS_TOL,EPS_IN,EPS_OUT,Float64;max_depth=MAX_DEPTH,tkm_kmax=kmax)
end

function main()
    src = monolayer_screened_sources(orbital_1=XSF_1, orbital_2=XSF_2, source_tol=SOURCE_TOL, z_center=Z_CENTER)
    build_interface(BI.GaussianVolumeSource((0.,0.,0.),0.3,8,1e-3),5.,5.,1.)  # warmup
    @info "building real interface"
    interface = build_interface(src.vs1, L, L, LZ)
    cfg = BI.AdaptiveConfig(LHS_TOL, sqrt(eps(Float64)), 0, 20)
    nl = BI.build_neighbor_list(interface, MAX_ORDER, LHS_TOL; range_factor=5.0, correct_edges=true, adaptive_cfg=cfg)
    upsample = nl.upsample
    panels = interface.panels

    n_edge_panels = count(p -> p.is_edge, panels)
    both=0; one=0; neither=0
    capboth=0; capone=0; capneither=0
    for ((i,j),nup) in upsample
        ei = panels[i].is_edge; ej = panels[j].is_edge
        n = ei + ej
        if n==2; both+=1; elseif n==1; one+=1; else; neither+=1; end
        if nup==MAX_ORDER
            if n==2; capboth+=1; elseif n==1; capone+=1; else; capneither+=1; end
        end
    end
    println("\n", "="^56)
    @printf("  total panels        : %d\n", length(panels))
    @printf("  is_edge panels      : %d  (%.1f%%)\n", n_edge_panels, 100*n_edge_panels/length(panels))
    @printf("  upsample pairs      : %d\n", length(upsample))
    println("-"^56)
    println("  ALL upsample pairs, by is_edge membership:")
    @printf("     both panels edge    : %6d  (%.1f%%)\n", both, 100*both/length(upsample))
    @printf("     exactly one edge    : %6d  (%.1f%%)\n", one, 100*one/length(upsample))
    @printf("     neither edge        : %6d  (%.1f%%)\n", neither, 100*neither/length(upsample))
    println("-"^56)
    ncap = capboth+capone+capneither
    @printf("  capped n_up=64 pairs (%d):\n", ncap)
    @printf("     both panels edge    : %6d  (%.1f%%)\n", capboth, 100*capboth/max(ncap,1))
    @printf("     exactly one edge    : %6d  (%.1f%%)\n", capone, 100*capone/max(ncap,1))
    @printf("     neither edge        : %6d  (%.1f%%)\n", capneither, 100*capneither/max(ncap,1))
    println("="^56)
end
main()

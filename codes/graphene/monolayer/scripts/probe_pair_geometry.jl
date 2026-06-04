# probe_pair_geometry.jl
#
# Classify the near-correction upsample pairs by geometry to understand WHY
# there are so many. For each upsample pair (i,j) we record which box face the
# source and target panels belong to (by normal direction: +z,-z,+x,-x,+y,-y),
# the center-to-center distance, and the panel sizes. We then tally:
#   - face-pair combinations (e.g. +z<->-z = across-slab, +z<->+z = coplanar)
#   - the same breakdown restricted to the n_up=64 (max-order-capped) pairs
#   - distance / panel-half-size ratio (the "d" of the fig4 Bernstein picture)

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

# label a panel by its outward normal -> face name
function face_of(panel)
    n = panel.normal
    ax = argmax(abs.(n))
    s = n[ax] > 0 ? "+" : "-"
    return s * ("xyz"[ax])
end

panel_half(panel) = begin
    a, b, c, d = panel.corners
    max(norm(b .- a), norm(d .- a)) / 2
end
panel_center(panel) = begin
    a, b, c, d = panel.corners
    (a .+ b .+ c .+ d) ./ 4
end

function main()
    @info "loading sources"
    src = monolayer_screened_sources(orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER)

    @info "warmup"; build_interface(BI.GaussianVolumeSource((0.,0.,0.),0.3,8,1e-3), 5.,5.,1.)
    @info "building real interface"
    interface = build_interface(src.vs1, L, L, LZ)
    adaptive_cfg = BI.AdaptiveConfig(LHS_TOL, sqrt(eps(Float64)), 0, 20)
    nl = BI.build_neighbor_list(interface, MAX_ORDER, LHS_TOL;
        range_factor = 5.0, correct_edges = true, adaptive_cfg = adaptive_cfg)
    upsample = nl.upsample

    panels = interface.panels
    faces  = [face_of(p) for p in panels]
    halves = [panel_half(p) for p in panels]
    cents  = [panel_center(p) for p in panels]

    # face-pair tally (unordered), all upsample pairs and capped-only
    facepair(i,j) = join(sort([faces[i], faces[j]]), " <-> ")
    tally_all = Dict{String,Int}()
    tally_cap = Dict{String,Int}()
    dser = Float64[]   # center-distance / mean-half-size for all pairs
    same_face_center_cluster = 0
    for ((i,j), nup) in upsample
        fp = facepair(i,j)
        tally_all[fp] = get(tally_all, fp, 0) + 1
        if nup == MAX_ORDER
            tally_cap[fp] = get(tally_cap, fp, 0) + 1
        end
        dist = norm(cents[i] .- cents[j])
        hbar = (halves[i] + halves[j]) / 2
        push!(dser, dist / hbar)
    end

    @info "panel size range" min_half=minimum(halves) max_half=maximum(halves)
    # how many panels per face
    fcount = Dict{String,Int}()
    for f in faces; fcount[f] = get(fcount,f,0)+1; end

    println("\n", "="^60)
    @printf("  total panels: %d   upsample pairs: %d\n", length(panels), length(upsample))
    println("  panels per face:")
    for f in sort(collect(keys(fcount))); @printf("     %-3s : %6d\n", f, fcount[f]); end
    println("-"^60)
    println("  upsample pairs by face-pair (ALL):")
    for fp in sort(collect(keys(tally_all)); by = x->-tally_all[x])
        @printf("     %-12s : %6d  (%.1f%%)\n", fp, tally_all[fp], 100*tally_all[fp]/length(upsample))
    end
    println("-"^60)
    ncap = sum(values(tally_cap))
    @printf("  capped (n_up=64) pairs by face-pair  (total %d):\n", ncap)
    for fp in sort(collect(keys(tally_cap)); by = x->-tally_cap[x])
        @printf("     %-12s : %6d  (%.1f%%)\n", fp, tally_cap[fp], 100*tally_cap[fp]/ncap)
    end
    println("-"^60)
    sort!(dser)
    @printf("  center-dist / panel-half ratio:  min=%.2f  median=%.2f  max=%.2f\n",
            dser[1], dser[cld(length(dser),2)], dser[end])
    println("="^60)
end

main()

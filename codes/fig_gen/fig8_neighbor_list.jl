#=
Figure 8: Neighbor-list visualization for a high-aspect-ratio box.

System: Lx × Ly × Lz with edge-corner refinement l_ec, panels carrying
p×p Gauss-Legendre quadrature. We pick a large panel near the center of
the +z face and highlight every panel listed as one of its neighbors by
`BoundaryIntegral.build_neighbor_list` (the upsample dict plus the
corner-touching adaptive dict). Co-planar same-face panels are excluded
by design: they are handled by direct same-plane quadrature, not by the
near-field correction list.
=#

using LinearAlgebra
using CairoMakie
using CairoMakie.Colors: red, green, blue
using GeometryBasics
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# System
# ---------------------------------------------------------------------------
const Lx, Ly, Lz       = 30.0, 30.0, 1.0
const n_quad           = 8
const l_ec             = 0.1
const eps_in, eps_out  = 10.0, 1.0
const range_factor     = 5.0

@info "Building panelization" Lx Ly Lz n_quad l_ec
const interface = single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_ec,
                                          eps_in, eps_out)
const n_panels  = length(interface.panels)
@info "  panels" n_panels

# ---------------------------------------------------------------------------
# Panel geometry helpers
# ---------------------------------------------------------------------------
@inline panel_center(panel) =
    (panel.corners[1] .+ panel.corners[2] .+
     panel.corners[3] .+ panel.corners[4]) ./ 4

@inline function panel_length(panel)
    return max(norm(panel.corners[1] .- panel.corners[2]),
               norm(panel.corners[2] .- panel.corners[3]))
end

@inline function panel_area(panel)
    a, b, c, d = panel.corners
    return norm(b .- a) * norm(d .- a)
end

# ---------------------------------------------------------------------------
# Pick a large panel near the center of the +z face
# ---------------------------------------------------------------------------
top_face_ids = Int[]
for (i, p) in enumerate(interface.panels)
    if all(abs(c[3] - Lz/2) < 1e-12 for c in p.corners)
        push!(top_face_ids, i)
    end
end
@info "  +z face panels" length(top_face_ids)

function pick_focal(panels, ids, target_pt)
    best_idx = 0
    best     = -Inf
    for i in ids
        p   = panels[i]
        c   = panel_center(p)
        d2  = (c[1]-target_pt[1])^2 + (c[2]-target_pt[2])^2
        A   = panel_area(p)
        s   = A / (1 + d2)
        if s > best
            best     = s
            best_idx = i
        end
    end
    return best_idx
end

const focal_idx = pick_focal(interface.panels, top_face_ids,
                             (0.0, 0.0, Lz/2))
const focal     = interface.panels[focal_idx]
const c_focal   = panel_center(focal)
const l_focal   = panel_length(focal)
const r_focal   = range_factor * l_focal / focal.n_quad
@info "  focal panel" focal_idx c_focal l_focal r_focal area=panel_area(focal)

# ---------------------------------------------------------------------------
# Neighbor list from BoundaryIntegral (the actual list used by the solver)
# ---------------------------------------------------------------------------
@info "Building neighbor list" range_factor
const (upsample, adaptive) = BI.build_neighbor_list(interface, 1, 1e-12;
                                                    distance_only = true,
                                                    range_factor   = range_factor,
                                                    correct_edges  = true)
@info "  near pairs (i,j)" length(upsample) length(adaptive)

function collect_neighbors(focal, pair_dicts...)
    s = Set{Int}()
    for d in pair_dicts
        for (i, j) in keys(d)
            if i == focal
                push!(s, j)
            elseif j == focal
                push!(s, i)
            end
        end
    end
    return s
end

const neighbor_set = collect_neighbors(focal_idx, upsample, adaptive)
@info "  neighbors of focal panel" length(neighbor_set)

# ---------------------------------------------------------------------------
# 3D plot
# ---------------------------------------------------------------------------
begin
    fig = Figure(size = (900, 500), fontsize = 18, backgroundcolor = :white)

    ax = Axis3(fig[1, 1];
            aspect      = :data,
            azimuth     = 0.30π,
            elevation   = 0.22π,
            protrusions = (0, 0, 0, 0))
    hidedecorations!(ax)
    hidespines!(ax)

    # Color style consistent with the fig2/fig4/fig5 viridis palette: focal
    # and neighbor panels pulled from a categorical viridis ramp; everything
    # else rendered as a semi-transparent gray skin so cross-face neighbors
    # (which sit on the opposite face of the slab) remain visible.
    const vir            = cgrad(:viridis, 3, categorical = true)
    const focal_color    = RGBAf(red(vir[3]), green(vir[3]), blue(vir[3]), 1.0)
    const neighbor_color = RGBAf(red(vir[2]), green(vir[2]), blue(vir[2]), 1.0)
    const other_color    = RGBAf(0.85, 0.85, 0.85, 0.18)

    verts_all  = Point3f[]
    faces_all  = TriangleFace{Int}[]
    colors_all = RGBAf[]
    edge_segs  = Point3f[]

    for (i, panel) in enumerate(interface.panels)
        a, b, c, d = panel.corners
        base = length(verts_all)
        push!(verts_all, Point3f(a[1], a[2], a[3]))
        push!(verts_all, Point3f(b[1], b[2], b[3]))
        push!(verts_all, Point3f(c[1], c[2], c[3]))
        push!(verts_all, Point3f(d[1], d[2], d[3]))
        push!(faces_all, TriangleFace{Int}(base+1, base+2, base+3))
        push!(faces_all, TriangleFace{Int}(base+1, base+3, base+4))

        col = if i == focal_idx
            focal_color
        elseif i in neighbor_set
            neighbor_color
        else
            other_color
        end
        for _ in 1:4
            push!(colors_all, col)
        end

        push!(edge_segs, Point3f(a[1], a[2], a[3])); push!(edge_segs, Point3f(b[1], b[2], b[3]))
        push!(edge_segs, Point3f(b[1], b[2], b[3])); push!(edge_segs, Point3f(c[1], c[2], c[3]))
        push!(edge_segs, Point3f(c[1], c[2], c[3])); push!(edge_segs, Point3f(d[1], d[2], d[3]))
        push!(edge_segs, Point3f(d[1], d[2], d[3])); push!(edge_segs, Point3f(a[1], a[2], a[3]))
    end

    mesh!(ax, GeometryBasics.Mesh(verts_all, faces_all);
        color = colors_all,
        shading = NoShading, transparency = true, rasterize = 4)
    linesegments!(ax, edge_segs; color = (:black, 0.18), linewidth = 0.25)

    elem_other    = PolyElement(color = RGBAf(0.85, 0.85, 0.85, 1.0),
                                strokecolor = :black, strokewidth = 0.5)
    elem_neighbor = PolyElement(color = neighbor_color,
                                strokecolor = :black, strokewidth = 0.5)
    elem_focal    = PolyElement(color = focal_color,
                                strokecolor = :black, strokewidth = 0.5)
    Legend(fig[1, 2],
        [elem_focal, elem_neighbor],
        ["focal panel", "neighbors ($(length(neighbor_set)))"];
        framevisible = false)

    outpath = joinpath(@__DIR__, "figs/fig8_neighbor_list.pdf")
    save(outpath, fig; px_per_unit = 4)
    png_out = replace(outpath, ".pdf" => ".png")
    save(png_out, fig; px_per_unit = 4)
    @info "Saved figure" outpath png_out

    fig
end
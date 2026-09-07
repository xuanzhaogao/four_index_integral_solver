#=
Figure 8: Neighbor-list visualization for a high-aspect-ratio box.

System: Lx × Ly × Lz with edge-corner refinement l_ec, panels carrying
p×p Gauss-Legendre quadrature. We pick a large panel near the center of
the +z face and highlight every panel listed as one of its neighbors by
`BoundaryIntegral.build_neighbor_list` (the upsample dict plus the
corner-touching adaptive dict). Neighbors are selected by the Bernstein-radius
criterion ρ_min(P,Q) ≤ ρ⋆ = ε^{-1/(2p)} (ε = up_atol). Co-planar same-face panels
are excluded by design: they are handled by direct same-plane quadrature, not by
the near-field correction list.
=#

using LinearAlgebra
using CairoMakie
using GeometryBasics
using BoundaryIntegral
const BI = BoundaryIntegral

include(joinpath(@__DIR__, "fig_style.jl"))

# ---------------------------------------------------------------------------
# System
# ---------------------------------------------------------------------------
const Lx, Ly, Lz       = 30.0, 30.0, 1.0
const n_quad           = 8
const l_ec             = 0.1
const eps_in, eps_out  = 10.0, 1.0
# Zoom window: the (+x,+y) corner, the one nearest the camera at azimuth 0.30pi. It projects to
# bottom-centre rather than to the inset side, so its leaders run further than those of the
# (-x,+y) corner, but it reads best as a corner.
const ZW = 2.5                       # window width in box units
const ZX = (Lx/2 - ZW, Lx/2 + 0.2)
const ZY = (Ly/2 - ZW, Ly/2 + 0.2)
const up_atol          = 1e-6       # tolerance ε setting the Bernstein threshold ρ⋆ = ε^{-1/(2p)}
const max_order        = 64         # cap on the upsample order (selection is independent of it)

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

@inline function in_zoom(panel)
    # Show only the THREE faces that meet at the corner (+z top and the +x/+y sides). The -z
    # bottom face also falls inside the window, but drawing all four reads as a hollow shell
    # rather than a corner.
    all(abs(c[3] + Lz/2) < 1e-12 for c in panel.corners) && return false
    xs = (panel.corners[1][1], panel.corners[2][1], panel.corners[3][1], panel.corners[4][1])
    ys = (panel.corners[1][2], panel.corners[2][2], panel.corners[3][2], panel.corners[4][2])
    # overlap, not centre-containment: a centre test drops panels that straddle the window
    # boundary and leaves holes in the zoom patch.
    return maximum(xs) >= ZX[1] && minimum(xs) <= ZX[2] &&
           maximum(ys) >= ZY[1] && minimum(ys) <= ZY[2]
end

@inline function centre_in_zoom(panel)
    ctr = panel_center(panel)
    return ZX[1] <= ctr[1] <= ZX[2] && ZY[1] <= ctr[2] <= ZY[2]
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

function pick_source(panels, ids, target_pt)
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

const source_idx = pick_source(interface.panels, top_face_ids,
                             (0.0, 0.0, Lz/2))
const source     = interface.panels[source_idx]
const c_source   = panel_center(source)
const l_source   = panel_length(source)
const rho_star  = up_atol ^ (-1 / (2 * source.n_quad))
@info "  source panel" source_idx c_source l_source rho_star area=panel_area(source)

# ---------------------------------------------------------------------------
# Neighbor list from BoundaryIntegral (the actual list used by the solver)
# ---------------------------------------------------------------------------
@info "Building neighbor list (Bernstein ρ_min ≤ ρ⋆ criterion)" up_atol rho_star
const (upsample, adaptive) = BI.build_neighbor_list(interface, max_order, up_atol;
                                                    correct_edges = true)
@info "  near pairs (i,j)" length(upsample) length(adaptive)

function collect_neighbors(source, pair_dicts...)
    s = Set{Int}()
    for d in pair_dicts
        for (i, j) in keys(d)
            if i == source
                push!(s, j)
            elseif j == source
                push!(s, i)
            end
        end
    end
    return s
end

const neighbor_set = collect_neighbors(source_idx, upsample, adaptive)
@info "  neighbors of source panel" length(neighbor_set)

# ---------------------------------------------------------------------------
# 3D plot
# ---------------------------------------------------------------------------
begin
    fig = Figure(size = (FIG_W, FIG_H_3D), backgroundcolor = :white)

    ax = Axis3(fig[1, 2];
            aspect      = :data,
            azimuth     = 0.30π,
            elevation   = 0.22π,
            protrusions = (0, 0, 0, 0))
    hidedecorations!(ax)
    hidespines!(ax)

    # (b) zoom window on one edge: the edge-corner refinement (l_ec) grades the panels
    # down toward the edge, which is only a dark band at print size in panel (a).
    axz = Axis3(fig[1, 3];
            aspect      = :data,
            azimuth     = 0.30π,
            elevation   = 0.22π,
            protrusions = (0, 0, 0, 0))
    hidedecorations!(axz)
    hidespines!(axz)

    # Color style (see fig_style.jl): source and neighbor panels are highlighted
    # with the bright qualitative QUAL palette; every other panel is a
    # semi-transparent gray skin so cross-face neighbors (which sit on the
    # opposite face of the slab) remain visible.
    const source_color    = QUAL.orange
    const neighbor_color = QUAL.blue
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

        col = if i == source_idx
            source_color
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

    # ---- zoom panel ----
    # Build a SEPARATE mesh from only the panels centred inside the window. Reusing the full
    # mesh and relying on limits! does not clip: panels straddling the window boundary are still
    # drawn in full and trail long stray edges across the panel.
    verts_z = Point3f[]; faces_z = TriangleFace{Int}[]; colors_z = RGBAf[]; edge_z = Point3f[]
    for (i, panel) in enumerate(interface.panels)
        in_zoom(panel) || continue
        # Clip the panel to the window. Every panel is an axis-aligned rectangle on a box face,
        # so clipping is an interval intersection per axis. Drawing overlapping panels WHOLE
        # instead leaves a ragged border in which a large protruding panel reads as detached.
        xs = (panel.corners[1][1], panel.corners[2][1], panel.corners[3][1], panel.corners[4][1])
        ys = (panel.corners[1][2], panel.corners[2][2], panel.corners[3][2], panel.corners[4][2])
        zs = (panel.corners[1][3], panel.corners[2][3], panel.corners[3][3], panel.corners[4][3])
        x0 = max(minimum(xs), ZX[1]); x1 = min(maximum(xs), ZX[2])
        y0 = max(minimum(ys), ZY[1]); y1 = min(maximum(ys), ZY[2])
        z0 = minimum(zs);             z1 = maximum(zs)
        (x1 >= x0 && y1 >= y0) || continue
        # `real` marks which of the quad's four sides is an actual panel boundary. A side that
        # the clip created is NOT a panel edge and must not be stroked, or the cut face reads as
        # a panel row that does not exist.
        xlo, xhi = minimum(xs), maximum(xs)
        ylo, yhi = minimum(ys), maximum(ys)
        quad, real = if z0 == z1                 # +z face: varies in x and y
            ((Point3f(x0,y0,z0), Point3f(x1,y0,z0), Point3f(x1,y1,z0), Point3f(x0,y1,z0)),
             (y0 == ylo, x1 == xhi, y1 == yhi, x0 == xlo))
        elseif x0 == x1                          # x-normal side: varies in y and z
            ((Point3f(x0,y0,z0), Point3f(x0,y1,z0), Point3f(x0,y1,z1), Point3f(x0,y0,z1)),
             (true, y1 == yhi, true, y0 == ylo))
        else                                     # y-normal side: varies in x and z
            ((Point3f(x0,y0,z0), Point3f(x1,y0,z0), Point3f(x1,y0,z1), Point3f(x0,y0,z1)),
             (true, x1 == xhi, true, x0 == xlo))
        end
        a, b, c, d = quad
        base = length(verts_z)
        for v in quad
            push!(verts_z, v)
        end
        push!(faces_z, TriangleFace{Int}(base+1, base+2, base+3))
        push!(faces_z, TriangleFace{Int}(base+1, base+3, base+4))
        col = i == source_idx ? source_color : (i in neighbor_set ? neighbor_color : other_color)
        for _ in 1:4
            push!(colors_z, col)
        end
        for (k, (p, q)) in enumerate(((a,b), (b,c), (c,d), (d,a)))
            real[k] || continue
            push!(edge_z, p); push!(edge_z, q)
        end
    end
    @info "  zoom-window panels" length(faces_z) ÷ 2
    # measured panel-size range inside the zoom window, for the caption
    let pl = [panel_length(pn) for pn in interface.panels if centre_in_zoom(pn)]
        @info "  zoom panel_length" min = minimum(pl) max = maximum(pl) ratio = maximum(pl)/minimum(pl)
    end
    mesh!(axz, GeometryBasics.Mesh(verts_z, faces_z);
        color = colors_z,
        shading = NoShading, transparency = true, rasterize = 4)
    # heavier strokes here: at this scale the panel grading is the content of the panel
    linesegments!(axz, edge_z; color = (:black, 0.55), linewidth = 0.5)

    # outline the zoom window on panel (a), on the +z face
    zrect = [Point3f(ZX[1], ZY[1], Lz/2), Point3f(ZX[2], ZY[1], Lz/2),
             Point3f(ZX[2], ZY[2], Lz/2), Point3f(ZX[1], ZY[2], Lz/2),
             Point3f(ZX[1], ZY[1], Lz/2)]
    lines!(ax, zrect; color = :black, linewidth = 1.6)

    elem_other    = PolyElement(color = RGBAf(0.85, 0.85, 0.85, 1.0),
                                strokecolor = :black, strokewidth = 0.5)
    elem_neighbor = PolyElement(color = neighbor_color,
                                strokecolor = :black, strokewidth = 0.5)
    elem_source    = PolyElement(color = source_color,
                                strokecolor = :black, strokewidth = 0.5)
    Legend(fig[1, 1],
        [elem_source, elem_neighbor],
        ["source panel", "neighbors ($(length(neighbor_set)))"];
        framevisible = false)

    colgap!(fig.layout, 4)
    # The inset is a magnifier, not a co-equal panel: give it a much narrower column than the
    # main view. Must precede update_state_before_display! below, since the callout geometry is
    # read back from the resulting viewports.
    colsize!(fig.layout, 2, Relative(0.58))
    colsize!(fig.layout, 3, Relative(0.24))

    # ---- magnifier callout ----
    # The zoom is not a standalone panel: frame it and run leaders back to the outlined window
    # in the main view. Positions come from projecting the window's 3D corners into figure pixel
    # space, so the leaders stay attached if the camera, window, or layout changes.
    Makie.update_state_before_display!(fig)
    vp = ax.scene.viewport[]; vz = axz.scene.viewport[]
    proj = map((Point3f(ZX[1], ZY[1], Lz/2), Point3f(ZX[2], ZY[1], Lz/2),
                Point3f(ZX[2], ZY[2], Lz/2), Point3f(ZX[1], ZY[2], Lz/2))) do p
        q = Makie.project(ax.scene, p)
        Point2f(Float32(vp.origin[1] + q[1]), Float32(vp.origin[2] + q[2]))
    end
    # Frame the CONTENT, not the cell. The inset axis keeps a tall cell while the corner block
    # it draws is short and wide, so a viewport-sized frame is mostly empty space. Project the
    # eight corners of the zoom volume through the inset camera and fit the frame to them.
    xlo = max(ZX[1], -Lx/2); xhi = min(ZX[2], Lx/2)
    ylo = max(ZY[1], -Ly/2); yhi = min(ZY[2], Ly/2)
    projz = [begin
                 q = Makie.project(axz.scene, Point3f(xx, yy, zz))
                 Point2f(Float32(vz.origin[1] + q[1]), Float32(vz.origin[2] + q[2]))
             end
             for xx in (xlo, xhi), yy in (ylo, yhi), zz in (-Lz/2, Lz/2)]
    pad = 10f0
    bx0 = minimum(p[1] for p in projz) - pad; bx1 = maximum(p[1] for p in projz) + pad
    by0 = minimum(p[2] for p in projz) - pad; by1 = maximum(p[2] for p in projz) + pad
    # square the frame about the content centre: the corner block projects wider than it is tall
    side = max(bx1 - bx0, by1 - by0)
    cx = (bx0 + bx1) / 2;  cy = (by0 + by1) / 2
    zx0 = cx - side/2; zx1 = cx + side/2
    zy0 = cy - side/2; zy1 = cy + side/2
    lines!(fig.scene, [Point2f(zx0,zy0), Point2f(zx1,zy0), Point2f(zx1,zy1),
                       Point2f(zx0,zy1), Point2f(zx0,zy0)];
        color = :black, linewidth = 1.2, space = :pixel)
    # leaders from the window's screen-extreme corners to the near side of the frame
    itop = argmax([p[2] for p in proj]); ibot = argmin([p[2] for p in proj])
    for (src, dst) in ((proj[itop], Point2f(zx0, zy1)), (proj[ibot], Point2f(zx0, zy0)))
        lines!(fig.scene, [src, dst]; color = (:black, 0.5), linewidth = 1.0,
            linestyle = :dash, space = :pixel)
    end

    outpath = joinpath(@__DIR__, "figs/fig8_neighbor_list.pdf")
    save(outpath, fig; px_per_unit = PX_PER_UNIT)
    @info "Saved figure" outpath

    fig
end
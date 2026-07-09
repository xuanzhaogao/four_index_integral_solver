#=
Figure 2 plot script. Loads pre-computed data from fig2_data.jls and
generates fig2_rhs_adaptive.png. Re-run this freely to iterate on the figure;
re-run fig2_data.jl only if you change the experiment.

Layout:
  (a) 3D view of the cube with the adaptive panelization (the three
      camera-facing faces, so the cube reads as a solid opaque body),
      colored by refinement level; source projection marked.
  (b) Relative global interpolation error E_f vs. number of boundary
      unknowns N for RHS-adaptive and uniform refinement.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using GeometryBasics
using Printf

include(joinpath(@__DIR__, "fig_style.jl"))

const datapath = joinpath(@__DIR__, "fig2_data.jls")
const data = open(deserialize, datapath, "r")

const Lx = data.geometry.Lx
const Ly = data.geometry.Ly
const Lz = data.geometry.Lz
const x0 = data.source.x0

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
begin
    fig = Figure(size = (FIG_W, FIG_H))

    # Panel (a): 3D cube panelization colored by refinement level
    panel_records = data.panel_records
    levels = [r.depth for r in panel_records]
    lvl_max = maximum(levels)
    cmap_a = cgrad(FIELD_CMAP, lvl_max + 1, categorical = true)

    # View direction: keep azimuth/elevation in sync with the back-face cull below.
    az_a = 0.35π
    el_a = 0.22π
    ax_a = Axis3(fig[1, 1];
                aspect = :data,
                # title = "(a) adaptive panelization (ε = $(@sprintf("%.0e", data.plot_ε)))",
                xlabel = "x", ylabel = "y", zlabel = "z",
                azimuth = az_a, elevation = el_a,
                protrusions = (40, 20, 25, 25))

    # Eye direction (from cube centre toward the camera). CairoMakie composites
    # the wireframe overlay on top of the whole mesh without depth-testing it, so
    # the back-face grid otherwise bleeds through and the cube reads as
    # transparent. Drawing only the camera-facing faces makes it a solid body.
    eye_dir = (cos(el_a) * cos(az_a), cos(el_a) * sin(az_a), sin(el_a))
    face_visible(n) = (n[1]*eye_dir[1] + n[2]*eye_dir[2] + n[3]*eye_dir[3]) > 1e-9

    # Build a mesh for each panel: 2 triangles, with edge lines on top
    verts_all  = Point3f[]
    faces_all  = TriangleFace{Int}[]
    colors_all = Float32[]
    edge_segs  = Point3f[]

    for rec in panel_records
        face_visible(rec.normal) || continue
        a, b, c, d = rec.corners
        base = length(verts_all)
        push!(verts_all, Point3f(a[1], a[2], a[3]))
        push!(verts_all, Point3f(b[1], b[2], b[3]))
        push!(verts_all, Point3f(c[1], c[2], c[3]))
        push!(verts_all, Point3f(d[1], d[2], d[3]))
        push!(faces_all, TriangleFace{Int}(base+1, base+2, base+3))
        push!(faces_all, TriangleFace{Int}(base+1, base+3, base+4))
        cval = Float32(rec.depth)
        push!(colors_all, cval); push!(colors_all, cval)
        push!(colors_all, cval); push!(colors_all, cval)

        push!(edge_segs, Point3f(a[1], a[2], a[3])); push!(edge_segs, Point3f(b[1], b[2], b[3]))
        push!(edge_segs, Point3f(b[1], b[2], b[3])); push!(edge_segs, Point3f(c[1], c[2], c[3]))
        push!(edge_segs, Point3f(c[1], c[2], c[3])); push!(edge_segs, Point3f(d[1], d[2], d[3]))
        push!(edge_segs, Point3f(d[1], d[2], d[3])); push!(edge_segs, Point3f(a[1], a[2], a[3]))
    end

    mesh!(ax_a, GeometryBasics.Mesh(verts_all, faces_all);
        color = colors_all, colormap = cmap_a,
        colorrange = (-0.5, lvl_max + 0.5),
        shading = NoShading, rasterize = 4)
    linesegments!(ax_a, edge_segs; color = (:black, 0.45), linewidth = 0.3)

    # Source marker (3D) and a thin stem dropped to its projection on +z face
    scatter!(ax_a, [Point3f(x0[1], x0[2], x0[3])];
            color = :red, marker = :xcross, markersize = 14, strokewidth = 2)
    lines!(ax_a, [Point3f(x0[1], x0[2], x0[3]), Point3f(x0[1], x0[2], Lz/2)];
        color = :red, linewidth = 1.0, linestyle = :dash)
    scatter!(ax_a, [Point3f(x0[1], x0[2], Lz/2)];
            color = :red, marker = :circle, markersize = 8)

    Colorbar(fig[1, 2]; colormap = cmap_a,
            limits = (-0.5, lvl_max + 0.5),
            label = "refinement level",
            ticks = 0:lvl_max, width = 12)

    # Panel (b): convergence
    ax_c = Axis(fig[1, 3];
                xscale = log10, yscale = log10,
                xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
                yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5),
                xlabel = L"N",
                ylabel = L"\mathcal{E}_f",
                yticks = ([10.0^i for i in -12:2:0],
                          [rich("10", superscript(string(i))) for i in -12:2:0]),
                # title  = "(c) RHS interpolation convergence"
                )

    p_palette = sweep_colors(length(data.sweeps))
    p_markers = sweep_markers(length(data.sweeps))
    # marker + color both track p (redundant, B&W-safe); adaptive vs uniform is
    # the second category, carried by solid vs dashed line.
    for (i, sw) in enumerate(data.sweeps)
        col = p_palette[i]
        mk  = p_markers[i]
        Ns_a = [r.N for r in sw.adaptive]
        Es_a = [r.Ef for r in sw.adaptive]
        Ns_u = [r.N for r in sw.uniform]
        Es_u = [r.Ef for r in sw.uniform]

        scatterlines!(ax_c, Ns_a, Es_a; color = col, marker = mk,
                      markersize = MS, linewidth = LW_DATA,
                      label = L"\text{adaptive}, p=%$(sw.p)")
        scatterlines!(ax_c, Ns_u, Es_u; color = col, marker = mk,
                      markersize = MS, linewidth = LW_DATA, linestyle = :dot,
                      label = L"\text{uniform}, p=%$(sw.p)")
    end
    axislegend(ax_c; position = :rt, nbanks = 1)
    xlims!(ax_c, 10^(1.8), 10^(8.2))
    ylims!(ax_c, 10^(-12.5), 10^(0.5))


    colgap!(fig.layout, 1, 6)
    colgap!(fig.layout, 2, 22)

    outpath = joinpath(@__DIR__, "figs/fig2_rhs_adaptive.pdf")
    save(outpath, fig; px_per_unit = PX_PER_UNIT)
    @info "Saved figure" outpath

    fig
end
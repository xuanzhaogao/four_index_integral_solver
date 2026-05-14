#=
Figure 2 plot script. Loads pre-computed data from fig2_data.jls and
generates fig2_rhs_adaptive.png. Re-run this freely to iterate on the figure;
re-run fig2_data.jl only if you change the experiment.

Layout:
  (a) 3D view of the cube with the adaptive panelization (all 6 faces),
      colored by refinement level; source projection marked.
  (b) Relative global interpolation error E_f vs. number of boundary
      unknowns N for RHS-adaptive and uniform refinement.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using GeometryBasics
using Printf

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
    fig = Figure(size = (1000, 450), fontsize = 18)

    # Panel (a): 3D cube panelization colored by refinement level
    panel_records = data.panel_records
    levels = [r.depth for r in panel_records]
    lvl_max = maximum(levels)
    cmap_a = cgrad(:viridis, lvl_max + 1, categorical = true)

    ax_a = Axis3(fig[1, 1];
                aspect = :data,
                # title = "(a) adaptive panelization (ε = $(@sprintf("%.0e", data.plot_ε)))",
                xlabel = "x", ylabel = "y", zlabel = "z",
                azimuth = 0.35π, elevation = 0.22π,
                protrusions = (40, 20, 25, 25))

    # Build a mesh for each panel: 2 triangles, with edge lines on top
    verts_all  = Point3f[]
    faces_all  = TriangleFace{Int}[]
    colors_all = Float32[]
    edge_segs  = Point3f[]

    for rec in panel_records
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
        shading = NoShading)
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
                xlabel = "DOF",
                ylabel = L"\mathcal{E}_f",
                # title  = "(c) RHS interpolation convergence"
                )

    Ns_a = [r.N for r in data.adaptive]
    Es_a = [r.Ef for r in data.adaptive]
    Ns_u = [r.N for r in data.uniform]
    Es_u = [r.Ef for r in data.uniform]

    scatterlines!(ax_c, Ns_a, Es_a; color = :crimson, marker = :circle,
                markersize = 12, linewidth = 2, label = "RHS-adaptive")
    scatterlines!(ax_c, Ns_u, Es_u; color = :royalblue, marker = :rect,
                markersize = 12, linewidth = 2, label = "uniform")
    axislegend(ax_c; position = :lb)

    colgap!(fig.layout, 1, 6)
    colgap!(fig.layout, 2, 22)

    outpath = joinpath(@__DIR__, "fig2_rhs_adaptive.png")
    save(outpath, fig, px_per_unit = 4)
    @info "Saved figure" outpath

    fig
end
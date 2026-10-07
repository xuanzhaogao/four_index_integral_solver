# Visualize the 6.1 slab system geometry: 10x10x1 slab (eps 10), Gaussian
# source at the center. RHS-adaptive mesh at sweep parameters (MakieExt
# rendering), with a top-view inset of the source-driven refinement.
# Output: figs/fig61_system_slab_geometry.{png,pdf}

include(joinpath(@__DIR__, "..", "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using CairoMakie, Printf

const SYS = slab_internal()
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

const ADAPTIVE = get(ENV, "VIZ_ADAPTIVE", "1") == "1"
const VIZ_EPS, VIZ_P, VIZ_R = 1e-4, 6, 2
iface = if ADAPTIVE
    vs = gaussian_source(SYS.src_center, SYS.src_sigma, VIZ_EPS)
    svs = Harness.screened_source(SYS, vs)
    kmax = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(svs))
    BI.multi_dielectric_box3d_rhs_adaptive(VIZ_P, Harness.l_ec_of(SYS, VIZ_R),
        SYS.boxes, SYS.epses, svs, 1.0, VIZ_EPS, SYS.eps_out;
        max_depth = 128, tkm_kmax = kmax)
else
    BI.multi_dielectric_box3d(2, Harness.l_ec_of(SYS, 1), SYS.boxes, SYS.epses, SYS.eps_out)
end
panels = iface.panels
@printf("panels: %d   points: %d\n", length(panels), BI.num_points(iface))

fig = Figure(size = (1000, 460))
ax = Axis3(fig[1, 1]; aspect = :data, azimuth = 1.72π, elevation = 0.18π,
    title = "6.1 slab system: 10×10×1, ε = 10, source (σ = 0.05) at center",
    titlealign = :left, xlabel = "x", ylabel = "y", zlabel = "z")

const MExt = Base.get_extension(BI, :MakieExt)
MExt.viz_3d!(ax, iface; show_points = false, highlight_edges = false,
    base_color = "#56B4E9")
scatter!(ax, [Point3f(SYS.src_center...)]; color = :red, markersize = 14, label = "source")
scatter!(ax, [Point3f(SYS.tgt_center...)]; color = :black, marker = :diamond,
    markersize = 10, label = "V target")
axislegend(ax; position = :rt, framevisible = false, labelsize = 11)

# top view of the slab top face (z = 0.5): source-driven refinement
ax2 = Axis(fig[1, 2]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
    title = ADAPTIVE ? "top face (z = 0.5): RHS-adaptive panels" : "top face (z = 0.5)")
seg = Point2f[]
for i in eachindex(panels)
    c = (panels[i].corners[1] .+ panels[i].corners[2] .+ panels[i].corners[3] .+ panels[i].corners[4]) ./ 4
    abs(c[3] - 0.5) < 1e-9 && abs(panels[i].normal[3]) > 0.999 || continue
    cs = panels[i].corners
    for (a, b) in ((1, 2), (2, 3), (3, 4), (4, 1))
        push!(seg, Point2f(cs[a][1], cs[a][2])); push!(seg, Point2f(cs[b][1], cs[b][2]))
    end
end
linesegments!(ax2, seg; color = "#56B4E9", linewidth = 0.7)
scatter!(ax2, [Point2f(SYS.src_center[1], SYS.src_center[2])]; color = :red, markersize = 10)
colsize!(fig.layout, 2, Auto(0.5))

save(joinpath(FIGS, "fig61_system_slab_geometry.png"), fig; px_per_unit = 2)
save(joinpath(FIGS, "fig61_system_slab_geometry.pdf"), fig)
println("wrote figs/fig61_system_slab_geometry.{png,pdf}")

# Visualize the Fig.-1 system geometry: 10x10x1 slab on two 10x10x10 cubes.
# Uses the package's MakieExt (BI.viz_3d!) so panels render exactly as the
# solver sees them; one viz call per interface type to color regions.
# Also prints the shared-face census (geometry sanity, 6.1.7-style).
# Output: figs/fig61_system_fig1_geometry.{png,pdf}

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using CairoMakie, Printf

const SYS = system_fig1()
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

# coarse, non-adaptive build: p = 2, one edge-refinement level
iface = BI.multi_dielectric_box3d(2, Harness.l_ec_of(SYS, 1), SYS.boxes, SYS.epses, SYS.eps_out)
panels = iface.panels
println("panels: ", length(panels), "   points: ", BI.num_points(iface))

pairkey(i) = (round(min(iface.eps_in[i], iface.eps_out[i]); digits = 6),
              round(max(iface.eps_in[i], iface.eps_out[i]); digits = 6))
LABEL = Dict(
    (1.0, 4.0) => "Ω₁ | vac",
    (1.0, 12.0) => "Ω₂ | vac",
    (4.0, 12.0) => "Ω₁ | Ω₂  (shared)",
    (1.0, 10.0) => "slab | vac",
    (4.0, 10.0) => "slab | Ω₁ (contact)",
    (10.0, 12.0) => "slab | Ω₂ (contact)",
)
COLOR = Dict(
    (1.0, 4.0) => "#E8C547",
    (1.0, 12.0) => "#9aa0a6",
    (4.0, 12.0) => "#CC79A7",
    (1.0, 10.0) => "#56B4E9",
    (4.0, 10.0) => "#D55E00",
    (10.0, 12.0) => "#009E73",
)

census = Dict{Tuple{Float64, Float64}, Int}()
for i in eachindex(panels)
    census[pairkey(i)] = get(census, pairkey(i), 0) + 1
end
println("\nface census (eps pair => panel count):")
for (k, v) in sort(collect(census))
    @printf("  %-22s  %5d panels\n", get(LABEL, k, string(k)), v)
end
let seen = Dict{NTuple{3, Float64}, Int}(), dup = 0
    for p in panels
        c = (p.corners[1] .+ p.corners[2] .+ p.corners[3] .+ p.corners[4]) ./ 4
        key = (round(c[1]; digits = 9), round(c[2]; digits = 9), round(c[3]; digits = 9))
        dup += (seen[key] = get(seen, key, 0) + 1) > 1 ? 1 : 0
    end
    println("duplicate-center panels: ", dup, " (expected 0)")
end

# split the interface by eps-pair and render each group via the package viz
function sub_interface(key)
    ids = [i for i in eachindex(panels) if pairkey(i) == key]
    BI.DielectricInterface(panels[ids], iface.eps_in[ids], iface.eps_out[ids])
end

fig = Figure(size = (1000, 620))
ax = Axis3(fig[1, 1]; aspect = :data, azimuth = 1.72π, elevation = 0.16π,
    title = "Fig.-1 system: 10×10×1 slab (ε=10) on two 10×10×10 cubes (ε=4, ε=12)",
    xlabel = "x", ylabel = "y", zlabel = "z")

const MExt = Base.get_extension(BI, :MakieExt)
for (k, _) in sort(collect(census))
    MExt.viz_3d!(ax, sub_interface(k);
        show_points = false, highlight_edges = false, base_color = COLOR[k])
end
scatter!(ax, [Point3f(SYS.src_center...)]; color = :red, markersize = 14)

for (k, v) in sort(collect(census))
    scatter!(ax, [Point3f(NaN, NaN, NaN)]; color = COLOR[k], marker = :rect,
             markersize = 14, label = get(LABEL, k, string(k)) * " ($v)")
end
scatter!(ax, [Point3f(NaN, NaN, NaN)]; color = :red, markersize = 10, label = "source")
axislegend(ax; position = :rt, framevisible = false, labelsize = 11)

save(joinpath(FIGS, "fig61_system_fig1_geometry.png"), fig; px_per_unit = 2)
save(joinpath(FIGS, "fig61_system_fig1_geometry.pdf"), fig)
println("\nwrote figs/fig61_system_fig1_geometry.{png,pdf}")

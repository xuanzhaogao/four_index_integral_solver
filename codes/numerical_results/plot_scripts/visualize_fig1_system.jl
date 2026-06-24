# Visualize the Fig.-1 system geometry: 10x10x1 slab on two 10x10x10 cubes.
# Uses the package's MakieExt (BI.viz_3d!) so panels render exactly as the
# solver sees them; a single uniform-color wireframe, source marked.
# Also prints the shared-face census (geometry sanity, 6.1.7-style).
# Output: figs/fig61_system_fig1_geometry.pdf

include(joinpath(@__DIR__, "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using CairoMakie, Printf
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

const SYS = system_fig1()
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

# RHS-adaptive build at sweep parameters (set VIZ_ADAPTIVE=0 for the coarse
# geometry-only mesh): panels refine near the source according to the screened
# RHS at tolerance eps, exactly as in the benchmark runs.
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
# single uniform wireframe color for the whole interface (no per-category coding)
const GEOM_COLOR = RGBAf(0.40, 0.52, 0.68, 1.0)

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

fig = Figure(size = (FIG_W, FIG_H_3D))
ax = Axis3(fig[1, 1]; aspect = :data, azimuth = 1.72π, elevation = 0.16π,
    xlabel = "x", ylabel = "y", zlabel = "z")

const MExt = Base.get_extension(BI, :MakieExt)
MExt.viz_3d!(ax, iface; show_points = false, highlight_edges = false,
    base_color = GEOM_COLOR)
scatter!(ax, [Point3f(SYS.src_center...)]; color = :red, markersize = 14)

save(joinpath(FIGS, "fig61_system_fig1_geometry.pdf"), fig; px_per_unit = PX_PER_UNIT)
println("\nwrote figs/fig61_system_fig1_geometry.pdf")

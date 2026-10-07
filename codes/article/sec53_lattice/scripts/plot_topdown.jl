# Top-down (x–y) view of a campaign's orbital lattice + slab footprint + junction line,
# to eyeball the supercell structure. Run:
#   julia --project=codes/article scripts/plot_topdown.jl [campaign.toml] [out.pdf]
using BoundaryIntegral, CairoMakie
const BI = BoundaryIntegral

toml = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml")
out = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(@__DIR__, "..", "figs", "fig_benchmark_het3x_topdown.pdf")

c = BI.load_campaign(toml)
orbs = c.orbitals

# slab = thinnest box (smallest Lz); junction = +x face of the lower-x thick cube
slab = c.boxes[argmin([b.Lz for b in c.boxes])]
cubes = [(b, e) for (b, e) in zip(c.boxes, c.epses) if b !== slab]
sort!(cubes, by = t -> t[1].center[1])
xjunc = cubes[1][1].center[1] + cubes[1][1].Lx / 2

fig = Figure(size = (760, 720))
ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x (Å)", ylabel = "y (Å)",
          title = "$(c.name) — top-down (z = $(orbs[1].pos[3]))")

# slab footprint rectangle
sx0, sx1 = slab.center[1] - slab.Lx / 2, slab.center[1] + slab.Lx / 2
sy0, sy1 = slab.center[2] - slab.Ly / 2, slab.center[2] + slab.Ly / 2
lines!(ax, [sx0, sx1, sx1, sx0, sx0], [sy0, sy0, sy1, sy1, sy0];
       color = :green, linestyle = :dash, linewidth = 2, label = "ε=10 slab footprint")
vlines!(ax, [xjunc]; color = :black, linestyle = :dash, linewidth = 2, label = "junction (Si|SiO₂)")

pal = Makie.wong_colors()
for (i, t) in enumerate(sort(unique(o.type for o in orbs)))
    pts = [o.pos for o in orbs if o.type == t]
    scatter!(ax, [p[1] for p in pts], [p[2] for p in pts];
             color = pal[i], markersize = 12, label = "type $t")
end
axislegend(ax; position = :rb, framevisible = true)
save(out, fig; px_per_unit = 2)
println("wrote $out")

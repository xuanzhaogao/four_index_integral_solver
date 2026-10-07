# 3D geometry plot for any campaign TOML, via BoundaryIntegral's MakieExt
# `plot_campaign_geometry` (dielectric cuboids colored by ε + orbitals colored by sublattice).
# Run: julia --project=codes/article scripts/plot_benchmark_geometry.jl [campaign.toml] [out.pdf]
using BoundaryIntegral, CairoMakie

toml = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml")
out = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(@__DIR__, "..", "figs", "fig_benchmark_het3x_geometry.pdf")

fig = plot_campaign_geometry(toml)
save(out, fig; px_per_unit = 2)
println("wrote $out")

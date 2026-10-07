# 3D surface plot of onsite U(x,y) from an onsite_grid campaign. The graphene lattice is
# scattered (hexagonal), so we interpolate (x,y,U) onto a regular grid by inverse-distance
# weighting (no external deps) and render an Axis3 surface.
# Run: julia --project=codes/article scripts/plot_onsite_surface.jl campaigns/onsite_grid.toml
using TOML, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))

toml = abspath(get(ARGS, 1, joinpath(@__DIR__, "..", "campaigns", "onsite_grid.toml")))
cfg = TOML.parsefile(toml)
xs = Float64[]; ys = Float64[]; U = Float64[]
for ln in eachline(joinpath(cfg["root"], "onsite_U.tsv"))
    (isempty(ln) || startswith(ln, "orbital")) && continue
    f = split(ln, '\t'); push!(xs, parse(Float64, f[2])); push!(ys, parse(Float64, f[3])); push!(U, parse(Float64, f[5]))
end
isempty(U) && error("no onsite U in $(cfg["root"])/onsite_U.tsv")

# inverse-distance interpolation onto a regular grid over the data's bounding box
const NX, NY, R2 = 80, 64, 6.0^2
xg = range(minimum(xs), maximum(xs); length = NX)
yg = range(minimum(ys), maximum(ys); length = NY)
Z = fill(NaN, NX, NY)
for i in 1:NX, j in 1:NY
    sw = 0.0; swu = 0.0
    @inbounds for k in eachindex(xs)
        d2 = (xs[k] - xg[i])^2 + (ys[k] - yg[j])^2
        d2 <= R2 || continue
        w = 1.0 / (d2 + 1e-6); sw += w; swu += w * U[k]
    end
    sw > 0 && (Z[i, j] = swu / sw)
end

fig = Figure(size = (FIG_W, FIG_H_3D))
ax = Axis3(fig[1, 1]; xlabel = "x (Å)", ylabel = "y (Å)", zlabel = "U (eV)",
           title = "onsite U(x, y)", azimuth = -0.55π, elevation = 0.22π)
s = surface!(ax, xg, yg, Z; colormap = FIELD_CMAP)
Colorbar(fig[1, 2], s; label = "U (eV)")
save(joinpath(@__DIR__, "..", "figs", "fig_$(cfg["name"])_surface.png"), fig; px_per_unit = 2)
println("wrote figs/fig_$(cfg["name"])_surface.png  (n=$(length(U)) sites)")

# Plot the onsite U produced by scripts/onsite_eval.jl (<root>/onsite_U.tsv).
# Auto-detects 1D line (small y-spread) -> U(x) curve, else 2D -> U(x,y) scatter.
# No BoundaryIntegral dependency (TOML + CairoMakie only).
#
# Run:
#   julia --project=codes/lattice_scale codes/lattice_scale/scripts/plot_onsite.jl campaigns/onsite_line.toml
#   julia --project=codes/lattice_scale codes/lattice_scale/scripts/plot_onsite.jl campaigns/onsite_grid.toml

using TOML, CairoMakie
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

const XJUNC = 5.547

toml_path = abspath(get(ARGS, 1, joinpath(@__DIR__, "..", "campaigns", "onsite_line.toml")))
cfg = TOML.parsefile(toml_path)
tsv = joinpath(cfg["root"], "onsite_U.tsv")

xs = Float64[]; ys = Float64[]; U = Float64[]
for ln in eachline(tsv)
    (isempty(ln) || startswith(ln, "orbital")) && continue
    f = split(ln, '\t')
    push!(xs, parse(Float64, f[2])); push!(ys, parse(Float64, f[3])); push!(U, parse(Float64, f[5]))
end
isempty(U) && error("no onsite U in $(tsv) — run onsite_eval.jl first")

const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

if (maximum(ys) - minimum(ys)) < 5.0          # ---- 1D line ----
    p = sortperm(xs)
    side = [x < XJUNC ? QUAL.blue : QUAL.orange for x in xs[p]]
    fig = Figure(size = (FIG_W, FIG_H))
    ax = Axis(fig[1, 1]; xlabel = "x (Å)", ylabel = "U (eV)",
              title = "onsite U across the junction (y ≈ $(round(ys[1]; digits=1)) Å)")
    lines!(ax, xs[p], U[p]; color = (:gray50, 0.8), linewidth = LW_GUIDE)
    scatter!(ax, xs[p], U[p]; color = side, markersize = MS)
    vlines!(ax, [XJUNC]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    text!(ax, XJUNC, maximum(U); text = "Si | SiO₂ junction", align = (:center, :top),
          offset = (6, -2), fontsize = FS_LEGEND)
    save(joinpath(FIGS, "fig_$(cfg["name"]).pdf"), fig; px_per_unit = PX_PER_UNIT)
    println("wrote figs/fig_$(cfg["name"]).pdf  (n=$(length(U)), U = $(round(minimum(U);digits=3))–$(round(maximum(U);digits=3)) eV)")
else                                           # ---- 2D map ----
    fig = Figure(size = (FIG_W, FIG_H_3D))
    ax = Axis(fig[1, 1]; xlabel = "x (Å)", ylabel = "y (Å)", aspect = DataAspect(),
              title = "onsite U(x, y)")
    sc = scatter!(ax, xs, ys; color = U, colormap = FIELD_CMAP, markersize = 7)
    vlines!(ax, [XJUNC]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    Colorbar(fig[1, 2], sc; label = "U (eV)")
    save(joinpath(FIGS, "fig_$(cfg["name"])_2d.pdf"), fig; px_per_unit = PX_PER_UNIT)
    println("wrote figs/fig_$(cfg["name"])_2d.pdf  (n=$(length(U)), U = $(round(minimum(U);digits=3))–$(round(maximum(U);digits=3)) eV)")
end

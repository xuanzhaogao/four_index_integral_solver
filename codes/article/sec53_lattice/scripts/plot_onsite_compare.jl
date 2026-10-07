# Overlay onsite U(x) from several 1D-line campaigns for comparison.
# Run: julia --project scripts/plot_onsite_compare.jl campaigns/onsite_line.toml campaigns/onsite_line_eps8.4.toml
using TOML, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
const XJUNC = 5.547

function load_line(toml)
    cfg = TOML.parsefile(abspath(toml))
    xs = Float64[]; U = Float64[]
    for ln in eachline(joinpath(cfg["root"], "onsite_U.tsv"))
        (isempty(ln) || startswith(ln, "orbital")) && continue
        f = split(ln, '\t'); push!(xs, parse(Float64, f[2])); push!(U, parse(Float64, f[5]))
    end
    p = sortperm(xs)
    m = match(r"_eps([\d.]+)", cfg["name"])
    return xs[p], U[p], "+x cube ε = " * (m === nothing ? "3.9" : m.captures[1])
end

fig = Figure(size = (FIG_W, FIG_H))
ax = Axis(fig[1, 1]; xlabel = "x (Å)", ylabel = "U (eV)",
          title = "onsite U across the junction (y ≈ 10.7 Å)")
cols = sweep_colors(max(length(ARGS), 2))
for (i, t) in enumerate(ARGS)
    x, U, lab = load_line(t)
    scatterlines!(ax, x, U; color = cols[i], marker = :circle, linewidth = LW_DATA,
                  markersize = MS, label = lab)
end
vlines!(ax, [XJUNC]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
axislegend(ax; position = :lt)
save(joinpath(@__DIR__, "..", "figs", "fig_onsite_line_compare.pdf"), fig; px_per_unit = PX_PER_UNIT)
println("wrote figs/fig_onsite_line_compare.pdf")

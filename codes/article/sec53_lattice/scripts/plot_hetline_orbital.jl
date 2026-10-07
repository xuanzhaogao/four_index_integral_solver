# Plot the ×3 heterojunction ORBITAL onsite-U sweep (real graphene orbitals, eV) from
# <root>/onsite_U.tsv — the clean real-orbital analog of scripts/plot_hetsweep.jl.
# Run: julia --project scripts/plot_hetline_orbital.jl
using CairoMakie, DelimitedFiles
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))

tsv = "/mnt/ceph/users/xgao1/four_index/onsite_hetline_3x/onsite_U.tsv"
data = readdlm(tsv, '\t'; skipstart = 1)
x = Float64.(data[:, 2]); U = Float64.(data[:, 5])
p = sortperm(x); x = x[p]; U = U[p]

fig = Figure(size = (FIG_W, FIG_H))
ax = Axis(fig[1, 1];
    xlabel = "x  (Å)      Si  ←   junction   →  SiO₂",
    ylabel = "onsite U  (eV)",
    title = "×3 heterojunction — real graphene-orbital onsite U across the junction")

hlines!(ax, [1.924]; color = (QUAL.blue, 0.5), linestyle = :dot, linewidth = LW_DATA,
        label = "single Si cube ×4 (1.924 eV)")
hlines!(ax, [2.131]; color = (QUAL.orange, 0.5), linestyle = :dot, linewidth = LW_DATA,
        label = "single SiO₂ cube ×4 (2.131 eV)")
vlines!(ax, [5.547]; color = (:gray, 0.6), linestyle = :dash, linewidth = LW_DATA)  # junction

scatterlines!(ax, x, U; color = QUAL.red, marker = :circle, linewidth = LW_DATA,
              markersize = MS, label = "orbital onsite U")

axislegend(ax; position = :lt, framevisible = true, labelsize = 12)
out = joinpath(@__DIR__, "..", "figs", "fig_hetline_orbital.pdf")
save(out, fig; px_per_unit = PX_PER_UNIT)
println("wrote $out")

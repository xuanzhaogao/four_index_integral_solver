# Plot the ×3 heterojunction point-charge sweep (figs/pointcharge_hetsweep_s3.0.tsv):
# induced self-energy U_proxy vs x, from the Si bulk through the junction into the SiO₂ bulk.
# Run: julia --project scripts/plot_hetsweep.jl
using CairoMakie, DelimitedFiles
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

figs = joinpath(@__DIR__, "..", "figs")
data = readdlm(joinpath(figs, "pointcharge_hetsweep_s3.0.tsv"), '\t'; skipstart = 1)
x = Float64.(data[:, 1]); U = Float64.(data[:, 4])

fig = Figure(size = (FIG_W, FIG_H))
ax = Axis(fig[1, 1];
    xlabel = "x  (Å)      Si  ←   junction   →  SiO₂",
    ylabel = "induced self-energy  U_proxy",
    title = "×3 heterojunction (Si | SiO₂, cubes 270³, slab 270×270×9) — point-charge sweep")

# single-cube converged bulk references
hlines!(ax, [0.17]; color = (QUAL.blue, 0.5), linestyle = :dot, linewidth = LW_DATA,
        label = "Si bulk (single cube ~0.17)")
hlines!(ax, [0.37]; color = (QUAL.orange, 0.5), linestyle = :dot, linewidth = LW_DATA,
        label = "SiO₂ bulk (single cube ~0.37)")
vlines!(ax, [0.0]; color = (:gray, 0.6), linestyle = :dash, linewidth = LW_DATA)  # junction

scatterlines!(ax, x, U; color = QUAL.red, marker = :circle, linewidth = LW_DATA,
              markersize = MS, label = "sweep (induced U_proxy)")

# flag the two slab-edge points (not bulk)
scatter!(ax, [x[1], x[end]], [U[1], U[end]]; color = :black, marker = :xcross,
         markersize = MS + 3, label = "slab-edge (contaminated)")

axislegend(ax; position = :lt, framevisible = true, labelsize = 12)
out = joinpath(figs, "fig_hetsweep.pdf")
save(out, fig; px_per_unit = PX_PER_UNIT)
println("wrote $out")

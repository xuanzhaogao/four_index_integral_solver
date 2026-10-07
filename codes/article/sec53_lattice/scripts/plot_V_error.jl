# Two-panel figure: (left) magnitude of the four-index tensor |V[ρ_ij,ρ_kl]| at l_ec=1.14 (the
# finer, reference run); (right) the l_ec discretization error |V(l_ec=2.27) − V(l_ec=1.14)|.
# Both log scale (eV), pairs ordered lexicographically by (i,j). Run:
#   julia --project=codes/article scripts/plot_V_error.jl
using Serialization, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
const CEPH = "/mnt/ceph/users/xgao1/four_index"

d3 = open(deserialize, joinpath(CEPH, "lattice_conv_l3", "V_full_eV.jls"))   # l_ec = 1.14 (reference)
d2 = open(deserialize, joinpath(CEPH, "lattice_conv_l2", "V_full_eV.jls"))   # l_ec = 2.27
perm3 = sortperm(d3.pair_ids); perm2 = sortperm(d2.pair_ids)
@assert d3.pair_ids[perm3] == d2.pair_ids[perm2] "pair sets differ between runs"
V3 = d3.V[perm3, perm3]; V2 = d2.V[perm2, perm2]
n = size(V3, 1)

Aref = abs.(V3)
Aerr = abs.(V2 .- V3)
reflo = 1e-4                                  # display floor for |V|
maskfloor = 1e-3                              # below this |V|, a per-element ratio is denom noise
Rel = Aerr ./ Aref
Rel[Aref .< maskfloor] .= NaN                 # mask insignificant elements
Mref = log10.(clamp.(Aref, reflo, Inf))
Mrel = log10.(Rel)

fig = Figure(size = (1180, 540))
ax1 = Axis(fig[1, 1]; xlabel = "pair (k,l)", ylabel = "pair (i,j)", aspect = DataAspect(),
           yreversed = true, title = "|V[ρ_ij, ρ_kl]|  (eV),  l_ec = 1.14")
hm1 = heatmap!(ax1, 1:n, 1:n, Mref; colormap = :viridis, colorrange = (log10(reflo), maximum(Mref)), rasterize = 4)
Colorbar(fig[1, 2], hm1, label = "log₁₀ |V|  (eV)")

ax2 = Axis(fig[1, 3]; xlabel = "pair (k,l)", ylabel = "pair (i,j)", aspect = DataAspect(),
           yreversed = true, title = "relative l_ec error  |V(2.27) − V(1.14)| / |V(1.14)|")
hm2 = heatmap!(ax2, 1:n, 1:n, Mrel; colormap = :viridis, colorrange = (-4, -1),
               nan_color = (:gray, 0.15), rasterize = 4)
Colorbar(fig[1, 4], hm2, label = "log₁₀ relative error", ticks = ([-4, -3, -2, -1], ["0.01%", "0.1%", "1%", "10%"]))

out = joinpath(@__DIR__, "..", "figs", "fig_benchmark_het3x_V_and_error.pdf")
save(out, fig; px_per_unit = 2)
maxrel = maximum(filter(!isnan, Rel))
println("wrote $out   |V|max=$(round(maximum(Aref);digits=3)) eV  max rel err (|V|>$(maskfloor))=$(round(maxrel*100;digits=2))%")

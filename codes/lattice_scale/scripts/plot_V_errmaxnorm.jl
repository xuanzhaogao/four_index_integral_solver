# Heatmap of the l_ec discretization error normalized by the global max element:
#   E[ij,kl] = |V(l_ec=2.27) − V(l_ec=1.14)| / max|V(l_ec=1.14)|
# (same normalization as the exchange-symmetry metric — well-defined for every element, no mask).
# Run: julia --project=codes/lattice_scale scripts/plot_V_errmaxnorm.jl
using Serialization, CairoMakie
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))
const CEPH = "/mnt/ceph/users/xgao1/four_index"

d3 = open(deserialize, joinpath(CEPH, "lattice_conv_l3", "V_full_eV.jls"))    # l_ec = 1.14 (ref)
d2 = open(deserialize, joinpath(CEPH, "lattice_conv_l2", "V_full_eV.jls"))    # l_ec = 2.27
p3 = sortperm(d3.pair_ids); p2 = sortperm(d2.pair_ids)
@assert d3.pair_ids[p3] == d2.pair_ids[p2] "pair sets differ"
V3 = d3.V[p3, p3]; V2 = d2.V[p2, p2]; n = size(V3, 1)

maxV = maximum(abs.(V3))
E = abs.(V2 .- V3) ./ maxV
flo = 1e-6
M = log10.(clamp.(E, flo, Inf))

fig = Figure(size = (660, 560))
ax = Axis(fig[1, 1]; xlabel = "pair (k,l)", ylabel = "pair (i,j)", aspect = DataAspect(),
          yreversed = true, title = "l_ec error / max|V|:  |V(2.27) − V(1.14)| / max|V|")
hm = heatmap!(ax, 1:n, 1:n, M; colormap = :viridis, colorrange = (log10(flo), maximum(M)), rasterize = 4)
Colorbar(fig[1, 2], hm; label = "|ΔV| / max|V|",
         ticks = ([-5, -4, -3, -2], ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²"]))
out = joinpath(@__DIR__, "..", "figs", "fig_V_errmaxnorm.pdf")
save(out, fig; px_per_unit = 2)
println("wrote $out   max|V|=$(round(maxV;digits=3)) eV,  max E=$(round(maximum(E);sigdigits=3)) (= $(round(100*maximum(E);digits=2))% of max|V|)")

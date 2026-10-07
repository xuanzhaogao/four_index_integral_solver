# Magnitude heatmap of the assembled four-index Coulomb tensor V[ρ_ij, ρ_kl] (eV, log scale)
# for the lattice_10x10_het3x benchmark, from V_full_eV.jls — shows the diagonal-dominant
# structure. Pairs ordered lexicographically by (i,j) so each anchor's neighbour block sits
# together near the diagonal.
# Run: julia --project=codes/article scripts/plot_Vheatmap_het3x.jl [campaign.toml]
using BoundaryIntegral, Serialization, CairoMakie
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))

toml = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml")
c = BI.load_campaign(toml)
d = open(deserialize, joinpath(c.root, "V_full_eV.jls"))
pair_ids, V = d.pair_ids, d.V
n = length(pair_ids)

perm = sortperm(pair_ids)                 # group by anchor i then neighbour j
A = abs.(V[perm, perm])
floor_eV = 1e-4
M = log10.(clamp.(A, floor_eV, Inf))
lo, hi = log10(floor_eV), maximum(M)

fig = Figure(size = (620, 560))
ax = Axis(fig[1, 1]; xlabel = "pair (k,l)", ylabel = "pair (i,j)", aspect = DataAspect(),
          title = "$(c.name):  |V[ρ_ij, ρ_kl]|  (eV),  n=$(n)", yreversed = true)
hm = heatmap!(ax, 1:n, 1:n, M; colormap = :viridis, colorrange = (lo, hi), rasterize = 4)
Colorbar(fig[1, 2], hm, label = "log₁₀ |V|  (eV)")
out = joinpath(@__DIR__, "..", "figs", "fig_benchmark_het3x_Vheatmap.pdf")
save(out, fig; px_per_unit = PX_PER_UNIT)
println("wrote $out   (|V| max $(round(maximum(A); digits = 3)) eV)")

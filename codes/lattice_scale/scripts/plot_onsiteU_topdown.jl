# On-site Hubbard U as a 2D top-down map: orbitals as filled markers coloured by U (x–y plane),
# junction marked. Real per-orbital values (no interpolation), l_ec=1.14 (level 3). Run:
#   julia --project=codes/lattice_scale scripts/plot_onsiteU_topdown.jl
using BoundaryIntegral, Serialization, CairoMakie
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))
const CEPH = "/mnt/ceph/users/xgao1/four_index"

c = BI.load_campaign(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml"))
d = open(deserialize, joinpath(CEPH, "lattice_conv_l3", "V_full_eV.jls"))
pair_ids, V = d.pair_ids, d.V
XJ = 5.547
px = [o.pos[1] - XJ for o in c.orbitals]; py = [o.pos[2] for o in c.orbitals]; norb = length(c.orbitals)
U = fill(NaN, norb); for (a,p) in enumerate(pair_ids); p[1]==p[2] && (U[p[1]] = V[a,a]); end

fig = Figure(size = (640, 560))
ax = Axis(fig[1,1]; xlabel = "Δx to junction (Å)", ylabel = "y (Å)", aspect = DataAspect(),
          title = "onsite Hubbard U (eV),  l_ec = 1.14")
vlines!(ax, [0.0]; color = :black, linestyle = :dash, linewidth = 2)
sc = scatter!(ax, px, py; color = U, colormap = :viridis, markersize = 17,
              colorrange = (minimum(U), maximum(U)))
text!(ax, minimum(px)+0.5, maximum(py)+0.6; text = "Si", fontsize = 14, color = :gray25)
text!(ax, 1.5, maximum(py)+0.6; text = "SiO₂", fontsize = 14, color = :gray25)
ylims!(ax, minimum(py)-1.5, maximum(py)+2.5)
Colorbar(fig[1,2], sc, label = "U (eV)")
out = joinpath(@__DIR__, "..", "figs", "fig_onsiteU_topdown.pdf")
save(out, fig; px_per_unit = 2)
println("wrote $out   U∈[$(round(minimum(U);digits=3)),$(round(maximum(U);digits=3))] eV")

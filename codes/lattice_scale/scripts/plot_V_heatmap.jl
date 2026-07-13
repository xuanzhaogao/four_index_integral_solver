# Two-panel view of the 10x10 four-index Coulomb tensor V_ijkl (eV), multicube campaign.
# Left:  log10|V| over the full 2611x2611 tensor (pairs in campaign order) — diagonal-dominant.
# Right: the onsite element V_iiii (Hubbard U_i) as a SURFACE over the orbital lattice (x,y),
#        IDW-interpolated from the 200 orbital sites; shows the heterojunction + edge variation.
#   julia --project=codes/lattice_scale codes/lattice_scale/scripts/plot_V_heatmap.jl
using CairoMakie, Serialization, TOML, Statistics

root = "/mnt/ceph/users/xgao1/four_index/lattice_10x10_multicube"
d = deserialize(joinpath(root, "V_full_eV.jls"))
V, pid = d.V, d.pair_ids
n = size(V, 1)

# orbital positions (id = order in the campaign TOML), and V_iiii per orbital
orbs = TOML.parsefile(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_multicube.toml"))["orbital"]
px = [Float64(o["x"]) for o in orbs]; py = [Float64(o["y"]) for o in orbs]
norb = length(orbs)
U = fill(NaN, norb)
for (a, (i, j)) in enumerate(pid); i == j && (U[i] = V[a, a]); end
@info "V_iiii" norb min=minimum(U) max=maximum(U) mean=round(mean(U); digits=3)

# inverse-distance-weighted interpolation of U onto a regular grid over the lattice footprint
function idw(gx, gy; p = 4.0)
    num = den = 0.0
    @inbounds for k in 1:norb
        d2 = (gx - px[k])^2 + (gy - py[k])^2
        d2 < 1e-9 && return U[k]
        w = 1 / d2^(p / 2)
        num += w * U[k]; den += w
    end
    return num / den
end
xr = range(minimum(px), maximum(px); length = 90)
yr = range(minimum(py), maximum(py); length = 70)
Ug = [idw(gx, gy) for gx in xr, gy in yr]
xj = 5.547                                            # Si | SiO2 junction (x = lattice centre)

fig = Figure(size = (1180, 500), fontsize = 16)

ax1 = Axis(fig[1, 1]; xlabel = "pair index (k,l)", ylabel = "pair index (i,j)",
    title = "log₁₀|V_ij,kl|  (eV),  $(n)×$(n)", aspect = 1, yreversed = true)
hm1 = heatmap!(ax1, 1:n, 1:n, log10.(abs.(V) .+ 1e-12); colormap = :viridis, colorrange = (-4, 0.6), rasterize = 4)
Colorbar(fig[1, 2], hm1, label = "log₁₀|V| (eV)")

ax3 = Axis3(fig[1, 3]; xlabel = "x (Å)", ylabel = "y (Å)", zlabel = "V_iiii (eV)",
    title = "onsite U_i over the lattice", azimuth = 1.05π, elevation = 0.22π)
sf = surface!(ax3, xr, yr, Ug; colormap = :viridis, rasterize = 4)
scatter!(ax3, px, py, U; color = :black, markersize = 5, rasterize = 4)   # actual orbital sites
Colorbar(fig[1, 4], sf, label = "V_iiii (eV)")

FIGS = joinpath(@__DIR__, "..", "figs")
save(joinpath(FIGS, "fig_V_ijkl_heatmap.png"), fig)
save(joinpath(FIGS, "fig_V_ijkl_heatmap.pdf"), fig)
println("saved figs/fig_V_ijkl_heatmap.{png,pdf}  V_iiii: min=$(round(minimum(U);digits=3)) max=$(round(maximum(U);digits=3)) eV")

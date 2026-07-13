# On-site Hubbard U as a SURFACE over the lattice (x,y), at l_ec=1.14 (level 3).
# U_i (diagonal of V) is IDW-interpolated onto a regular grid over the orbital footprint; the
# actual orbital sites are overlaid as points. Run:
#   julia --project=codes/lattice_scale scripts/plot_onsiteU_surface.jl
using BoundaryIntegral, Serialization, CairoMakie
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))
const CEPH = "/mnt/ceph/users/xgao1/four_index"

c = BI.load_campaign(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml"))
d = open(deserialize, joinpath(CEPH, "lattice_conv_l3", "V_full_eV.jls"))
pair_ids, V = d.pair_ids, d.V
px = [o.pos[1] for o in c.orbitals]; py = [o.pos[2] for o in c.orbitals]; norb = length(c.orbitals)
U = fill(NaN, norb); for (a,p) in enumerate(pair_ids); p[1]==p[2] && (U[p[1]] = V[a,a]); end

function idw(gx, gy; p = 4.0)                      # inverse-distance-weighted interpolation
    num = den = 0.0
    @inbounds for k in 1:norb
        d2 = (gx-px[k])^2 + (gy-py[k])^2
        d2 < 1e-9 && return U[k]
        w = 1/d2^(p/2); num += w*U[k]; den += w
    end
    num/den
end
xr = range(minimum(px), maximum(px); length = 100)
yr = range(minimum(py), maximum(py); length = 90)
Ug = [idw(gx, gy) for gx in xr, gy in yr]

fig = Figure(size = (780, 560))
ax = Axis3(fig[1,1]; xlabel = "x (Å)", ylabel = "y (Å)", zlabel = "U (eV)",
           title = "onsite Hubbard U  (l_ec = 1.14)", azimuth = 0.72π, elevation = 0.24π)
sf = surface!(ax, xr, yr, Ug; colormap = :viridis)
scatter!(ax, px, py, U; color = :black, markersize = 4)      # actual orbital sites
Colorbar(fig[1,2], sf, label = "U (eV)")
out = joinpath(@__DIR__, "..", "figs", "fig_onsiteU_surface.pdf")
save(out, fig; px_per_unit = 2)
println("wrote $out   U∈[$(round(minimum(U);digits=3)),$(round(maximum(U);digits=3))] eV")

# Combined §6.3 figure (3 panels):
#   (a) system geometry, zoomed to the orbital region — graphene honeycomb (A/B sublattices) inside
#       the ε=10 host slab, with the Si|SiO₂ junction plane marked;
#   (b) heatmap of the four-index Coulomb tensor |V[ρ_ij,ρ_kl]| (eV, log), from the l_ec=1.14 run;
#   (c) on-site Hubbard U_i as a 3D column chart over the lattice (x,y).
# Run: julia --project=codes/lattice_scale scripts/plot_fig_lattice.jl
using BoundaryIntegral, Serialization, CairoMakie
const Point3f = CairoMakie.Point3f
const Vec3f = CairoMakie.Vec3f
const Rect3f = CairoMakie.Rect3f
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))
const CEPH = "/mnt/ceph/users/xgao1/four_index"

c = BI.load_campaign(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml"))
d = open(deserialize, joinpath(CEPH, "lattice_conv_l3", "V_full_eV.jls"))   # l_ec = 1.14 (converged ref)
pair_ids, V = d.pair_ids, d.V
n = length(pair_ids)

xs = [o.pos[1] for o in c.orbitals]; ys = [o.pos[2] for o in c.orbitals]; ts = [o.type for o in c.orbitals]
norb = length(c.orbitals)
U = fill(NaN, norb)
for (a, p) in enumerate(pair_ids); p[1] == p[2] && (U[p[1]] = V[a, a]); end
XJ = 5.547

slab = c.boxes[argmin([b.Lz for b in c.boxes])]        # tight ε=10 slab
sx = (slab.center[1] - slab.Lx/2, slab.center[1] + slab.Lx/2)
sy = (slab.center[2] - slab.Ly/2, slab.center[2] + slab.Ly/2)
sz = (slab.center[3] - slab.Lz/2, slab.center[3] + slab.Lz/2)
function box_edges(x0, x1, y0, y1, z0, z1)
    c8 = [(x,y,z) for x in (x0,x1), y in (y0,y1), z in (z0,z1)]
    segs = NTuple{2,NTuple{3,Float64}}[]
    idx(a,b,cc) = c8[a,b,cc]
    for b in 1:2, cc in 1:2; push!(segs,(idx(1,b,cc),idx(2,b,cc))); end
    for a in 1:2, cc in 1:2; push!(segs,(idx(a,1,cc),idx(a,2,cc))); end
    for a in 1:2, b in 1:2; push!(segs,(idx(a,b,1),idx(a,b,2))); end
    segs
end

fig = Figure(size = (1560, 470))

# ---- (a) 3D graphene sheet (C atoms + C–C bonds) in the ε=10 host slab ----
xr = (minimum(xs)-3, maximum(xs)+3); yr = (minimum(ys)-3, maximum(ys)+3); zc = 7.5
axa = Axis3(fig[1,1]; aspect = :data, title = "(a) graphene sheet",
            xlabel = "x", ylabel = "y", zlabel = "z", elevation = 0.30π, azimuth = 0.62π,
            zticks = ([3.0, 12.0], ["3", "12"]))
# ε=10 host slab (wireframe) + buried Si|SiO₂ junction plane
for (p,q) in box_edges(sx..., sy..., sz...)
    lines!(axa, [p[1],q[1]], [p[2],q[2]], [p[3],q[3]]; color = (:seagreen, 0.5), linewidth = 1.0)
end
jverts = [Point3f(XJ, sy[1], sz[1]), Point3f(XJ, sy[2], sz[1]), Point3f(XJ, sy[2], sz[2]), Point3f(XJ, sy[1], sz[2])]
mesh!(axa, jverts, [1 2 3; 1 3 4]; color = (:gray, 0.18), transparency = true)
# C–C bonds (nearest neighbours < 1.6 Å) as one linesegments! call
bx = Float64[]; by = Float64[]; bz = Float64[]
for i in 1:norb, j in i+1:norb
    hypot(xs[i]-xs[j], ys[i]-ys[j]) < 1.6 || continue
    append!(bx, (xs[i], xs[j])); append!(by, (ys[i], ys[j])); append!(bz, (zc, zc))
end
linesegments!(axa, bx, by, bz; color = (:gray30, 0.85), linewidth = 1.4)
# C atoms coloured by sublattice
for (t,col,lab) in ((1,QUAL.blue,"A"),(2,QUAL.orange,"B"))
    m = ts .== t
    scatter!(axa, xs[m], ys[m], fill(zc, count(m)); color = col, markersize = 9, label = "C ($lab)")
end
limits!(axa, xr, yr, (sz[1]-1, sz[2]+1))
axislegend(axa; position = :rt, labelsize = 10, framevisible = true)

# ---- (b) V heatmap ----
perm = sortperm(pair_ids); A = abs.(V[perm,perm]); flo = 1e-4
axb = Axis(fig[1,2]; xlabel = "pair (k,l)", ylabel = "pair (i,j)", aspect = DataAspect(),
           yreversed = true, title = "(b) |V[ρ_ij, ρ_kl]|  (eV)")
hmb = heatmap!(axb, 1:n, 1:n, log10.(clamp.(A, flo, Inf)); colormap = :viridis,
               colorrange = (log10(flo), maximum(log10.(clamp.(A,flo,Inf)))), rasterize = 4)
Colorbar(fig[1,3], hmb, label = "log₁₀ |V| (eV)")

# ---- (c) onsite U 3D columns ----
base = 1.9; dx = 1.0
h = U .- base
axc = Axis3(fig[1,4]; aspect = (1,1,0.7), title = "(c) onsite Hubbard U (eV)",
            xlabel = "x (Å)", ylabel = "y (Å)", zlabel = "U (eV)", elevation = 0.28π, azimuth = 0.72π)
mc = meshscatter!(axc, [Point3f(xs[i], ys[i], base + h[i]/2) for i in 1:norb];
                  marker = Rect3f(Vec3f(-0.5), Vec3f(1)),
                  markersize = [Vec3f(dx, dx, h[i]) for i in 1:norb],
                  color = U, colormap = :viridis)
limits!(axc, xr, yr, (base, maximum(U) + 0.02))
Colorbar(fig[1,5], mc, label = "U (eV)")

colsize!(fig.layout, 2, Relative(0.24))
out = joinpath(@__DIR__, "..", "figs", "fig_lattice_63.pdf")
save(out, fig; px_per_unit = 2)
println("wrote $out   U∈[$(round(minimum(U);digits=3)), $(round(maximum(U);digits=3))] eV")

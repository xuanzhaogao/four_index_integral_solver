# Standalone system-geometry figure: an OVERVIEW of the ×3 Si|SiO₂ heterojunction (both cubes +
# ε=10 slab, wireframe coloured by ε, no axes) with a ZOOM INSET showing the graphene sheet
# (C atoms + C–C bonds). A highlight box in the overview marks the zoomed region.
# Run: julia --project=codes/article scripts/plot_system_fig.jl
using BoundaryIntegral, CairoMakie
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))

c = BI.load_campaign(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml"))
xs = [o.pos[1] for o in c.orbitals]; ys = [o.pos[2] for o in c.orbitals]; ts = [o.type for o in c.orbitals]
norb = length(c.orbitals); zc = 7.5
xr = (minimum(xs), maximum(xs)); yr = (minimum(ys), maximum(ys))
slab = c.boxes[argmin([b.Lz for b in c.boxes])]

function box_segs(cx, cy, cz, Lx, Ly, Lz)
    hx,hy,hz = Lx/2, Ly/2, Lz/2
    C = [(cx+sx*hx, cy+sy*hy, cz+sz*hz) for sx in (-1,1), sy in (-1,1), sz in (-1,1)]
    s = NTuple{2,NTuple{3,Float64}}[]
    for b in 1:2, cc in 1:2; push!(s,(C[1,b,cc],C[2,b,cc])); end
    for a in 1:2, cc in 1:2; push!(s,(C[a,1,cc],C[a,2,cc])); end
    for a in 1:2, b in 1:2; push!(s,(C[a,b,1],C[a,b,2])); end
    s
end
box_segs(b) = box_segs(b.center..., b.Lx, b.Ly, b.Lz)
noaxis!(ax) = (hidedecorations!(ax); ax.xypanelvisible[]=false; ax.yzpanelvisible[]=false;
    ax.xzpanelvisible[]=false;
    try; ax.xspinesvisible[]=false; ax.yspinesvisible[]=false; ax.zspinesvisible[]=false; catch; end)

grad = cgrad(:viridis); εlo, εhi = extrema(c.epses)
εcol(ε) = grad[εhi > εlo ? (ε-εlo)/(εhi-εlo) : 0.5]
bonds() = begin
    bx=Float64[]; by=Float64[]; bz=Float64[]
    for i in 1:norb, j in i+1:norb
        hypot(xs[i]-xs[j], ys[i]-ys[j]) < 1.6 || continue
        append!(bx,(xs[i],xs[j])); append!(by,(ys[i],ys[j])); append!(bz,(zc,zc))
    end
    bx,by,bz
end

fig = Figure(size = (860, 700))

# ---- overview ----
axm = Axis3(fig[1,1]; aspect = :data, elevation = 0.17π, azimuth = 0.58π, protrusions = 0)
noaxis!(axm)
for (b, ε) in zip(c.boxes, c.epses), (p,q) in box_segs(b)
    lines!(axm, [p[1],q[1]], [p[2],q[2]], [p[3],q[3]]; color = εcol(ε), linewidth = 1.4)
end
scatter!(axm, xs, ys, fill(zc, norb); color = :orangered, markersize = 3)
# highlight box marking the zoomed sheet region
hb = box_segs((xr[1]+xr[2])/2, (yr[1]+yr[2])/2, zc, xr[2]-xr[1]+6, yr[2]-yr[1]+6, 12.0)
for (p,q) in hb
    lines!(axm, [p[1],q[1]], [p[2],q[2]], [p[3],q[3]]; color = :red, linewidth = 1.6)
end

# ---- zoom inset: graphene honeycomb ----
axi = Axis3(fig[1,1]; width = Relative(0.46), height = Relative(0.46), halign = :left, valign = :top,
            aspect = :data, elevation = 0.30π, azimuth = 0.60π, protrusions = 0)
noaxis!(axi)
bx,by,bz = bonds()
linesegments!(axi, bx, by, bz; color = (:gray25, 0.9), linewidth = 1.3)
for (t,col) in ((1,QUAL.blue),(2,QUAL.orange))
    m = ts .== t
    scatter!(axi, xs[m], ys[m], fill(zc, count(m)); color = col, markersize = 7)
end
limits!(axi, (xr[1]-2, xr[2]+2), (yr[1]-2, yr[2]+2), (zc-3, zc+3))

Colorbar(fig[1,2], limits = (εlo, εhi), colormap = :viridis, label = "ε", height = Relative(0.5))
out = joinpath(@__DIR__, "..", "figs", "fig_system_standalone.pdf")
save(out, fig; px_per_unit = 2)
println("wrote $out")

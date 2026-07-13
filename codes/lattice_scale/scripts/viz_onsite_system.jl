# Pre-flight visualization of the onsite-U probe system (Phase 2 full-density lattice):
# the Si/SiO2 heterojunction substrate + the margin-clipped graphene lattice of probe
# orbital positions (z = 7.5 slab mid-plane). Two panels:
#   (a) top view (x,y): slab footprint, margin box, junction, lattice colored by side
#   (b) side view (x,z): slab vs cube tops, orbital plane + support band (clearance check)
#
# CairoMakie only (no solver) -> fast. Run:
#   julia --project=codes/lattice_scale codes/lattice_scale/scripts/viz_onsite_system.jl
#
# r_supp here is an ESTIMATE of the orbital (phi^2) support radius; the real value will
# be measured from the Wannier template at support_rtol and used to set the margin.

using CairoMakie
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

# ---- substrate boxes [cx, cy, cz, Lx, Ly, Lz, eps] (lattice_10x10_multicube.toml) ----
const SI   = (cx = -39.453, cy = 10.318, cz = -42.0, Lx = 90.0, Ly = 90.0, Lz = 90.0)
const SIO2 = (cx =  50.547, cy = 10.318, cz = -42.0, Lx = 90.0, Ly = 90.0, Lz = 90.0)
const SLAB = (cx =   5.547, cy = 10.318, cz =   7.5, Lx = 90.0, Ly = 90.0, Lz =  9.0)
const XJUNC = 5.547                      # buried Si|SiO2 junction plane
const ZORB  = 7.5                        # orbital plane (slab mid-plane)

# slab footprint and z-extent
slabx = (SLAB.cx - SLAB.Lx/2, SLAB.cx + SLAB.Lx/2)   # (-39.453, 50.547)
slaby = (SLAB.cy - SLAB.Ly/2, SLAB.cy + SLAB.Ly/2)   # (-34.682, 55.318)
slabz = (SLAB.cz - SLAB.Lz/2, SLAB.cz + SLAB.Lz/2)   # (3.0, 12.0)
ztop  = SI.cz + SI.Lz/2                                # cube tops = 3.0 = slab bottom

# ---- graphene lattice (same registry as lattice_10x10_multicube.toml) ----
const A1 = (2.465, 0.0)
const A2 = (-1.2325, 2.1347526)
const ORIG = ((0.00093991, 0.00092506),       # sublattice A
              (-0.00013665, 1.42319335))       # sublattice B

# measured from graphene_00001.xsf at support_rtol=1e-4 (centroid z=7.5):
# in-plane phi^2 support ~10.0 Å, vertical half-extent ~3.83 Å.
const R_XY    = 10.0      # in-plane phi^2 support radius (Å)
const R_Z     = 3.83      # vertical phi^2 half-extent (Å)
const BUFFER  = 1.0
const MARGIN  = R_XY + BUFFER

clipx = (slabx[1] + MARGIN, slabx[2] - MARGIN)
clipy = (slaby[1] + MARGIN, slaby[2] - MARGIN)

inclip(x, y) = clipx[1] <= x <= clipx[2] && clipy[1] <= y <= clipy[2]

xs_si = Float64[]; ys_si = Float64[]; xs_ox = Float64[]; ys_ox = Float64[]
for Rx in -60:60, Ry in -60:60, o in ORIG
    x = o[1] + Rx*A1[1] + Ry*A2[1]
    y = o[2] + Rx*A1[2] + Ry*A2[2]
    inclip(x, y) || continue
    if x < XJUNC; push!(xs_si, x); push!(ys_si, y)
    else;         push!(xs_ox, x); push!(ys_ox, y); end
end
ntot = length(xs_si) + length(xs_ox)
println("kept $(ntot) probe orbitals  (Si side $(length(xs_si)), SiO2 side $(length(xs_ox)))",
        "  margin=$(MARGIN) Å")

# rectangle helper (corners for poly/lines)
rect(xr, yr) = Point2f[(xr[1],yr[1]), (xr[2],yr[1]), (xr[2],yr[2]), (xr[1],yr[2]), (xr[1],yr[1])]

fig = Figure(size = (FIG_W, FIG_H))

# ----- (a) top view -----
ax1 = Axis(fig[1, 1]; xlabel = "x (Å)", ylabel = "y (Å)", aspect = DataAspect(),
           title = "(a) top view, z = 7.5  ($(ntot) orbitals)")
poly!(ax1, rect(slabx, slaby); color = (QUAL.green, 0.12), strokecolor = QUAL.green, strokewidth = 1.5)
lines!(ax1, rect(clipx, clipy); color = :gray40, linestyle = :dash, linewidth = LW_GUIDE)
vlines!(ax1, [XJUNC]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
scatter!(ax1, xs_si, ys_si; color = QUAL.blue,   markersize = 4, label = "over Si")
scatter!(ax1, xs_ox, ys_ox; color = QUAL.orange, markersize = 4, label = "over SiO₂")
text!(ax1, XJUNC, slaby[2]; text = "junction", align = (:center, :bottom),
      offset = (0, 2), fontsize = FS_LEGEND)
axislegend(ax1; position = :rb)

# ----- (b) side view (x,z), zoomed on the slab region -----
ax2 = Axis(fig[1, 2]; xlabel = "x (Å)", ylabel = "z (Å)",
           title = "(b) side view: orbital support clearance")
# cube tops (drawn down to z=0 for context; cubes extend to z=-87)
poly!(ax2, rect(slabx, (0.0, ztop)); color = (QUAL.blue, 0.0))  # frame only baseline
poly!(ax2, Point2f[(slabx[1],0),(XJUNC,0),(XJUNC,ztop),(slabx[1],ztop)]; color = (QUAL.blue, 0.18))
poly!(ax2, Point2f[(XJUNC,0),(slabx[2],0),(slabx[2],ztop),(XJUNC,ztop)]; color = (QUAL.orange, 0.18))
# slab
poly!(ax2, rect(slabx, slabz); color = (QUAL.green, 0.18), strokecolor = QUAL.green, strokewidth = 1.5)
# orbital support band [z - R_Z, z + R_Z]  (vertical half-extent at support_rtol=1e-4)
band!(ax2, [clipx[1], clipx[2]], fill(ZORB - R_Z, 2), fill(ZORB + R_Z, 2);
      color = (:red, 0.18))
hlines!(ax2, [ZORB]; xmin = 0, xmax = 1, color = :red, linewidth = LW_DATA)
vlines!(ax2, [XJUNC]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
text!(ax2, slabx[1] + 3, slabz[2]; text = "slab top z=12", align = (:left, :bottom), fontsize = FS_LEGEND)
text!(ax2, slabx[1] + 3, slabz[1]; text = "slab bottom / cube tops z=3", align = (:left, :top), fontsize = FS_LEGEND)
text!(ax2, slabx[1] + 3, ZORB; text = "orbital ± z-extent ($(R_Z)); clears by $(round(slabz[2]-(ZORB+R_Z);digits=2))",
      align = (:left, :bottom), offset = (0, 3), color = :red, fontsize = FS_LEGEND)
ylims!(ax2, 0, 14)

save(joinpath(@__DIR__, "..", "figs", "fig_onsite_system.pdf"), fig; px_per_unit = PX_PER_UNIT)
println("wrote figs/fig_onsite_system.pdf")

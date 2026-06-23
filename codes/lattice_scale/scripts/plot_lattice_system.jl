# System figure for the 10x10 multicube campaign: the two-cube Si/SiO2 heterojunction +
# host slab, and the 200 graphene orbital (atom) positions. Geometry read straight from the
# campaign TOML (boxes + orbital sites) — no solve needed.
# Left:  3D view (boxes translucent, view clipped to the atom region; cubes are 90^3, shown
#        partially). Right: top-down (x-y) lattice over the Si|SiO2 junction.
#   julia --project=codes/lattice_scale codes/lattice_scale/scripts/plot_lattice_system.jl
using CairoMakie, TOML

cfg  = TOML.parsefile(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_multicube.toml"))
rows = cfg["dielectrics"]["boxes"]                       # each: [cx,cy,cz,Lx,Ly,Lz,eps]
orbs = cfg["orbital"]
ox = Float64[o["x"] for o in orbs]; oy = Float64[o["y"] for o in orbs]
oz = Float64[o["z"] for o in orbs]; ty = Int[o["type"] for o in orbs]
A, B = ty .== 1, ty .== 2
xj = rows[1][1] + rows[1][4] / 2                         # Si cube max-x = junction plane
boxcol = [:steelblue, :orangered, :seagreen]            # Si cube | SiO2 cube | slab

fig = Figure(size = (1240, 560), fontsize = 15)

# --- panel A: side cross-section (x-z): heterojunction stack + slab + graphene plane ---
ax1 = Axis(fig[1, 1]; xlabel = "x (Å)", ylabel = "z (Å)",
    title = "side view (x–z): Si | SiO₂ substrate + slab", aspect = DataAspect())
xL, xR, zB, zT = -22.0, 32.0, -24.0, 16.0      # view window (cubes extend to z=-87, clipped)
ztop = rows[3][3] - rows[3][6] / 2              # slab bottom = cube tops (z=3)
zsl  = rows[3][3] + rows[3][6] / 2              # slab top (z=12)
poly!(ax1, Point2f[(xL, zB), (xj, zB), (xj, ztop), (xL, ztop)]; color = (:steelblue, 0.45))  # Si cube
poly!(ax1, Point2f[(xj, zB), (xR, zB), (xR, ztop), (xj, ztop)]; color = (:orangered, 0.45))  # SiO2 cube
poly!(ax1, Point2f[(xL, ztop), (xR, ztop), (xR, zsl), (xL, zsl)]; color = (:seagreen, 0.30))  # slab
scatter!(ax1, ox[A], oz[A]; color = :black,   markersize = 5)        # atoms (projected) at z=7.5
scatter!(ax1, ox[B], oz[B]; color = :magenta, markersize = 5)
vlines!(ax1, [xj]; color = :black, linestyle = :dash, linewidth = 2)
text!(ax1, -18, -10; text = "Si\nε=11.9", color = :steelblue, fontsize = 13)
text!(ax1, 18, -10; text = "SiO₂\nε=3.9", color = :orangered, fontsize = 13)
text!(ax1, -18, 13; text = "slab ε=10", color = :seagreen, fontsize = 12)
text!(ax1, 6, 9; text = "graphene", color = :black, fontsize = 11)
xlims!(ax1, xL, xR); ylims!(ax1, zB, zT)

# --- panel B: top-down (x-y) lattice over the junction ---
ax2 = Axis(fig[1, 2]; xlabel = "x (Å)", ylabel = "y (Å)",
    title = "orbital lattice (top view): Si | SiO₂", aspect = DataAspect())
xlo, xhi = minimum(ox) - 2, maximum(ox) + 2
ylo, yhi = minimum(oy) - 2, maximum(oy) + 2
poly!(ax2, Point2f[(xlo, ylo), (xj, ylo), (xj, yhi), (xlo, yhi)]; color = (:steelblue, 0.12))
poly!(ax2, Point2f[(xj, ylo), (xhi, ylo), (xhi, yhi), (xj, yhi)]; color = (:orangered, 0.12))
vlines!(ax2, [xj]; color = :black, linestyle = :dash, linewidth = 2,
    label = "junction x=$(round(xj; digits=2))")
scatter!(ax2, ox[A], oy[A]; color = :black,   markersize = 10, label = "sublattice A")
scatter!(ax2, ox[B], oy[B]; color = :magenta, markersize = 10, label = "sublattice B")
axislegend(ax2; position = :rb, labelsize = 10)

Label(fig[2, :], "boxes: Si ε=11.9 (blue) | SiO₂ ε=3.9 (orange) | slab ε=10 (green); cubes are 90³ Å (3D view clipped to the flake region)";
    fontsize = 12, color = :gray30)

FIGS = joinpath(@__DIR__, "..", "figs")
save(joinpath(FIGS, "fig_system_atoms.png"), fig)
save(joinpath(FIGS, "fig_system_atoms.pdf"), fig)
println("saved figs/fig_system_atoms.{png,pdf}  ($(length(orbs)) atoms, junction x=$(round(xj; digits=3)))")

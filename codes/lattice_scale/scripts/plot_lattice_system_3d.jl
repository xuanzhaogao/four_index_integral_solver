# 3D system figure via BI.jl's Makie extension (viz_3d): the multicube heterojunction
# interface (two 90^3 cubes Si|SiO2 + slab, coarsely panelized, one colour per eps region)
# with the 200 graphene orbital (atom) sites scattered on top. View clipped to the flake region.
#   julia --project=codes/lattice_scale codes/lattice_scale/scripts/plot_lattice_system_3d.jl
using CairoMakie, TOML       # loading a Makie backend activates BoundaryIntegral's MakieExt
import BoundaryIntegral as BI

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const ORB_DG = BI.read_xsf(joinpath(REF_DIR, "graphene_00001.xsf"))[2]
const ORB_CENTROID = ntuple(d -> Float64(BI.density_centroid(ORB_DG)[d]), 3)

cfg = TOML.parsefile(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_multicube.toml"))
rows = cfg["dielectrics"]["boxes"]
boxes = BI.BoxGeom[(center = (r[1], r[2], r[3]), Lx = r[4], Ly = r[5], Lz = r[6]) for r in rows]
epses = Float64[r[7] for r in rows]
eps_out = Float64(cfg["dielectrics"]["eps_out"])
orbs = cfg["orbital"]
ox = Float64[o["x"] for o in orbs]; oy = Float64[o["y"] for o in orbs]
oz = Float64[o["z"] for o in orbs]; ty = Int[o["type"] for o in orbs]
A, B = ty .== 1, ty .== 2

# coarse interface (visualization only): large l_ec -> clean box wireframe, minimal refinement
steps = BI.snap_orbital(ORB_DG, ORB_CENTROID, ORB_CENTROID)
bt = BI.assemble_lattice_batch([ORB_DG], Dict(1 => BI.OrbitalInstance(1, 1, steps)), [(1, 1)]; support_rtol = 1e-2)
env = BI.envelope_volume_source(bt); kmax = BI._estimate_tkm3dc_kmax(env)
l_ec = minimum(bx.Lz for bx in boxes) * 1.01
interface = BI.multi_dielectric_box3d_rhs_adaptive(2, l_ec, boxes, epses, env, 1e-1;
    eps_out = eps_out, max_depth = 2, tkm_kmax = kmax)
@info "viz interface" panels = length(interface.panels)

M = Base.get_extension(BI, :MakieExt)
@assert M !== nothing "MakieExt not loaded — need a Makie backend"
ein = interface.eps_in
sub(mask) = BI.DielectricInterface(interface.panels[mask], interface.eps_in[mask], interface.eps_out[mask])

fig = Figure(size = (1000, 900), fontsize = 15)
ax = Axis3(fig[1, 1]; aspect = :data, azimuth = 1.15π, elevation = 0.16π,
    xlabel = "x (Å)", ylabel = "y (Å)", zlabel = "z (Å)",
    title = "multicube system (BI viz_3d) + graphene orbitals")
M.viz_3d!(ax, sub(ein .== epses[1]); base_color = :steelblue, show_points = false)   # Si cube
M.viz_3d!(ax, sub(ein .== epses[2]); base_color = :orangered, show_points = false)   # SiO2 cube
M.viz_3d!(ax, sub(ein .== epses[3]); base_color = :seagreen,  show_points = false)   # slab
scatter!(ax, ox[A], oy[A], oz[A]; color = :black,   markersize = 8, rasterize = 4)
scatter!(ax, ox[B], oy[B], oz[B]; color = :magenta, markersize = 8, rasterize = 4)
# full system (no clipping) — aspect=:data shows the whole 90^3 substrate, as in exp66
Label(fig[2, :], "Si cube ε=11.9 (blue) | SiO₂ cube ε=3.9 (orange) | slab ε=10 (green); orbitals A=black, B=magenta. Full 90³ Å system.";
    fontsize = 11, color = :gray30)

FIGS = joinpath(@__DIR__, "..", "figs")
save(joinpath(FIGS, "fig_system_3d.png"), fig)
save(joinpath(FIGS, "fig_system_3d.pdf"), fig)
println("saved figs/fig_system_3d.{png,pdf}  ($(length(interface.panels)) panels, 200 atoms)")

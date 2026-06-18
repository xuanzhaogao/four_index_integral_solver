# 6.6 system figure: the multicube geometry (two 90^3 substrate cubes + slab) together with
# the graphene p_z orbital source, rendered with BI.jl's Makie extension (viz_3d). Builds a
# COARSE interface (visualization only) so it is fast and the panel wireframe is legible.
#
# Run:  julia --project=codes/numerical_results exp66_multicube/scripts/plot_system.jl

using CairoMakie                 # loading a Makie backend activates BoundaryIntegral's MakieExt
import BoundaryIntegral as BI

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const ORB_DG = BI.read_xsf(joinpath(REF_DIR, "graphene_00001.xsf"))[2]
const ORB_CENTROID = ntuple(d -> Float64(BI.density_centroid(ORB_DG)[d]), 3)

const L = 90.0
const SLAB_THICK = L / 10
const EPS1, EPS2, EPS_SLAB, EPS_OUT = 11.9, 3.9, 10.0, 1.0

cx, cy, cz = ORB_CENTROID
h = L / 2; tz = SLAB_THICK; czi = cz - tz / 2 - h
boxes = BI.BoxGeom[
    (center = (cx - h, cy, czi), Lx = L, Ly = L, Lz = L),    # Omega_1 (eps1, Si)
    (center = (cx + h, cy, czi), Lx = L, Ly = L, Lz = L),    # Omega_2 (eps2, SiO2)
    (center = (cx, cy, cz),      Lx = L, Ly = L, Lz = tz)]   # slab (eps_slab)
epses = Float64[EPS1, EPS2, EPS_SLAB]

# central onsite orbital density rho = phi^2 (coarser support for a lighter figure)
steps = BI.snap_orbital(ORB_DG, ORB_CENTROID, ORB_CENTROID)
b = BI.assemble_lattice_batch([ORB_DG], Dict(1 => BI.OrbitalInstance(1, 1, steps)), [(1, 1)]; support_rtol = 1e-2)
orbital = BI.batch_volume_sources(b)[1]
env = BI.envelope_volume_source(b)
kmax = BI._estimate_tkm3dc_kmax(env)

# COARSE interface (viz only): loose tolerances -> few panels, clean wireframe.
l_ec = tz / 2^1 * 1.01
interface = BI.multi_dielectric_box3d_rhs_adaptive(2, l_ec, boxes, epses, env, 1e-1;
    eps_out = EPS_OUT, max_depth = 6, tkm_kmax = kmax)
@info "viz: interface panels = $(length(interface.panels)), orbital pts = $(size(orbital.positions, 2))"

fig = BI.viz_3d(; interfaces = interface, sources = orbital,
    base_color = :gray72, highlight_edges = true, edge_color = :steelblue,
    show_points = false, markersize = 5, alpha = 0.95, max_points = 15000,
    size = (950, 820))

ax = fig.content[1]
ax.azimuth[] = 1.15π
ax.elevation[] = 0.16π
ax.xlabel = "x (bohr)"; ax.ylabel = "y (bohr)"; ax.zlabel = "z (bohr)"
ax.title = "exp66: two 90³ cubes (ε 11.9 | 3.9) + slab (ε 10), graphene p_z orbital source"

mkpath(joinpath(@__DIR__, "..", "figs"))
save(joinpath(@__DIR__, "..", "figs", "fig66_system.png"), fig)
save(joinpath(@__DIR__, "..", "figs", "fig66_system.pdf"), fig)
println("saved figs/fig66_system.{png,pdf}")

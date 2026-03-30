using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

# ── Geometry (same as convergence.jl) ──────────────────────────────────────────
const SUBSTRATE_CENTER = (0.0, 0.0, -15.0)
const SUBSTRATE_SIZE   = (30.0, 30.0, 30.0)
const FILM_CENTER      = (0.0, 0.0,   1.4)
const FILM_SIZE        = (20.0, 20.0,  2.8)

const EPSES  = [10.0, 2.4]
const EPS_OUT = 1.0
const EPS_SRC = 2.4

const BOXES = [
    (center = SUBSTRATE_CENTER, Lx = SUBSTRATE_SIZE[1], Ly = SUBSTRATE_SIZE[2], Lz = SUBSTRATE_SIZE[3]),
    (center = FILM_CENTER,      Lx = FILM_SIZE[1],      Ly = FILM_SIZE[2],      Lz = FILM_SIZE[3]),
]

# ── Orbital source ──────────────────────────────────────────────────────────────
# Wannier orbital for graphene (5×5×1 supercell, shifted so graphene is at z = 8.0 Å)
const XSF_PATH = "/mnt/home/xgao1/work/four_index_integral_solver/density_data/graphene_00001_5x5x1_shifted.xsf"
const GRAPHENE_Z_XSF = 8.00015750   # z of C atom in XSF
const GRAPHENE_Z_SIM = 1.4          # z of graphene film in simulation
const Z_SHIFT = GRAPHENE_Z_SIM - GRAPHENE_Z_XSF  # ≈ −6.6 Å

# ── Mesh parameters ─────────────────────────────────────────────────────────────
const N_QUAD      = 4
const EDGE_REFINE = 3
const RHS_ATOL    = 1e-4

println("Reading XSF orbital: $(basename(XSF_PATH))")
_, datagrid = BI.read_xsf(XSF_PATH)

datagrid.values .*= datagrid.values

vs = BI.VolumeSource(datagrid; shift = (0.0, 0.0, Z_SHIFT), tol = 1e-4)
println("  VolumeSource: $(size(vs.positions, 2)) active points")

# Screen density by local permittivity (ρ → ρ/ε_local)
tol_screen = sqrt(eps(Float64)) * maximum(b -> max(b.Lx, b.Ly, b.Lz), BOXES)
rho_screened = BI._screened_volume_density_multibox(vs, BOXES, EPSES, EPS_OUT, BI.SharpScreening(), tol_screen)
vs_screened = BI.VolumeSource(copy(vs.positions), copy(vs.weights), rho_screened)
println("  Screened VolumeSource: $(size(vs_screened.positions, 2)) points")

l_ec = FILM_SIZE[3] / 2.0^EDGE_REFINE * 1.01
println("\nBuilding RHS-adaptive interface")
println("  n_quad=$N_QUAD  edge_refine=$EDGE_REFINE  l_ec=$(round(l_ec; digits=4))  rhs_atol=$RHS_ATOL")

interface = BI.multi_dielectric_box3d_rhs_adaptive(
    N_QUAD, l_ec, BOXES, EPSES, vs_screened, 1.0, RHS_ATOL, EPS_OUT;
    max_depth = 100,
)
println("  n_pts = $(BI.num_points(interface))")

# ── Visualization ───────────────────────────────────────────────────────────
function make_figure(interface, src_z, n_quad, edge_refine, rhs_atol)
    fig = BI.viz_3d(interface;
                    show_points     = false,
                    highlight_edges = true,
                    base_color      = :steelblue,
                    edge_color      = :orange,
                    size = (800, 700))
    ax = contents(fig[1, 1])[1]
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.zlabel = "z"
    ax.title  = "RHS-adaptive mesh  (n_quad=$n_quad, edge_refine=$edge_refine, atol=$rhs_atol)"
    return fig
end

function make_figure_vs(interface, vs, src_z, n_quad, edge_refine, rhs_atol)
    fig = BI.viz_3d(;interfaces = [interface],
                    sources = [vs],
                    show_points     = false,
                    highlight_edges = true,
                    base_color      = :steelblue,
                    edge_color      = :orange,
                    size = (800, 700))
    ax = contents(fig[1, 1])[1]
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.zlabel = "z"
    ax.title  = "RHS-adaptive mesh  (n_quad=$n_quad, edge_refine=$edge_refine, atol=$rhs_atol)"
    return fig
end

mkpath("figs")

fig = make_figure(interface, GRAPHENE_Z_SIM, N_QUAD, EDGE_REFINE, RHS_ATOL)
save("figs/graphene_slab_rhs_mesh.png", fig)

fig_2 = make_figure_vs(interface, vs_screened, GRAPHENE_Z_SIM, N_QUAD, EDGE_REFINE, RHS_ATOL)
save("figs/graphene_slab_rhs_mesh_vs.png", fig_2)
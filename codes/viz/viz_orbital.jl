using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using GLMakie

orbital_file = joinpath(@__DIR__, "../../density_data/graphene_00002_cropped.xsf")
structure, datagrid = BI.read_xsf(orbital_file)

datagrid.values .*= datagrid.values

vs = BoundaryIntegral.VolumeSource(datagrid, shift = (0.0, 0.0, - 7.920155482424242), tol = 1e-4)
vs_trg = BoundaryIntegral.VolumeSource(datagrid, shift = (20.0, 20.0, 0.0), tol = 1e-4)

L = 90.0
Lx = L
Ly = L
Lz = 9.0

l_panel = 1.0
p = 6
eps_out = 1.0
eps_in = 6.0

max_order = 128
l_ec = 10.0 / 2^4 * 1.01

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, vs, 1.0, l_ec, 1e-6, eps_in, eps_out, Float64)

fig = BI.viz_3d(; interfaces = [interface], sources = [vs], show_points = false, highlight_edges = true)

targets = zeros(3, length(vs_trg.density))
for trg in BI.eachpoint(vs_trg)
    targets[1, trg.global_idx] = trg.point[1]
    targets[2, trg.global_idx] = trg.point[2]
    targets[3, trg.global_idx] = trg.point[3]
end

panel_size_limit = minimum(BI._panel_max_length(panel) for panel in interface.panels)
interface_refined, _, _ = BI._refine_interface_for_targets(interface, targets, panel_size_limit; range_factor = 5.0)

fig = BI.viz_3d(; interfaces = [interface_refined], sources = [vs, vs_trg], show_points = false, highlight_edges = true)

save(joinpath(@__DIR__, "figs/graphene_orbital_2.png"), fig)

rhs = BI.Rhs_dielectric_box3d_fmm3d(interface, vs, eps_in, 1e-6)


fig = BoundaryIntegral.viz_3d_interface_solution(
    interface,
    rhs;
    n_sample = 100,
    log_abs = true,
    add_colorbar = true,
)

save(joinpath(@__DIR__, "figs/graphene_orbital_rhs.png"), fig)

lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(interface, 1e-6, 1e-6, 128, include_edges_src = false, include_edges_trg = false)

sigma, status = Krylov.gmres(lhs, rhs, atol=1e-6, verbose = 1)
total_flux = dot(sigma, BI.all_weights(interface))

fig = BoundaryIntegral.viz_3d_interface_solution(
    interface,
    sigma;
    n_sample = 200,
    log_abs = true,
    add_colorbar = true,
)

save(joinpath(@__DIR__, "figs/graphene_orbital_solution.png"), fig)
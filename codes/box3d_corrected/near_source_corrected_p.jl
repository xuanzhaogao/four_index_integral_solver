include("utils.jl")

df = joinpath(@__DIR__, "data/near_source_corrected_p.csv")
CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], dz = Float64[], potential_u = Float64[], potential_cff = Float64[]))

ps = [4, 5, 6, 7]
rs = [0, 1, 2]

eps_in = 4.0
eps_out = 1.0

fmm_tol = 1e-6
up_tol = 1e-6
rhs_tol = 1e-6
max_order = 128

target = (5.5, 5.5, 0.2)

Lx = 20.0
Ly = 20.0
Lz = 0.5

source = PointSource((5.0, 6.0, Lz / 2 + 0.01), 1.0)

for (p, r) in [(a, b) for a in ps, b in rs]

    l_ec = 1.0 / 2^r * 1.01

    tbox_u, sigma_u, total_flux_u, n_val_u, n_iter_u = solve_single_box3d_adaptive_mesh(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, source)
    potential_u = BI.laplace3d_pottrg_near(tbox_u, target, sigma_u, 1e-10)

    tbox_cff, sigma_cff, total_flux_cff, n_val_cff, n_iter_cff = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, rhs_tol, up_tol, max_order, source, include_edges_src = false, include_edges_trg = false)
    potential_cff = BI.laplace3d_pottrg_near(tbox_cff, target, sigma_cff, 1e-10)

    CSV.write(df, DataFrame(p = p, r = r, L = Lx, l_ec = l_ec, dz = source.point[3] - Lz / 2, potential_u = potential_u, potential_cff = potential_cff), append = true)
end
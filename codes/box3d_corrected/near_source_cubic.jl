include("utils.jl")

df = joinpath(@__DIR__, "data/near_source_cubic_L1.csv")
# CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], total_flux_u = Float64[], total_flux_cff = Float64[], total_flux_ctf = Float64[], total_flux_cft = Float64[], total_flux_ctt = Float64[], n_val_u = Int[], n_val_cff = Int[], n_val_ctf = Int[], n_val_cft = Int[], n_val_ctt = Int[], n_iter_u = Int[], n_iter_cff = Int[], n_iter_ctf = Int[], n_iter_cft = Int[], n_iter_ctt = Int[], potential_u = Float64[], potential_cff = Float64[], potential_ctf = Float64[], potential_cft = Float64[], potential_ctt = Float64[]))

CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], potential_u = Float64[], flux_u = Float64[], n_val_u = Int[], n_iter_u = Int[], potential_n = Float64[], flux_n = Float64[], n_val_n = Int[], n_iter_n = Int[]))

begin
    ps = [3, 4, 5, 6]
    rs = 0:2:8

    eps_in = 4.0
    eps_out = 1.0

    fmm_tol = 1e-6
    up_tol = 1e-7
    max_order = 128

    source = PointSource((0.2, 0.3, 0.51), 1.0)
    target = (0.1, 0.2, 0.4)

    for (p, r) in [(a, b) for a in ps, b in rs] ∪ [(7, 8)]
        for L in [1.0]
            Lx = L
            Ly = L
            Lz = L

            l_ec = 1.0 / 2^r * 1.01

            # uncorrected case
            tbox_u, sigma_u, total_flux_u, n_val_u, n_iter_u = solve_single_box3d_adaptive_mesh(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, source)
            potential_u = BI.laplace3d_pottrg_near(tbox_u, target, sigma_u, up_tol)

            tbox_n, sigma_n, total_flux_n, n_val_n, n_iter_n = solve_single_box3d_nonadaptive_mesh(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, source)
            potential_n = BI.laplace3d_pottrg_near(tbox_n, target, sigma_n, up_tol)


            @show p, r, L, l_ec, potential_u, total_flux_u, n_val_u, n_iter_u, potential_n, total_flux_n, n_val_n, n_iter_n
            CSV.write(df, DataFrame(p = p, r = r, L = L, l_ec = l_ec, potential_u = potential_u, flux_u = total_flux_u, n_val_u = n_val_u, n_iter_u = n_iter_u, potential_n = potential_n, flux_n = total_flux_n, n_val_n = n_val_n, n_iter_n = n_iter_n), append = true)
        end
    end
end
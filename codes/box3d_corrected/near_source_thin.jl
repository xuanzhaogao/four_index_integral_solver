include("utils.jl")

df = joinpath(@__DIR__, "data/near_source_far_trg_thin_box.csv")

CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], potential_a = Float64[], flux_a = Float64[], n_val_a = Int[], n_iter_a = Int[], potential_n = Float64[], flux_n = Float64[], n_val_n = Int[], n_iter_n = Int[], potential_ac = Float64[], flux_ac = Float64[], n_val_ac = Int[], n_iter_ac = Int[]))

begin
    ps = [4, 5, 6, 7]
    rs = 0:2:8

    eps_in = 4.0
    eps_out = 1.0

    fmm_tol = 1e-6
    up_tol = 1e-6
    max_order = 128

    source = PointSource((0.2, 0.3, 0.51), 1.0)
    target = (9.0, 10.0, 11.0)

    for (p, r) in [(a, b) for a in ps, b in rs] ∪ [(7, 9)]
        for L in [10.0]
            Lx = L
            Ly = L
            Lz = 1.0

            l_ec = 1.0 / 2^r * 1.01

            # uncorrected case
            tbox_a, sigma_a, total_flux_a, n_val_a, n_iter_a = solve_single_box3d_adaptive_mesh(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, source)
            potential_a = BI.laplace3d_pottrg_near(tbox_a, target, sigma_a, up_tol)

            tbox_n, sigma_n, total_flux_n, n_val_n, n_iter_n = solve_single_box3d_nonadaptive_mesh(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, source)
            potential_n = BI.laplace3d_pottrg_near(tbox_n, target, sigma_n, up_tol)

            tbox_ac, sigma_ac, total_flux_ac, n_val_ac, n_iter_ac = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, up_tol, max_order, source)
            potential_ac = BI.laplace3d_pottrg_near(tbox_ac, target, sigma_ac, up_tol)

            CSV.write(df, DataFrame(p = p, r = r, L = L, l_ec = l_ec, potential_a = potential_a, flux_a = total_flux_a, n_val_a = n_val_a, n_iter_a = n_iter_a, potential_n = potential_n, flux_n = total_flux_n, n_val_n = n_val_n, n_iter_n = n_iter_n, potential_ac = potential_ac, flux_ac = total_flux_ac, n_val_ac = n_val_ac, n_iter_ac = n_iter_ac), append = true)
        end
    end
end
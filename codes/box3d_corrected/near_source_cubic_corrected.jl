include("utils.jl")

df = joinpath(@__DIR__, "data/near_source_cubic_L1_corrected.csv")
# CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], total_flux_u = Float64[], total_flux_cff = Float64[], total_flux_ctf = Float64[], total_flux_cft = Float64[], total_flux_ctt = Float64[], n_val_u = Int[], n_val_cff = Int[], n_val_ctf = Int[], n_val_cft = Int[], n_val_ctt = Int[], n_iter_u = Int[], n_iter_cff = Int[], n_iter_ctf = Int[], n_iter_cft = Int[], n_iter_ctt = Int[], potential_u = Float64[], potential_cff = Float64[], potential_ctf = Float64[], potential_cft = Float64[], potential_ctt = Float64[]))

CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], potential = Float64[], flux = Float64[], n_val = Int[], n_iter = Int[]))

begin
    ps = [3, 4, 5, 6]
    rs = 0:2:8

    eps_in = 4.0
    eps_out = 1.0

    fmm_tol = 1e-6
    up_tol = 1e-7
    rhs_tol = 1e-6
    max_order = 128

    source = PointSource((0.2, 0.3, 0.51), 1.0)
    target = (0.1, 0.2, 0.4)

    for (p, r) in [(a, b) for a in ps, b in rs] ∪ [(7, 8)]
        for L in [1.0]
            Lx = L
            Ly = L
            Lz = L

            l_ec = 1.0 / 2^r * 1.01

            tbox, sigma, total_flux, n_val, n_iter = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, rhs_tol, up_tol, max_order, source)
            potential = BI.laplace3d_pottrg_near(tbox, target, sigma, up_tol)

            @show p, r, L, l_ec, potential, total_flux, n_val, n_iter
            CSV.write(df, DataFrame(p = p, r = r, L = L, l_ec = l_ec, potential = potential, flux = total_flux, n_val = n_val, n_iter = n_iter), append = true)
        end
    end
end
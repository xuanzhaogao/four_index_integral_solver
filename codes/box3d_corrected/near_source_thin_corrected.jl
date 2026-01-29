include("utils.jl")

df = joinpath(@__DIR__, "data/near_source_thin_corrected.csv")

CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], potential = Float64[], flux = Float64[], n_val = Int[], n_iter = Int[]))

begin
    ps = [4, 5, 6, 7]
    rs = 0:2:8

    eps_in = 4.0
    eps_out = 1.0

    fmm_tol = 1e-6
    up_tol = 1e-6
    max_order = 128

    source = PointSource((0.2, 0.3, 0.51), 1.0)
    target = (0.1, 0.2, 0.4)

    for (p, r) in [(a, b) for a in ps, b in rs] ∪ [(7, 9)]
        for L in [10.0]
            Lx = L
            Ly = L
            Lz = 1.0

            l_ec = 1.0 / 2^r * 1.01

            tbox, sigma, total_flux, n_val, n_iter = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, up_tol, max_order, source)
            potential = BI.laplace3d_pottrg_near(tbox, target, sigma, up_tol)

            @show p, r, L, l_ec, potential, total_flux, n_val, n_iter
            CSV.write(df, DataFrame(p = p, r = r, L = L, l_ec = l_ec, potential = potential, flux = total_flux, n_val = n_val, n_iter = n_iter), append = true)
        end
    end
end
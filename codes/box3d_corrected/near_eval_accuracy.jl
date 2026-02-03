include("utils.jl")

df = joinpath(@__DIR__, "data/near_eval_accuracy.csv")
CSV.write(df, DataFrame(trg_z = Float64[], tol = Float64[], pot = Float64[], pot_ref = Float64[], err = Float64[]))


begin
    p = 6
    r = 4

    eps_in = 4.0
    eps_out = 1.0

    fmm_tol = 1e-4
    up_tol = 1e-5
    max_order = 128

    L = 1.0

    Lx = L
    Ly = L
    Lz = 1.0

    l_ec = 1.0 / 2^r * 1.01

    tbox_cff, sigma_cff, total_flux_cff, n_val_cff, n_iter_cff = solve_single_box3d_nonadaptive_mesh(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, PointSource((0.1, 0.2, 0.3), eps_in))
end

begin
    trg_zs = 0.5 .+ [0.001, 0.01, 0.1, 0.2, 0.5, 1.0]

    for trg_z in trg_zs
        trg = (0.1, 0.2, trg_z)

        pot_ref = BI.laplace3d_pottrg_near(tbox_cff, trg, sigma_cff, 1e-10, range_factor = Inf)

        pots = []
        errs = []
        for tol in [1e-2, 1e-4, 1e-6, 1e-8]
            pot_cff = BI.laplace3d_pottrg_near(tbox_cff, trg, sigma_cff, tol, range_factor = 5.0)
            CSV.write(df, DataFrame(trg_z = trg_z - 0.5, tol = tol, pot = pot_cff, pot_ref = pot_ref, err = abs(pot_cff - pot_ref)), append = true)
        end
    end
end
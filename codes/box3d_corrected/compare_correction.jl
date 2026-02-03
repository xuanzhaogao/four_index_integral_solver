include("utils.jl")

df = joinpath(@__DIR__, "data/compare_correction.csv")

CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], dz = [], l_ec = Float64[], phi_nn = [], phi_nf = [], phi_cn = [], phi_cf =[]))

p = 4
rs = 0:1:6
Lx = 10.0
Ly = 10.0
Lz = 0.1

fmm_tol = 1e-6
rhs_tol = 1e-6
up_tol = 1e-6
eval_tol = 1e-8

for dz in [0.1, 1.0, 10.0]
    source = PointSource((5.0, 6.0, L_z / 2 + dz), 1.0)
    for r in rs
        l_ec = 1.0 / 2^r * 1.01
        box_n, sigma_n, tf_n, n_val_n, n_iter_n = solve_single_box3d_adaptive_mesh(Lx, Ly, Lz, p, l_ec, 4.0, 1.0, fmm_tol, rhs_tol, source)
        phi_nn = BI.laplace3d_pottrg_near(box_n, (0.1, 0.2, L_z / 2 + 0.1), sigma_n, eval_tol)
        phi_nf = BI.laplace3d_pottrg_near(box_n, (0.1, 0.2, 1.5), sigma_n, eval_tol)

        box_c, sigma_c, tf_c, n_val_c, n_iter_c = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p, l_ec, 4.0, 1.0, fmm_tol, rhs_tol, up_tol, 128, source)
        phi_cn = BI.laplace3d_pottrg_near(box_c, (0.1, 0.2, L_z / 2 + 0.1), sigma_c, eval_tol)
        phi_cf = BI.laplace3d_pottrg_near(box_c, (0.1, 0.2, 1.5), sigma_c, eval_tol)

        println(p, r, Lx, dz, l_ec, phi_nn, phi_nf, phi_cn, phi_cf)
        CSV.write(df, DataFrame(p = p, r = r, L = Lx, dz = dz, l_ec = l_ec, phi_nn = phi_nn, phi_nf = phi_nf, phi_cn = phi_cn, phi_cf = phi_cf), append = true)

    end

    # for reference: p = 5, r = 7
    box_ref, sigma_ref, tf_ref, n_val_ref, n_iter_ref = solve_single_box3d_adaptive_mesh(Lx, Ly, Lz, 5, 1.0 / 2^7 * 1.01, 4.0, 1.0, fmm_tol, rhs_tol, source)
    phi_ref_n = BI.laplace3d_pottrg_near(box_ref, (0.1, 0.2, L_z / 2 + 0.1), sigma_ref, eval_tol)
    phi_ref_f = BI.laplace3d_pottrg_near(box_ref, (0.1, 0.2, 1.5), sigma_ref, eval_tol)
    box_c_ref, sigma_c_ref, tf_c_ref, n_val_c_ref, n_iter_c_ref = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, 5, 1.0 / 2^7 * 1.01, 4.0, 1.0, fmm_tol, rhs_tol, up_tol, 128, source)
    phi_c_ref_n = BI.laplace3d_pottrg_near(box_c_ref, (0.1, 0.2, L_z / 2 + 0.1), sigma_c_ref, eval_tol)
    phi_c_ref_f = BI.laplace3d_pottrg_near(box_c_ref, (0.1, 0.2, 1.5), sigma_c_ref, eval_tol)

    CSV.write(df, DataFrame(p = 5, r = 7, L = Lx, dz = "ref", l_ec = 1.0 / 2^7 * 1.01, phi_nn = phi_ref_n, phi_nf = phi_ref_f, phi_cn = phi_c_ref_n, phi_cf = phi_c_ref_f), append = true)
end
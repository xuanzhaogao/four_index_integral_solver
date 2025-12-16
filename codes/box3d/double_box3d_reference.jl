include(joinpath(@__DIR__, "single_box3d_utils.jl"))

trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

n_quad = 4
n_edge = 12

dbox, sigma, gi, n_val, n_iter = solve_double_box3d(3.0, 4.0, n_quad, n_quad, n_edge, n_edge, (0.2, 0.3, 0.4))

S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
pot_trg = S_trg * sigma

res_ref = Dict("dbox" => dbox, "sigma" => sigma, "gi" => gi, "pot_trg" => pot_trg)
save(joinpath(@__DIR__, "data/double_box3d_ref.jld2"), res_ref)
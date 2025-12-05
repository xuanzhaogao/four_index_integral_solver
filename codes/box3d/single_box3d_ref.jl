using JLD2

include(joinpath(@__DIR__, "single_box3d_utils.jl"))
trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

# dbox, sigma, gi, n_val, n_iter = solve_single_box3d(4.0, 1, 12, false, 10, 10, (0.2, 0.3, 0.4), rtol = 1e-8)

# S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-8)
# pot_trg = S_trg * sigma

# res_ref = Dict("dbox" => dbox, "sigma" => sigma, "gi" => gi, "pot_trg" => pot_trg)
# save(joinpath(@__DIR__, "data/single_box3d_ref.jld2"), res_ref)


dbox, sigma, gi, n_val, n_iter = solve_single_box3d(4.0, 1, 12, true, 12, 12, (0.2, 0.3, 0.4), rtol = 1e-8)

S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-8)
pot_trg = S_trg * sigma

res_ref = Dict("dbox" => dbox, "sigma" => sigma, "gi" => gi, "pot_trg" => pot_trg)
save(joinpath(@__DIR__, "data/single_box3d_ref_reduce.jld2"), res_ref)
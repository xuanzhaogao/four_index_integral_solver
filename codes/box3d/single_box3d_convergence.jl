include(joinpath(@__DIR__, "single_box3d_utils.jl"))

res_ref = load(joinpath(@__DIR__, "data/single_box3d_ref_non_reduce.jld2"))
pot_ref = res_ref["pot_trg"]
trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

df = joinpath(@__DIR__, "data/single_box3d_convergence.csv")
# CSV.write(df, DataFrame(n_quad = [], n_edge = [], gi = [], n_val = [], n_iter = [], pot_abserr = [], pot_relerr = []))

src = (0.2, 0.3, 0.4)

# n_quads = [2, 4, 6, 8]
# n_edges = 0:1:10

# n_quads = [10:2:24...]
# n_edges = [3, 4]

# n_quads = [2]
# n_edges = 11:1:16

# n_quads = [3]
# n_edges = 0:1:16

n_quads = [26:2:32...]
n_edges = [3]

for n_quad in n_quads
    for n_edge in n_edges
        dbox, sigma, gi, n_val, n_iter = solve_single_box3d(4.0, 1, n_quad, n_quad, n_edge, n_edge, (0.2, 0.3, 0.4))
        S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
        pot_trg = S_trg * sigma
        pot_abserr = maximum(abs.(pot_trg - pot_ref))
        pot_relerr = maximum(abs.((pot_trg - pot_ref) ./ pot_ref))

        println("n_quad = $n_quad, n_edge = $n_edge, gi = $gi, n_val = $n_val, n_iter = $n_iter, pot_abserr = $pot_abserr, pot_relerr = $pot_relerr")

        CSV.write(df, DataFrame(n_quad = [n_quad], n_edge = [n_edge], gi = [gi], n_val = [n_val], n_iter = [n_iter], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr]), append = true)
    end
end
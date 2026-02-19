using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Krylov
using Random
using CSV, DataFrames

Random.seed!(1234)

L = 90.0
Lx = L
Ly = L
Lz = 9.0

orbital_file = joinpath(@__DIR__, "../../../density_data/graphene_00002_cropped.xsf")
structure, datagrid = BI.read_xsf(orbital_file)
datagrid.values .*= datagrid.values
vs = BoundaryIntegral.VolumeSource(datagrid, shift = (0.0, 0.0, - 7.920155482424242))

eps_out = 1.0
eps_in = 6.0
max_order = 128

p = 4
r = 4

l_ec = 9.0 / 2^r * 1.01

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, vs, eps_in, l_ec, 1e-6, eps_in, eps_out, Float64)

rhs = BI.Rhs_dielectric_box3d_fmm3d(interface, vs, eps_in, 1e-6)
lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(interface, 1e-6, 1e-6, 128, include_edges_src = false, include_edges_trg = false)
sigma, status = Krylov.gmres(lhs, rhs, atol=1e-6, verbose = 1)

n_trg = 100
pos_trg = [(5.0 * randn(), 5.0 * randn(), 5.0 * randn()) for _ in 1:n_trg]
pos_trg_mat = zeros(3, n_trg)
for i in 1:n_trg
    for j in 1:3
        pos_trg_mat[j, i] = pos_trg[i][j]
    end
end

@time for (i, pos) in enumerate(pos_trg)
    pot_trg[i] = BI.laplace3d_pottrg_near(interface, pos, sigma, 1e-12, range_factor = Inf)
end


df = joinpath(@__DIR__, "../data/batch_near_eval.csv")
CSV.write(df, DataFrame(tol = Float64[], range = Float64[], l2_abs_err = Float64[], l2_rel_err = Float64[]))

for tol in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
    for range in [4.0, 6.0, 8.0]
        @time begin
            pot_near_batch_op = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, pos_trg_mat, tol, tol, range)
            pot_near_batch = pot_near_batch_op * sigma
        end

        l2_abs_err = norm(pot_trg .- pot_near_batch)
        l2_rel_err = l2_abs_err / norm(pot_trg)

        CSV.write(df, DataFrame(tol = tol, range = range, l2_abs_err = l2_abs_err, l2_rel_err = l2_rel_err), append=true)
        @info "tol = $tol, range = $range, l2_abs_err = $l2_abs_err, l2_rel_err = $l2_rel_err"
    end
end
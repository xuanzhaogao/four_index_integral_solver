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
vs_trg = BoundaryIntegral.VolumeSource(datagrid, shift = (20.0, 20.0, - 7.920155482424242))

eps_out = 1.0
eps_in = 2.4
max_order = 128

p = 4
r = 4

l_ec = 9.0 / 2^r * 1.01

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, vs, eps_in, l_ec, 1e-6, eps_in, eps_out, Float64)

rhs = BI.Rhs_dielectric_box3d_fmm3d(interface, vs, eps_in, 1e-6)
lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(interface, 1e-6, 1e-6, 128, include_edges_src = false, include_edges_trg = false)
sigma, status = Krylov.gmres(lhs, rhs, atol=1e-6, verbose = 1)

targets = zeros(3, length(vs_trg.density))
for trg in BI.eachpoint(vs_trg)
    targets[1, trg.global_idx] = trg.point[1]
    targets[2, trg.global_idx] = trg.point[2]
    targets[3, trg.global_idx] = trg.point[3]
end

ix = 1
iy = 1
pos_evaluated = Int[]
n_x, n_y, n_z = length.(vs.axes)
for iz in 1:n_z
    push!(pos_evaluated, ix + (iy - 1) * n_x + (iz - 1) * n_x * n_y)
end

pot_trg = zeros(length(pos_evaluated))
@time for (i, pos) in enumerate(pos_evaluated)
    pot_trg[i] = BI.laplace3d_pottrg_near(interface, Tuple(targets[:, pos]), sigma, 1e-12, range_factor = Inf)
end


df = joinpath(@__DIR__, "../data/batch_near_eval.csv")
CSV.write(df, DataFrame(tol = Float64[], range = Float64[], l2_abs_err = Float64[], l2_rel_err = Float64[]))

for tol in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
    for range in [4.0, 6.0, 8.0]
        @time begin
            pot_near_batch_op = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, targets, tol, tol, range)
            pot_near_batch = pot_near_batch_op * sigma
        end

        pot_near_batch_compare = pot_near_batch[pos_evaluated]

        l2_abs_err = norm(pot_trg .- pot_near_batch_compare)
        l2_rel_err = l2_abs_err / norm(pot_trg)

        CSV.write(df, DataFrame(tol = tol, range = range, l2_abs_err = l2_abs_err, l2_rel_err = l2_rel_err), append=true)
        @info "tol = $tol, range = $range, l2_abs_err = $l2_abs_err, l2_rel_err = $l2_rel_err"
    end
end
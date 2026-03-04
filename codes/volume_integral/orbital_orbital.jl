using FBCPoisson, BoxDMK, BoundaryIntegral
using FINUFFT
import BoundaryIntegral as BI

graphene_1 = joinpath(@__DIR__, "../../density_data/graphene_00001_5x5x1_shifted.xsf")
graphene_2 = joinpath(@__DIR__, "../../density_data/graphene_00002_5x5x1_shifted.xsf")

structure_1, datagrid_1 = BI.read_xsf(graphene_1)
datagrid_1.values .*= datagrid_1.values
vs_src_1 = BoundaryIntegral.VolumeSource(datagrid_1, shift = (0.0, 0.0, - 7.920155482424242))
vs_src_1_truncated = BoundaryIntegral.VolumeSource(datagrid_1, shift = (0.0, 0.0, - 7.920155482424242), tol = 1e-6)

structure_2, datagrid_2 = BI.read_xsf(graphene_2)
datagrid_2.values .*= datagrid_2.values
vs_src_2 = BoundaryIntegral.VolumeSource(datagrid_2, shift = (0.0, 0.0, - 7.920155482424242))
vs_src_2_truncated = BoundaryIntegral.VolumeSource(datagrid_2, shift = (0.0, 0.0, - 7.920155482424242), tol = 1e-6)

trg_pot_24 = lfbc3d(24, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_2_truncated.positions, 1e-6, 1)
trg_pot_36 = lfbc3d(36, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_2_truncated.positions, 1e-6, 1)
trg_pot_48 = lfbc3d(48, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_2_truncated.positions, 1e-6, 1)
trg_pot_64 = lfbc3d(64, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_2_truncated.positions, 1e-6, 1)
trg_pot_128 = lfbc3d(128, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_2_truncated.positions, 1e-6, 1)
trg_pot_256 = lfbc3d(256, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_2_truncated.positions, 1e-6, 1)

E_24 = sum(trg_pot_24 .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)
E_36 = sum(trg_pot_36 .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)
E_48 = sum(trg_pot_48 .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)
E_64 = sum(trg_pot_64 .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)
E_128 = sum(trg_pot_128 .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)
E_256 = sum(trg_pot_256 .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)

@show E_24
@show E_36
@show E_48
@show E_64
@show E_128
@show E_256

# BoxDMK verification
# datagrid_1.values already contains orbital_1^2 (squared in-place above)
# Build a shifted datagrid whose origin matches the shift passed to VolumeSource
import BoundaryIntegral: _datagrid_trilinear_value, _datagrid_affine

let o = datagrid_1.origin
    global dg1_shifted = merge(datagrid_1, (origin = typeof(o)(o[1], o[2], o[3] - 7.920155482424242),))
end
_, _, Minv_1, _, _, _ = _datagrid_affine(dg1_shifted)

density_1_cb(x, _) = begin
    v = _datagrid_trilinear_value(dg1_shifted, Minv_1, (x[1], x[2], x[3]))
    isnan(v) ? 0.0 : v
end

# boxlen must contain all source positions of orbital 1
boxlen_bdmk = 2.0 * maximum(abs, vs_src_1.positions) * 1.01

prob_bdmk = LaplaceProblem(density = density_1_cb, nd = 1, ndim = 3, boxlen = boxlen_bdmk)
opts_bdmk = BDMKOptions(eps = 1e-4, norder = 12)
tree_bdmk, _ = solve_problem(prob_bdmk; compute = :potential, opts = opts_bdmk)
result_bdmk = evaluate_targets(prob_bdmk, tree_bdmk, vs_src_2_truncated.positions;
                               compute = :potential, eps = 1e-4)

trg_pot_bdmk = vec(result_bdmk.pote)
E_bdmk = sum(trg_pot_bdmk .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)

# BoxDMK uses kernel 1/|r-r'|; lfbc3d uses 1/(4π|r-r'|), so E_bdmk/(4π) ≈ E_256
println("E_256           = $(E_256)")
println("E_bdmk / (4π)   = $(E_bdmk / (4π))")
println("relative diff   = $(abs(E_bdmk / (4π) - E_256) / abs(E_256))")

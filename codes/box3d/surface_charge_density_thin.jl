using JLD2
using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra

res_ref = load(joinpath(@__DIR__, "cache/single_thin_box_convergence_L10_r8_p6_r8.jld2"))
tbox = res_ref["tbox"]
sigma = res_ref["sigma"]

sigma_check = BI.nystrom_interpolation_dielectric_box3d(tbox, tbox, (0.2, 0.3, 0.4), 4.0, sigma, 1e-6)
println("check error = $(LinearAlgebra.norm(sigma_check - sigma, Inf))")

nm = 500
targets_x = range(-5.0, 5.0, length = nm)
targets_y = range(-5.0, 5.0, length = nm)
targets_z = 0.5

n = nm^2

points = Vector{NTuple{3, Float64}}(undef, n)
for i in 1:nm
    for j in 1:nm
        points[(i - 1) * nm + j] = (targets_x[i], targets_y[j], targets_z)
    end
end
norm = (0.0, 0.0, 1.0)
weights = ones(n)
corners = [(-5.0, -5.0, 0.5), (5.0, -5.0, 0.5), (5.0, 5.0, 0.5), (-5.0, 5.0, 0.5)]

panels = BI.Panel(n, points, norm, weights, corners)
interface = BI.Interface(1, [panels])
target_interface = BI.DielectricInterfaces(1, [(interface, 4.0, 1.0)])

sigma_f = BI.nystrom_interpolation_dielectric_box3d(tbox, target_interface, (0.2, 0.3, 0.4), 4.0, sigma, 1e-6)

save(joinpath(@__DIR__, "data/surface_charge_density_thin_L10_a0.25_p6_r8.jld2"), Dict("sigma_f" => sigma_f, "points" => points))
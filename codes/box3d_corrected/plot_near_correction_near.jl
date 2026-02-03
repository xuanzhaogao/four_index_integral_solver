include("utils.jl")

using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

Lx = 20.0
Ly = 20.0
Lz = 0.5

p = 4
r = 2
eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^r * 1.01

# tbox = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, PointSource((5.0, 6.0, Lz / 2 + 0.01), 1.0), 1.0, l_ec, 1e-6, eps_in, eps_out)

tbox, sigma, total_flux, n_val, n_iter = solve_single_box3d_adaptive_mesh(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, 1e-6, 1e-6, PointSource((5.0, 6.0, Lz / 2 + 0.01), 1.0))

# sigma_func = (x) -> sin(3 * x[1]) * cos(2 * x[2]) * exp(-x[3]^2)
# sigma = [sigma_func(pt.panel_point.point) for pt in BI.eachpoint(tbox)]

fmm_tol = 1e-10
hcubature_atol = 1e-8
range_factor = 5.0

DT_direct = BI.laplace3d_DT_fmm3d(tbox, fmm_tol)
DT_corrected_4 = BI.laplace3d_DT_fmm3d_corrected(tbox, fmm_tol, 1e-4, max_order, include_edges_src = false, include_edges_trg = false)
DT_corrected_8 = BI.laplace3d_DT_fmm3d_corrected(tbox, fmm_tol, 1e-8, 256, include_edges_src = false, include_edges_trg = false, range_factor = 6.0)

DT_hcubature_8 = BI.laplace3d_DT_fmm3d_corrected_hcubature(tbox, fmm_tol, hcubature_atol, 6.0, include_edges_src = false, include_edges_trg = false)

res_direct = DT_direct * sigma
res_corrected_4 = DT_corrected_4 * sigma
res_corrected_8 = DT_corrected_8 * sigma
res_hcubature_8 = DT_hcubature_8 * sigma

res_direct_approx = BI.interface_approx(tbox, res_direct)
res_corrected_4_approx = BI.interface_approx(tbox, res_corrected_4)
res_corrected_8_approx = BI.interface_approx(tbox, res_corrected_8)
res_hcubature_8_approx = BI.interface_approx(tbox, res_hcubature_8)

ll = Lx / 2

xs = range(-ll, ll, length = 200)
ys = range(-ll, ll, length = 200)
z_up = - Lz / 2

res_direct_val = [res_direct_approx((x, y, z_up)) for x in xs, y in ys]
res_corrected_4_val = [res_corrected_4_approx((x, y, z_up)) for x in xs, y in ys]
res_corrected_8_val = [res_corrected_8_approx((x, y, z_up)) for x in xs, y in ys]
res_hcubature_8_val = [res_hcubature_8_approx((x, y, z_up)) for x in xs, y in ys]

println("near correction (DT * sigma) max abs error vs hcubature apply:")
println("  direct:         ", maximum(abs.(res_direct_val .- res_hcubature_8_val)))
println("  corrected 1e-4: ", maximum(abs.((res_corrected_4_val .- res_hcubature_8_val))))
println("  corrected 1e-8: ", maximum(abs.((res_corrected_8_val .- res_hcubature_8_val))))

begin
    fig = Figure(size = (1000, 1000), fontsize = 20)
    ax = Axis(fig[1, 1], title = "hcubature apply value", aspect = DataAspect())
    hm1 = heatmap!(ax, xs, ys, log10.(abs.(res_hcubature_8_val)))
    Colorbar(fig[1, 2], hm1)

    ax = Axis(fig[1, 3], title = "error direct", aspect = DataAspect())
    hm2 = heatmap!(ax, xs, ys, log10.(abs.(res_direct_val .- res_hcubature_8_val)))
    Colorbar(fig[1, 4], hm2)

    ax = Axis(fig[2, 1], title = "error corrected 1e-4", aspect = DataAspect())
    hm3 = heatmap!(ax, xs, ys, log10.(abs.(res_corrected_4_val .- res_hcubature_8_val)))
    Colorbar(fig[2, 2], hm3)

    ax = Axis(fig[2, 3], title = "error corrected 1e-8", aspect = DataAspect())
    hm4 = heatmap!(ax, xs, ys, log10.(abs.(res_corrected_8_val .- res_hcubature_8_val)))
    Colorbar(fig[2, 4], hm4)

    fig
end

save("figs/near_correction_near_source.png", fig)

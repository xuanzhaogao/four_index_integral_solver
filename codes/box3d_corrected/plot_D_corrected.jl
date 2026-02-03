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


tbox = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, PointSource((5.0, 6.0, Lz / 2 + 0.01), 1.0), 1.0, l_ec, 1e-6, eps_in, eps_out)

fmm_tol = 1e-10
up_tol = 1e-8

D_direct = BI.laplace3d_DT_fmm3d(tbox, fmm_tol)
D_corrected_4 = BI.laplace3d_D_fmm3d_corrected(tbox, fmm_tol, 1e-4, max_order, include_edges_src = false, include_edges_trg = false)
D_corrected_8 = BI.laplace3d_D_fmm3d_corrected(tbox, fmm_tol, 1e-8, 256, include_edges_src = false, include_edges_trg = false)

res_direct = D_direct * ones(BI.num_points(tbox))
res_corrected_4 = D_corrected_4 * ones(BI.num_points(tbox))
res_corrected_8 = D_corrected_8 * ones(BI.num_points(tbox))

res_direct_approx = BI.interface_approx(tbox, res_direct)
res_corrected_4_approx = BI.interface_approx(tbox, res_corrected_4)
res_corrected_8_approx = BI.interface_approx(tbox, res_corrected_8)

ll = Lx / 2 

xs = range(-ll, ll, length = 200)
ys = range(-ll, ll, length = 200)
z_up = - Lz / 2

res_direct_val = [res_direct_approx((x, y, z_up)) for x in xs, y in ys]
res_corrected_4_val = [res_corrected_4_approx((x, y, z_up)) for x in xs, y in ys]
res_corrected_8_val = [res_corrected_8_approx((x, y, z_up)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 1000), fontsize = 20)
    ax = Axis(fig[1, 1], title = "direct apply value", aspect = DataAspect())
    hm1 = heatmap!(ax, xs, ys, log10.(abs.(res_direct_val .- 0.5)))
    Colorbar(fig[1, 2], hm1)

    ax = Axis(fig[1, 3], title = "error corrected 1e-4", aspect = DataAspect())
    hm2 = heatmap!(ax, xs, ys, log10.(abs.(res_corrected_4_val .- 0.5)))
    Colorbar(fig[1, 4], hm2)

    ax = Axis(fig[2, 1], title = "error corrected 1e-8", aspect = DataAspect())
    hm3 = heatmap!(ax, xs, ys, log10.(abs.(res_corrected_8_val .- 0.5)))
    Colorbar(fig[2, 2], hm3)

    fig
end

save("figs/D_mul_ones_err.png", fig)
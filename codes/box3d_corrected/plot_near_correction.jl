include("utils.jl")

using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

Lx = 10.0
Ly = 10.0
Lz = 1.0

p = 4
r = 3
eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^r * 1.01

na_box = BI.single_dielectric_box3d(Lx, Ly, Lz, p, l_ec, eps_in, eps_out)
na_D_direct = BI.laplace3d_D_fmm3d(na_box, 1e-6)
na_D_corrected_4 = BI.laplace3d_D_fmm3d_corrected(na_box, 1e-6, 1e-4, max_order, include_edges_src = false, include_edges_trg = false)
na_D_corrected_6 = BI.laplace3d_D_fmm3d_corrected(na_box, 1e-6, 1e-6, max_order, include_edges_src = false, include_edges_trg = false)
na_D_corrected_9 = BI.laplace3d_D_fmm3d_corrected(na_box, 1e-6, 1e-9, max_order, include_edges_src = false, include_edges_trg = false)

na_res_direct = na_D_direct * ones(BI.num_points(na_box))
na_res_corrected_4 = na_D_corrected_4 * ones(BI.num_points(na_box))
na_res_corrected_6 = na_D_corrected_6 * ones(BI.num_points(na_box))
na_res_corrected_9 = na_D_corrected_9 * ones(BI.num_points(na_box))

na_res_direct_approx = BI.interface_approx(na_box, na_res_direct)
na_res_corrected_4_approx = BI.interface_approx(na_box, na_res_corrected_4)
na_res_corrected_6_approx = BI.interface_approx(na_box, na_res_corrected_6)
na_res_corrected_9_approx = BI.interface_approx(na_box, na_res_corrected_9)


l = 4.0
xs = range(-l, l, length = 200)
ys = range(-l, l, length = 200)
z_up = 0.5

na_res_direct_val = [na_res_direct_approx((x, y, z_up)) for x in xs, y in ys]
na_res_corrected_4_val = [na_res_corrected_4_approx((x, y, z_up)) for x in xs, y in ys]
na_res_corrected_6_val = [na_res_corrected_6_approx((x, y, z_up)) for x in xs, y in ys]
na_res_corrected_9_val = [na_res_corrected_9_approx((x, y, z_up)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 1000), fontsize = 20)
    ax = Axis(fig[1, 1], title = "error direct", aspect = DataAspect())
    hm1 = heatmap!(ax, xs, ys, log10.(abs.(na_res_direct_val .- 0.5)))
    Colorbar(fig[1, 2], hm1)

    ax = Axis(fig[1, 3], title = "error corrected 1e-4", aspect = DataAspect())
    hm2 = heatmap!(ax, xs, ys, log10.(abs.(na_res_corrected_4_val .- 0.5)))
    Colorbar(fig[1, 4], hm2)

    ax = Axis(fig[2, 1], title = "error corrected 1e-6", aspect = DataAspect())
    hm3 = heatmap!(ax, xs, ys, log10.(abs.(na_res_corrected_6_val .- 0.5)))
    Colorbar(fig[2, 2], hm3)

    ax = Axis(fig[2, 3], title = "error corrected 1e-9", aspect = DataAspect())
    hm4 = heatmap!(ax, xs, ys, log10.(abs.(na_res_corrected_9_val .- 0.5)))
    Colorbar(fig[2, 4], hm4)

    fig
end

save("figs/plot_near_correction_D_mul_ones.png", fig)


ps = PointSource((0.1, 0.2, 10.0), 100.0)

tbox, sigma, total_flux, n_val, n_iter = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, 1e-6, 1e-6, 1e-6, max_order, ps)

DT_direct = BI.laplace3d_DT_fmm3d(tbox, 1e-6)
DT_corrected_4 = BI.laplace3d_DT_fmm3d_corrected(tbox, 1e-6, 1e-4, max_order, include_edges_src = false, include_edges_trg = false)
DT_corrected_6 = BI.laplace3d_DT_fmm3d_corrected(tbox, 1e-6, 1e-6, max_order, include_edges_src = false, include_edges_trg = false)

DT_hcubature = BI.laplace3d_DT_fmm3d_corrected_hcubature(tbox, 1e-8, 1e-8, 5.0)

res_direct = DT_direct * sigma
res_corrected_4 = DT_corrected_4 * sigma
res_corrected_6 = DT_corrected_6 * sigma
res_hcubature = DT_hcubature * sigma

res_direct_approx = BI.interface_approx(tbox, res_direct)
res_corrected_4_approx = BI.interface_approx(tbox, res_corrected_4)
res_corrected_6_approx = BI.interface_approx(tbox, res_corrected_6)
res_hcubature_approx = BI.interface_approx(tbox, res_hcubature)

xs = range(-4.5, 4.5, length = 200)
ys = range(-4.5, 4.5, length = 200)
z_up = 0.5

res_direct_val = [res_direct_approx((x, y, z_up)) for x in xs, y in ys]
res_corrected_4_val = [res_corrected_4_approx((x, y, z_up)) for x in xs, y in ys]
res_corrected_6_val = [res_corrected_6_approx((x, y, z_up)) for x in xs, y in ys]
res_hcubature_val = [res_hcubature_approx((x, y, z_up)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 1000), fontsize = 20)
    ax = Axis(fig[1, 1], title = "corrected value", aspect = DataAspect())
    hm1 = heatmap!(ax, xs, ys, log10.(abs.(res_hcubature_val)))
    Colorbar(fig[1, 2], hm1)

    ax = Axis(fig[1, 3], title = "error direct", aspect = DataAspect())
    hm2 = heatmap!(ax, xs, ys, log10.(abs.(res_direct_val .- res_hcubature_val)))
    Colorbar(fig[1, 4], hm2)

    ax = Axis(fig[2, 1], title = "error corrected 1e-4", aspect = DataAspect())
    hm3 = heatmap!(ax, xs, ys, log10.(abs.(res_corrected_4_val .- res_hcubature_val)))
    Colorbar(fig[2, 2], hm3)

    ax = Axis(fig[2, 3], title = "error corrected 1e-6", aspect = DataAspect())
    hm4 = heatmap!(ax, xs, ys, log10.(abs.(res_corrected_6_val .- res_hcubature_val)))
    Colorbar(fig[2, 4], hm4)

    fig
end

save("figs/plot_near_correction.png", fig)
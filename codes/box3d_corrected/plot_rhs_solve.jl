using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

L = 1.0
Lx = L
Ly = L
Lz = L

l_panel = 1.0
p = 6
r = 2
eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^4 * 1.01

ps = PointSource((0.1, 0.2, 1.5), 1.0)

interface_non = BI.single_dielectric_box3d(Lx, Ly, Lz, p, l_ec, eps_in, eps_out)
interface_refined_1e_6 = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, ps, 1.0, l_ec, 1e-6, eps_in, eps_out, Float64, max_depth = 12)
interface_refined_1e_9 = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, ps, 1.0, l_ec, 1e-9, eps_in, eps_out, Float64, max_depth = 12)

rhs_approx_func_non = BI.rhs_approx(interface_non, ps, 1.0)
rhs_approx_func_refined_1e_6 = BI.rhs_approx(interface_refined_1e_6, ps, 1.0)
rhs_approx_func_refined_1e_9 = BI.rhs_approx(interface_refined_1e_9, ps, 1.0)
rhs_func(p, n) =  -ps.charge * BI.laplace3d_grad(ps.point, p, n)

xs = range(-0.5, 0.5, length = 200)
ys = range(-0.5, 0.5, length = 200)
z = 0.5
n = (0.0, 0.0, 1.0)

vals_exact = [rhs_func((x, y, z), n) for x in xs, y in ys]
vals_approx_non = [rhs_approx_func_non((x, y, z)) for x in xs, y in ys]
vals_approx_refined_1e_6 = [rhs_approx_func_refined_1e_6((x, y, z)) for x in xs, y in ys]
vals_approx_refined_1e_9 = [rhs_approx_func_refined_1e_9((x, y, z)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 800), fontsize = 20)
    ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "exact rhs")
    ax2 = Axis(fig[1, 3], aspect = DataAspect(), title = "non-refined")
    ax3 = Axis(fig[2, 1], aspect = DataAspect(), title = "refined, atol = 1e-6")
    ax4 = Axis(fig[2, 3], aspect = DataAspect(), title = "refined, atol = 1e-9")

    hm1 = heatmap!(ax1, xs, ys, log10.(abs.(vals_exact)))
    Colorbar(fig[1, 2], hm1)

    hm2 = heatmap!(ax2, xs, ys, log10.(abs.(vals_approx_non .- vals_exact)))
    Colorbar(fig[1, 4], hm2)

    hm3 = heatmap!(ax3, xs, ys, log10.(abs.(vals_approx_refined_1e_6 .- vals_exact)))
    Colorbar(fig[2, 2], hm3)

    hm4 = heatmap!(ax4, xs, ys, log10.(abs.(vals_approx_refined_1e_9 .- vals_exact)))
    Colorbar(fig[2, 4], hm4)

    fig
end

save(joinpath(@__DIR__, "figs/rhs_solve_refined_d1.png"), fig)
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

ps = PointSource((0.1, 0.2, 0.6), 1.0)

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, ps, 1.0, l_ec, 1e-10, eps_in, eps_out, Float64, max_depth = 12)

rhs_approx_func = BI.rhs_approx(interface, ps, 1.0)
rhs_func(p, n) =  -ps.charge * BI.laplace3d_grad(ps.point, p, n)

xs = range(-0.5, 0.5, length = 200)
ys = range(-0.5, 0.5, length = 200)
z = 0.5
n = (0.0, 0.0, 1.0)

vals_exact = [rhs_func((x, y, z), n) for x in xs, y in ys]
vals_approx = [rhs_approx_func((x, y, z)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], aspect = DataAspect(), title = "log10(abs(exact rhs))")
    ax2 = Axis(fig[1, 3], aspect = DataAspect(), title = "log10(abs(error of approximate rhs))")

    hm1 = heatmap!(ax, xs, ys, log10.(abs.(vals_exact)))
    Colorbar(fig[1, 2], hm1)

    hm2 = heatmap!(ax2, xs, ys, log10.(abs.(vals_approx .- vals_exact)))
    Colorbar(fig[1, 4], hm2)

    fig
end

save(joinpath(@__DIR__, "figs/rhs_solve_10_digits.png"), fig)
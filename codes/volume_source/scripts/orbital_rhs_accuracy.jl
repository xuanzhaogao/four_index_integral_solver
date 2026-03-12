using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using CairoMakie

orbital_file = joinpath(@__DIR__, "../../../density_data/graphene_00001_5x5x1_shifted.xsf")
structure, datagrid = BI.read_xsf(orbital_file)

datagrid.values .*= datagrid.values

vs = BoundaryIntegral.VolumeSource(datagrid, shift = (0.0, 0.0, - 7.920155482424242), tol = 1e-4)


L = 90.0
Lx = L
Ly = L
Lz = 2.4

l_panel = 1.0
p = 6
eps_out = 1.0
eps_in = 6.0

max_order = 128
l_ec = 10.0 / 2^4 * 1.01

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, vs, 1.0, l_ec, 1e-4, eps_in, eps_out, Float64)

xs = range(-5 / 2, stop = 5 / 2, length = 400)
ys = range(-5 / 2, stop = 5 / 2, length = 400)
z = Lz / 2

rhs_hybrid = BI.Rhs_dielectric_box3d_hybrid(interface, vs, 1.0, 1e-4)
rhs_hybrid_approx = BI.interface_approx(interface, rhs_hybrid)
rhs_hybrid_values = [rhs_hybrid_approx((x, y, z)) for x in xs, y in ys]

rhs_hybrid_high = BI.Rhs_dielectric_box3d_hybrid(interface, vs, 1.0, 1e-6)
rhs_hybrid_high_approx = BI.interface_approx(interface, rhs_hybrid_high)
rhs_hybrid_high_values = [rhs_hybrid_high_approx((x, y, z)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 400), fontsize = 20)

    ax_1 = Axis(fig[1, 1], title = "value", xlabel = "x", ylabel = "y", aspect = DataAspect(),)
    hm_1 = heatmap!(ax_1, xs, ys,rhs_hybrid_values, colormap = :bluesreds)
    Colorbar(fig[1, 2], hm_1)

    ax_2 = Axis(fig[1, 3], title = "abs error", xlabel = "x", ylabel = "y", aspect = DataAspect(),)
    hm_2 = heatmap!(ax_2, xs, ys, abs.(rhs_hybrid_values .- rhs_hybrid_high_values), colormap = :viridis)
    Colorbar(fig[1, 4], hm_2)
    
    fig
end

save("figs/orbital_rhs_accuracy.png", fig)
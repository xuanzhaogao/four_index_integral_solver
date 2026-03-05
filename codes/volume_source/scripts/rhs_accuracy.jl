using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, SpecialFunctions, ForwardDiff
using CairoMakie

L = 20.0
Lx = L
Ly = L
Lz = 1.0

l_panel = 1.0
p = 6
r = 3
eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^r * 1.01

# Gaussian volume source on a tensor grid
center = (1.0, 2.0, 0.51)
σ = 0.5

# tol controls truncation of the Gaussian support
source = BI.GaussianVolumeSource(center, σ, 60, 1e-8)

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, source, 1.0, l_ec, 1e-6, eps_in, eps_out, Float64)

xs = range(-Lx / 2, stop = Lx / 2, length = 400)
ys = range(-Ly / 2, stop = Ly / 2, length = 400)
z = Lz / 2


rhs_exact = BI.Rhs_dielectric_box3d_gaussian(interface, center, σ, 1.0)
rhs_exact_approx = BI.interface_approx(interface, rhs_exact)
rhs_exact_values = [rhs_exact_approx((x, y, z)) for x in xs, y in ys]

rhs_fmm = BI.Rhs_dielectric_box3d_fmm3d(interface, source, 1.0, 1e-8)
rhs_fmm_approx = BI.interface_approx(interface, rhs_fmm)
rhs_fmm_values = [rhs_fmm_approx((x, y, z)) for x in xs, y in ys]

rhs_hybrid = BI.Rhs_dielectric_box3d_hybrid(interface, source, 1.0, 1e-9, fbc_N = 128)
rhs_hybrid_approx = BI.interface_approx(interface, rhs_hybrid)
rhs_hybrid_values = [rhs_hybrid_approx((x, y, z)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 400), fontsize = 20)

    ax_1 = Axis(fig[1, 1], title = "abs error hybrid", xlabel = "x", ylabel = "y", aspect = DataAspect(),)
    hm_1 = heatmap!(ax_1, xs, ys, log10.(abs.(rhs_hybrid_values .- rhs_exact_values) ./ abs.(rhs_exact_values)), colormap = :viridis)
    Colorbar(fig[1, 2], hm_1)

    ax_2 = Axis(fig[1, 3], title = "abs error fmm", xlabel = "x", ylabel = "y", aspect = DataAspect(),)
    hm_2 = heatmap!(ax_2, xs, ys, log10.(abs.(rhs_fmm_values .- rhs_exact_values) ./ abs.(rhs_exact_values)), colormap = :viridis)
    Colorbar(fig[1, 4], hm_2)

    save("figs/rhs_accuracy.png", fig)

    fig
end
using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, SpecialFunctions, ForwardDiff
using CairoMakie

# analytic solution of the potential is erf(x / (sqrt(2) * sigma)) / x
function analytic_potential(center, target, sigma)
    r = norm(target .- center)
    return erf(r / (sqrt(2) * sigma)) / r
end

function analytic_DT(center, target, sigma, normal)
    f = r -> erf(r[1] / (sqrt(2) * sigma)) / r[1]
    r_vec = target .- center
    r = norm(r_vec)
    grad_f = ForwardDiff.gradient(f, [r])

    return dot(grad_f .* (r_vec ./ r), normal)
end

function rhs_analytic(interface, center, sigma)
    rhs_analytic = zeros(BI.num_points(interface))
    for (i, p) in enumerate(BI.eachpoint(interface))
        target = p.panel_point.point
        normal = p.panel_point.normal
        rhs_analytic[i] = analytic_DT(center, target, sigma, normal)
    end
    return rhs_analytic ./ 4π
end

L = 20.0
Lx = L
Ly = L
Lz = 1.0

l_panel = 1.0
p = 4
r = 3
eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^r * 1.01

# Gaussian volume source on a tensor grid
center = (1.0, 2.0, 2.0)
σ = 0.1

# tol controls truncation of the Gaussian support
source = BI.GaussianVolumeSource(center, σ, 1e-6)

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, source, 1.0, l_ec, 1e-6, eps_in, eps_out, Float64)

# rhs = BI.Rhs_dielectric_box3d(interface, source, 1.0)
rhs_fmm = BI.Rhs_dielectric_box3d_fmm3d(interface, source, 1.0, 1e-12)
rhs_exact = rhs_analytic(interface, center, σ)

rhs_fmm_approx = BI.interface_approx(interface, rhs_fmm)
rhs_exact_approx = BI.interface_approx(interface, rhs_exact)

xs = range(-Lx / 2, stop = Lx / 2, length = 200)
ys = range(-Ly / 2, stop = Ly / 2, length = 200)
z = Lz / 2

rhs_fmm_values = [rhs_fmm_approx((x, y, z)) for x in xs, y in ys]
rhs_exact_values = [rhs_exact_approx((x, y, z)) for x in xs, y in ys]

begin
    fig = Figure(size = (1000, 400), fontsize = 20)

    ax_1 = Axis(fig[1, 1], title = "Analytic value", xlabel = "x", ylabel = "y", aspect = DataAspect(),)
    hm_1 = heatmap!(ax_1, xs, ys, log10.(abs.(rhs_exact_values)), colormap = :plasma)
    Colorbar(fig[1, 2], hm_1)

    ax_2 = Axis(fig[1, 3], title = "error", xlabel = "x", ylabel = "y", aspect = DataAspect(),)
    hm_2 = heatmap!(ax_2, xs, ys, log10.(abs.(rhs_fmm_values .- rhs_exact_values)), colormap = :plasma)
    Colorbar(fig[1, 4], hm_2)

    save("figs/rhs_accuracy.png", fig)

    fig
end
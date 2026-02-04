# first generate the using BoundaryIntegral
include("utils.jl")
using CairoMakie

L = 5.0
Lx = L
Ly = L
Lz = 1.0

p = 6
r = 2
eps_out = 1.0
eps_in = 4.0

fmm_tol = 1e-4
up_tol = 1e-5
max_order = 128
l_ec = 1 / 2^4 * 1.01

ps = PointSource((0.1, 0.2, 0.51), 1.0)

tbox_adaptive, sigma_adaptive, total_flux_adaptive, n_val_adaptive, n_iter_adaptive = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, up_tol, max_order, ps)

tbox_fine, sigma_fine, total_flux_fine, n_val_fine, n_iter_fine = solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, p + 2, l_ec, eps_in, eps_out, fmm_tol, up_tol, up_tol, max_order, ps)

rhs_function = BI.rhs_approx(tbox_adaptive, ps, 1.0)
sigma_function = BI.interface_approx(tbox_adaptive, sigma_adaptive)
sigma_fine_function = BI.interface_approx(tbox_fine, sigma_fine)

ll = 1.0

xs = range(-ll, ll, length = 200)
ys = range(-ll, ll, length = 200)

z_up = 0.5
z_down = - 0.5

rhs_up = [rhs_function((x, y, z_up)) for x in xs, y in ys]
rhs_down = [rhs_function((x, y, z_down)) for x in xs, y in ys]

sigma_up = [sigma_function((x, y, z_up)) for x in xs, y in ys]
sigma_down = [sigma_function((x, y, z_down)) for x in xs, y in ys]

sigma_fine_up = [sigma_fine_function((x, y, z_up)) for x in xs, y in ys]
sigma_fine_down = [sigma_fine_function((x, y, z_down)) for x in xs, y in ys]

begin
    fig = Figure(size = (1200, 800), fontsize = 20)
    ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "rhs up")
    ax2 = Axis(fig[2, 1], aspect = DataAspect(), title = "rhs down")
    ax3 = Axis(fig[1, 3], aspect = DataAspect(), title = "sigma up")
    ax4 = Axis(fig[2, 3], aspect = DataAspect(), title = "sigma down")
    ax5 = Axis(fig[1, 5], aspect = DataAspect(), title = "error sigma up")
    ax6 = Axis(fig[2, 5], aspect = DataAspect(), title = "error sigma down")

    hm1 = heatmap!(ax1, xs, ys, log10.(abs.(rhs_up)))
    Colorbar(fig[1, 2], hm1)

    hm2 = heatmap!(ax2, xs, ys, log10.(abs.(rhs_down)))
    Colorbar(fig[2, 2], hm2)

    hm3 = heatmap!(ax3, xs, ys, log10.(abs.(sigma_up)))
    Colorbar(fig[1, 4], hm3)

    hm4 = heatmap!(ax4, xs, ys, log10.(abs.(sigma_down)))
    Colorbar(fig[2, 4], hm4)

    hm5 = heatmap!(ax5, xs, ys, log10.(abs.(sigma_up .- sigma_fine_up)))
    Colorbar(fig[1, 6], hm5)

    hm6 = heatmap!(ax6, xs, ys, log10.(abs.((sigma_down .- sigma_fine_down))))
    Colorbar(fig[2, 6], hm6)

    fig
end

save("figs/sigma_thin_box.png", fig)

ll = L / 2

xs = range(-ll, ll, length = 400)
ys = range(-ll, ll, length = 400)

z_up = 0.5
z_down = - 0.5

rhs_up = [rhs_function((x, y, z_up)) for x in xs, y in ys]
rhs_down = [rhs_function((x, y, z_down)) for x in xs, y in ys]

sigma_up = [sigma_function((x, y, z_up)) for x in xs, y in ys]
sigma_down = [sigma_function((x, y, z_down)) for x in xs, y in ys]

sigma_fine_up = [sigma_fine_function((x, y, z_up)) for x in xs, y in ys]
sigma_fine_down = [sigma_fine_function((x, y, z_down)) for x in xs, y in ys]


# make 3d plot of rhs and sigma
begin
    fig3d = Figure(size = (1500, 800), fontsize = 16)
    
    # RHS plots
    ax3d_rhs_up = Axis3(fig3d[1, 1], title = "rhs up (z = 0.5)", xlabel = "x", ylabel = "y", zlabel = "log10(|rhs|)")
    ax3d_rhs_down = Axis3(fig3d[2, 1], title = "rhs down (z = -0.5)", xlabel = "x", ylabel = "y", zlabel = "log10(|rhs|)")
    
    # Sigma plots
    ax3d_sigma_up = Axis3(fig3d[1, 2], title = "sigma up (z = 0.5)", xlabel = "x", ylabel = "y", zlabel = "log10(|σ|)")
    ax3d_sigma_down = Axis3(fig3d[2, 2], title = "sigma down (z = -0.5)", xlabel = "x", ylabel = "y", zlabel = "log10(|σ|)")

    ax3d_error_sigma_up = Axis3(fig3d[1, 3], title = "error sigma up (z = 0.5)", xlabel = "x", ylabel = "y", zlabel = "log10(|σ - σ_fine|)")
    ax3d_error_sigma_down = Axis3(fig3d[2, 3], title = "error sigma down (z = -0.5)", xlabel = "x", ylabel = "y", zlabel = "log10(|σ - σ_fine|)")
    
    # Surface plots for RHS
    surface!(ax3d_rhs_up, xs, ys, log10.(abs.(rhs_up)), colormap = :viridis)
    surface!(ax3d_rhs_down, xs, ys, log10.(abs.(rhs_down)), colormap = :viridis)
    
    # Surface plots for Sigma
    surface!(ax3d_sigma_up, xs, ys, log10.(abs.(sigma_up)), colormap = :viridis)
    surface!(ax3d_sigma_down, xs, ys, log10.(abs.(sigma_down)), colormap = :viridis)

    # Surface plots for Error Sigma
    surface!(ax3d_error_sigma_up, xs, ys, log10.(abs.(sigma_up .- sigma_fine_up)), colormap = :viridis)
    surface!(ax3d_error_sigma_down, xs, ys, log10.(abs.(sigma_down .- sigma_fine_down)), colormap = :viridis)

    fig3d
end

save("figs/sigma_thin_box_3d.png", fig3d)

using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

function error_heatmap(nbox)
    nt = 300
    x = 0.0
    y = 0.0
    x_end = 1.0 - 1e-6
    y_end = 1.0 - 1e-6
    locs = zeros(3, nt * nt)
    xs = range(x, x_end, length = nt)
    ys = range(y, y_end, length = nt)
    for i in 1:nt
        for j in 1:nt
            locs[:, (i - 1) * nt + j] = [xs[i], ys[j], 1.0]
        end
    end

    D_trgs = BI.laplace3d_D_trg_fmm3d(nbox, locs, 1e-14)
    gi_trgs = D_trgs * ones(size(D_trgs, 2))
    log10_err_trgs = log10.(abs.(gi_trgs .+ 0.5))
    # log10_err_trgs = log10.(abs.(gi_trgs))
    return log10_err_trgs
end

# nbox = BI.dielectric_box3d(10.0, 1.0, 2, ns, true, 8, 8)
box_nonadapt_64 = BI.dielectric_box3d(10.0, 1.0, 1, 32, 32, 0, 0)
box_nonadapt_128 = BI.dielectric_box3d(10.0, 1.0, 1, 128, 128, 0, 0)
box_adapt_nr_4 = BI.dielectric_box3d(10.0, 1.0, 1, 16, 16, 4, 4)
box_adapt_nr_8 = BI.dielectric_box3d(10.0, 1.0, 1, 16, 16, 8, 8)
box_adapt_r_4 = BI.dielectric_box3d(10.0, 1.0, 1, 16, 4, 4, 4)
box_adapt_r_8 = BI.dielectric_box3d(10.0, 1.0, 1, 16, 4, 8, 8)

BI.num_points(box_nonadapt_64)
BI.num_points(box_nonadapt_128)
BI.num_points(box_adapt_nr_4)
BI.num_points(box_adapt_nr_8)
BI.num_points(box_adapt_r_4)
BI.num_points(box_adapt_r_8)

log10_err_trgs_nonadapt_64 = error_heatmap(box_nonadapt_64)
log10_err_trgs_nonadapt_128 = error_heatmap(box_nonadapt_128)
log10_err_trgs_adapt_nr_4 = error_heatmap(box_adapt_nr_4)
log10_err_trgs_adapt_nr_8 = error_heatmap(box_adapt_nr_8)
log10_err_trgs_adapt_r_4 = error_heatmap(box_adapt_r_4)
log10_err_trgs_adapt_r_8 = error_heatmap(box_adapt_r_8)

D_trg_nonadapt_64_direct = BI.laplace3d_D_trg(box_nonadapt_64, locs)
gi_nonadapt_64_direct = D_trg_nonadapt_64_direct * ones(size(D_trg_nonadapt_64_direct, 2))
log10_err_trgs_nonadapt_64_direct = log10.(abs.(gi_nonadapt_64_direct .- 0.5))

begin
    nt = 300
    x = 0.0
    y = 0.0
    x_end = 1.0 - 1e-6
    y_end = 1.0 - 1e-6
    locs = zeros(3, nt * nt)
    xs = range(x, x_end, length = nt)
    ys = range(y, y_end, length = nt)
    for i in 1:nt
        for j in 1:nt
            locs[:, (i - 1) * nt + j] = [xs[i], ys[j], 1.0]
        end
    end

    fig = Figure(size = (1200, 1800), fontsize = 20, dpi = 500)
    ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "non-adaptive p = 32", aspect=DataAspect())
    ax2 = Axis(fig[1, 3], xlabel = "x", ylabel = "y", title = "non-adaptive p = 128", aspect=DataAspect())
    ax3 = Axis(fig[2, 1], xlabel = "x", ylabel = "y", title = "adaptive non-reduced p = 16 r = 4", aspect=DataAspect())
    ax4 = Axis(fig[2, 3], xlabel = "x", ylabel = "y", title = "adaptive non-reduced p = 16 r = 8", aspect=DataAspect())
    ax5 = Axis(fig[3, 1], xlabel = "x", ylabel = "y", title = "adaptive reduced p = 16 r = 4", aspect=DataAspect())
    ax6 = Axis(fig[3, 3], xlabel = "x", ylabel = "y", title = "adaptive reduced p = 16 r = 8", aspect=DataAspect())
    hm1 = heatmap!(ax1, locs[1, :], locs[2, :], log10_err_trgs_nonadapt_64; colormap = :viridis, interpolate = 
    false)
    hm2 = heatmap!(ax2, locs[1, :], locs[2, :], log10_err_trgs_nonadapt_128; colormap = :viridis, interpolate = false)
    hm3 = heatmap!(ax3, locs[1, :], locs[2, :], log10_err_trgs_adapt_nr_4; colormap = :viridis, interpolate = false)
    hm4 = heatmap!(ax4, locs[1, :], locs[2, :], log10_err_trgs_adapt_nr_8; colormap = :viridis, interpolate = false)
    hm5 = heatmap!(ax5, locs[1, :], locs[2, :], log10_err_trgs_adapt_r_4; colormap = :viridis, interpolate = false)
    hm6 = heatmap!(ax6, locs[1, :], locs[2, :], log10_err_trgs_adapt_r_8; colormap = :viridis, interpolate = false)

    Colorbar(fig[1, 2], hm1, ticks = -15:1:2)
    Colorbar(fig[1, 4], hm2, ticks = -15:1:2)
    Colorbar(fig[2, 2], hm3, ticks = -15:1:2)
    Colorbar(fig[2, 4], hm4, ticks = -15:1:2)
    Colorbar(fig[3, 2], hm5, ticks = -15:1:2)
    Colorbar(fig[3, 4], hm6, ticks = -15:1:2)
    # xlims!(ax1, x, x_end)
    # ylims!(ax1, y, y_end)

    save(joinpath(@__DIR__, "figs/single_box3d_gi_heatmap.png"), fig)

    fig
end

begin
    fig = Figure(size = (1200, 500), fontsize = 20, dpi = 500)
    ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "Direct eval", aspect=DataAspect())
    ax2 = Axis(fig[1, 3], xlabel = "x", ylabel = "y", title = "FMM3D eval", aspect=DataAspect())

    hm1 = heatmap!(ax1, locs[1, :], locs[2, :], log10_err_trgs_nonadapt_64_direct; colormap = :viridis, interpolate = false)
    hm2 = heatmap!(ax2, locs[1, :], locs[2, :], log10_err_trgs_nonadapt_64; colormap = :viridis, interpolate = false)

    Colorbar(fig[1, 2], hm1, ticks = -15:1:2)
    Colorbar(fig[1, 4], hm2, ticks = -15:1:2)

    save(joinpath(@__DIR__, "figs/single_box3d_gi_heatmap_direct.png"), fig)

    fig
end
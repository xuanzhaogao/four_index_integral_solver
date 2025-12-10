using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

ns = 256

# nbox = BI.dielectric_box3d(10.0, 1.0, 2, ns, true, 8, 8)
nbox = BI.dielectric_box3d(10.0, 1.0, 1, 8, false, 12, 12)
n = BI.num_points(nbox)

# BI.viz_3d_dielectric_interfaces(nbox)

D = BI.laplace3d_D_fmm3d(nbox, 1e-6)
o = ones(size(D, 1))
gi = D * o

begin
    log10_err = log10.(abs.(gi .+ 0.5))
    point_xs = zeros(Int(n / 6))
    point_ys = zeros(Int(n / 6))
    vals = zeros(Int(n / 6))
    offset = 1
    for (i, p) in enumerate(BI.eachpoint(nbox))
        if p.point[3] == 1.0
            point_xs[offset] = p.point[1]
            point_ys[offset] = p.point[2]
            vals[offset] = log10_err[i]
            offset += 1
        end
    end

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "log10 error")
    hm = heatmap!(ax, point_xs, point_ys, vals; colormap = :viridis, interpolate = true)
    Colorbar(fig[1, 2], hm)
    fig
end


begin
    nt = 100
    x = 0.999
    y = 0.999
    locs = zeros(3, nt * nt)
    xs = range(x, 1.0, length = nt)
    ys = range(y, 1.0, length = nt)
    for i in 1:nt
        for j in 1:nt
            locs[:, (i - 1) * nt + j] = [xs[i], ys[j], 1.0]
        end
    end

    D_trgs = BI.laplace3d_D_trg_fmm3d(nbox, locs, 1e-14)
    gi_trgs = D_trgs * ones(size(D_trgs, 2))
    log10_err_trgs = log10.(abs.(gi_trgs .+ 0.5))
    # log10_err_trgs = log10.(abs.(gi_trgs))
end

begin
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "log10 error", aspect=DataAspect())
    hm = heatmap!(ax, locs[1, :], locs[2, :], log10_err_trgs; colormap = :viridis, interpolate = true)
    Colorbar(fig[1, 2], hm, ticks = -15:1:2)
    xlims!(ax, x, 1.0)
    ylims!(ax, y, 1.0)
    fig
end
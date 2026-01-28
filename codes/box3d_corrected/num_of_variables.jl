# to create a heatmap of the number of variables for different p and r

using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

function count_variables(p, r, up_tol)
    l_ec = 1.0 / 2^r * 1.01

    tbox = BI.single_dielectric_box3d_rhs_adaptive(1.0, 1.0, 1.0, p, PointSource((0.2, 0.3, 0.55), 1.0), 1.0, l_ec, up_tol, 4.0, 1.0, max_depth = 1000)

    return BI.num_points(tbox)
end

ps = collect(4:10)
rs = collect(0:1:8)

num_vars = [count_variables(p, r, 1e-6) for p in ps, r in rs]

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "p", ylabel = "l_min", title = "number of variables", yscale = log2)
    hm = heatmap!(ax, ps, 1 ./(2 .^ rs), log10.(num_vars))
    Colorbar(fig[1, 2], hm)

    fig
end

save(joinpath(@__DIR__, "figs/num_of_variables.png"), fig)
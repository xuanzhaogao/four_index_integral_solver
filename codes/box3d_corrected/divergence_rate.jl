# use a cubic box excited by a far away point source
# see the divergence rate of the edge singularities

include("utils.jl")
# using CairoMakie
using LsqFit

df = joinpath(@__DIR__, "data/divergence_rate.csv")
CSV.write(df, DataFrame(eps_in = Float64[], x = Float64[], slope = Float64[]))

L = 1.0
l = 0.0101
p = 6
atol = 1e-6
max_order = 128
source = PointSource((10.0, 10.0, 10.0), 10000.0)

for eps_in in range(2.0, 20.0, length = 10)
    box, sigma, _, _, _ = solve_single_box3d_adaptive_mesh_corrected(L, L, L, p, l, eps_in, 1.0, atol, atol, atol, max_order, source, include_edges_src = false, include_edges_trg = false)

    sigma_approx = BI.interface_approx(box, sigma)

    xs = range(-0.49, 0, length = 10)
    ys = range(0.48, 0.49, length = 100)
    z = 0.5

    vals = [[sigma_approx((x, y, z)) for y in ys] for x in xs]

    for (i, x) in enumerate(xs)
        val = vals[i]
        # fit a line to log-log data
        log_x = log10.(0.5 .- ys)
        log_val = log10.(abs.(val))
        model(x, p) = p[1] * x .+ p[2]
        fit = curve_fit(model, log_x, log_val, [1.0, 0.0])

        CSV.write(df, DataFrame(eps_in = eps_in, x = x, slope = fit.param[1]), append = true)
    end
end

# begin
#     fig = Figure(size = (600, 400), fontsize = 16)
#     ax = Axis(fig[1, 1], xlabel = "y", yscale = log10, xscale = log10)

#     for (i, x) in enumerate(xs)
#         lines!(ax, abs.(0.5 .- ys), abs.(vals[i]), label = "x=$(round(x, digits=2))")
#     end

#     Legend(fig[1, 2], ax, orientation = :vertical)

#     fig
# end

# save(joinpath(@__DIR__, "figs/divergence_rate.svg"), fig)
using CSV, DataFrames
using CairoMakie, LaTeXStrings
using LinearAlgebra

df = CSV.read("data/orbital_source_convergence.csv", DataFrame)


ps = [4, 6]
rs = collect(0:2:6)

l_ecs = 9.0 ./ 2 .^ rs * 1.01

ref_p = 6
ref_r = 7


ids = unique(df.loc_id)

colors = [:red, :blue, :green, :orange, :purple]
markers = [:circle, :diamond, :utriangle, :dtriangle]

begin
    fig = Figure(size = (500, 400), fontsize = 16)
    ax = Axis(fig[1, 1], xlabel = "r", ylabel = "relative error", yscale = log10, xscale = log10, xreversed = true)

    ref_res = df[(df.p .== ref_p) .& (df.r .== ref_r), :potential]

    for (i, p) in enumerate(ps)
        err_p = Float64[]
        for r in rs
            pot_r = df[(df.p .== p) .& (df.r .== r), :potential]
            err = norm(pot_r .- ref_res) / norm(ref_res)
            @show err
            push!(err_p, err)
        end
        @show err_p

        scatterlines!(ax, l_ecs, err_p, label = "p = $p", color = colors[i], marker = markers[i])
    end

    axislegend(position = :rt)

    fig

end

save("figs/orbital_source_convergence.png", fig)
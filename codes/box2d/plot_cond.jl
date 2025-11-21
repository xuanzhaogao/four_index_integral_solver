using CSV, DataFrames, CairoMakie, LaTeXStrings

df = CSV.read("data/condition_number_threaded.csv", DataFrame)

n_adapts = unique(df.n_adapt)
epses = unique(df.eps)[2:end]

conds = [zeros(length(epses)) for _ in eachindex(n_adapts)]
iters = [zeros(length(epses)) for _ in eachindex(n_adapts)]
for i in eachindex(n_adapts)
    for j in eachindex(epses)
        conds[i][j] = df[(df.n_adapt .== n_adapts[i]) .&& (df.eps .== epses[j]), :cond][1]
        iters[i][j] = df[(df.n_adapt .== n_adapts[i]) .&& (df.eps .== epses[j]), :niter][1]
    end
end

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = L"$\epsilon$", ylabel = "Condition number")
    for i in eachindex(n_adapts)
        lines!(ax, epses, conds[i], label = "$(n_adapts[i])")
    end

    ax2 = Axis(fig[1, 2], xlabel = L"$\epsilon$", ylabel = "GMRES iterations", xscale = log10)
    for i in eachindex(n_adapts)
        lines!(ax2, epses, iters[i], label = "$(n_adapts[i])")
    end

    Legend(fig[1, 3], ax, "N adpat", position = :lt, nbanks = 2)
end

fig

save("figs/condition_number.svg", fig)
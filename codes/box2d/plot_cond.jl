using CSV, DataFrames, CairoMakie, LaTeXStrings

df = CSV.read("data/condition_number_threaded.csv", DataFrame)

n_adapts = unique(df.n_adapt)
gammas = unique(df.gamma)

conds = [zeros(length(gammas)) for _ in eachindex(n_adapts)]
iters = [zeros(length(gammas)) for _ in eachindex(n_adapts)]
for i in eachindex(n_adapts)
    for j in eachindex(gammas)
        conds[i][j] = df[(df.n_adapt .== n_adapts[i]) .&& (df.gamma .== gammas[j]), :cond][1]
        iters[i][j] = df[(df.n_adapt .== n_adapts[i]) .&& (df.gamma .== gammas[j]), :niter][1]
    end
end

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], yscale = log10, xlabel = L"$\gamma$", ylabel = "Condition number")
    for i in eachindex(n_adapts)
        lines!(ax, gammas, conds[i], label = "$(n_adapts[i])")
    end

    ax2 = Axis(fig[1, 2], xlabel = L"$\gamma$", ylabel = "GMRES iterations")
    for i in eachindex(n_adapts)
        lines!(ax2, gammas, iters[i], label = "$(n_adapts[i])")
    end

    Legend(fig[1, 3], ax, "N adpat", position = :lt, nbanks = 2)
end

fig

save("figs/condition_number.svg", fig)
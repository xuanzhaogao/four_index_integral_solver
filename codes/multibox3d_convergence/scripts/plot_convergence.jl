using CairoMakie
using DelimitedFiles
using LaTeXStrings

const EXPECTED_COLUMNS = (
    :n_quad,
    :edge_refine,
    :n_pts,
    :n_iter,
    :residual,
    :U_00,
    :U_01,
    :U_02,
    :U_03,
)

const U_COLUMNS = (:U_00, :U_01, :U_02, :U_03)

parse_number(x) = x isa Number ? Float64(x) : parse(Float64, strip(String(x)))

function read_convergence_table(path::AbstractString)
    isfile(path) || error("Missing input file: $path")

    raw, header = readdlm(path, ',', Any, '\n'; header = true)
    names = Symbol.(strip.(String.(vec(header))))
    missing_columns = setdiff(EXPECTED_COLUMNS, names)
    isempty(missing_columns) || error("Missing expected columns: $(join(string.(missing_columns), ", "))")

    indices = Dict(name => findfirst(==(name), names) for name in names)
    table = Dict{Symbol, Vector{Float64}}()
    for name in EXPECTED_COLUMNS
        j = indices[name]
        table[name] = [parse_number(raw[i, j]) for i in axes(raw, 1)]
    end

    return table
end

function sorted_positive_xy(x::AbstractVector, y::AbstractVector)
    mask = (x .> 0) .& (y .> 0)
    xs = x[mask]
    ys = y[mask]
    order = sortperm(xs)
    return xs[order], ys[order]
end

function main()
    root_dir = @__DIR__
    data_path = joinpath(root_dir, "data", "convergence.csv")
    fig_dir = joinpath(root_dir, "figs")
    mkpath(fig_dir)

    table = read_convergence_table(data_path)

    n_quad = Int.(round.(table[:n_quad]))
    edge_refine = Int.(round.(table[:edge_refine]))
    n_pts = Int.(round.(table[:n_pts]))

    ref_idx = length(n_pts)
    if n_quad[ref_idx] != 6 || edge_refine[ref_idx] != 4
        @warn "Last row is not (n_quad=6, edge_refine=4); using the last row as reference anyway."
    end

    errors = Dict(name => abs.(table[name] .- table[name][ref_idx]) for name in U_COLUMNS)
    max_errors = [maximum(errors[name][i] for name in U_COLUMNS) for i in eachindex(n_pts)]

    CairoMakie.activate!()

    fig = with_theme(theme_latexfonts()) do
        Figure(size = (1000, 500), fontsize = 16)
    end

    ax1 = Axis(
        fig[1, 1],
        xscale = log10,
        yscale = log10,
        xlabel = L"N_{\text{pts}}",
        ylabel = "Absolute Error",
        title = "Convergence vs discretization size",
    )

    colors = Makie.wong_colors()[1:4]
    for (i, name) in enumerate(U_COLUMNS)
        xs, ys = sorted_positive_xy(n_pts, errors[name])
        scatterlines!(
            ax1,
            xs,
            ys;
            marker = :circle,
            markersize = 10,
            linewidth = 2,
            color = colors[i],
            label = String(name),
        )
    end
    axislegend(ax1, position = :rb)

    ax2 = Axis(
        fig[1, 2],
        xscale = log10,
        yscale = log10,
        xlabel = L"N_{\text{pts}}",
        ylabel = "Max Absolute Error",
        title = "Convergence by quadrature order",
    )

    markers = Dict(2 => :circle, 4 => :rect, 6 => :utriangle)
    quad_colors = Dict(2 => colors[1], 4 => colors[2], 6 => colors[3])

    for p in (2, 4, 6)
        mask = (n_quad .== p) .& (max_errors .> 0)
        xs, ys = sorted_positive_xy(n_pts[mask], max_errors[mask])
        scatterlines!(
            ax2,
            xs,
            ys;
            marker = markers[p],
            markersize = 10,
            linewidth = 2,
            color = quad_colors[p],
            label = "p=$p",
        )
    end
    axislegend(ax2, position = :rb)

    colgap!(fig.layout, 24)

    save(joinpath(fig_dir, "convergence.svg"), fig)
    save(joinpath(fig_dir, "convergence.png"), fig)
    println("Saved to figs/convergence.{svg,png}")
end

main()

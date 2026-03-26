using BoundaryIntegral
import BoundaryIntegral as BI

using CairoMakie
using LaTeXStrings

include(joinpath(@__DIR__, "..", "src", "ScreenedDensityAnalysis.jl"))
using .ScreenedDensityAnalysis

const ORBITAL_FILE = default_orbital_file()
const TOL = 1e-4
const L = 90.0
const LZ = 2.8
const EPS_IN = 2.4
const EPS_OUT = 1.0
const BANDWIDTHS = [0.05, 0.1, 0.2, 0.4]
const QZ_UPSAMPLE_FACTOR = 8
const QZ_GLOBAL_PAD_FACTOR = 4

function analysis_cases()
    cases = Tuple{String, BI.AbstractScreeningMode}[("Sharp", BI.SharpScreening())]
    for bandwidth in BANDWIDTHS
        push!(cases, ("SoftMix b = $(bandwidth)", BI.SoftMixInversePermittivity(bandwidth)))
    end
    return cases
end

function save_slice_figure(datagrid, bounds, cases)
    slice_data = [
        (label, yz_slice_at_x(datagrid, 0.0, mode; bounds = bounds, eps_in = EPS_IN, eps_out = EPS_OUT, ny = 220, nz = 260, tol = TOL))
        for (label, mode) in cases
    ]

    positive_values = Float64[]
    for (_, (_, _, slice)) in slice_data
        append!(positive_values, filter(>(0.0), filter(isfinite, vec(slice))))
    end
    floor_value = minimum(positive_values)
    log_slice_data = [
        (label, ys, zs, log10_clamped(slice; floor_value = floor_value))
        for (label, (ys, zs, slice)) in slice_data
    ]
    finite_values = Float64[]
    for (_, _, _, slice) in log_slice_data
        append!(finite_values, filter(isfinite, vec(slice)))
    end
    color_limits = (minimum(finite_values), maximum(finite_values))

    fig = Figure(size = (420 * length(slice_data) + 120, 420))
    heatmap_plot = nothing
    for (i, (label, ys, zs, slice)) in enumerate(log_slice_data)
        ax = Axis(
            fig[1, i],
            title = label,
            xlabel = "y",
            ylabel = i == 1 ? "z" : "",
            aspect = DataAspect(),
        )
        heatmap_plot = heatmap!(ax, ys, zs, slice; colormap = :magma, colorrange = color_limits)
    end

    Colorbar(fig[1, length(slice_data) + 1], heatmap_plot, label = "log10(screened density at x = 0)")
    Label(fig[0, 1:length(slice_data)], "Screened density slice on the plane x = 0", fontsize = 22)

    output = joinpath(@__DIR__, "..", "figs", "screened_density_xslice.png")
    save(output, fig)
    return output
end

function save_kz_decay_figure(vs, bounds, cases)
    spectra = [
        begin
            rho = screened_density_vector(vs, bounds, EPS_IN, EPS_OUT, mode; tol = TOL)
            kz_values, decay = normalized_kz_decay(vs, rho; num_points = 450)
            (label, kz_values, decay)
        end
        for (label, mode) in cases
    ]

    fig = Figure(size = (900, 560))
    ax = Axis(
        fig[1, 1],
        xlabel = L"k_z",
        ylabel = L"|\hat{\rho}(0,0,k_z)| / |\hat{\rho}(0,0,0)|",
        yscale = log10,
        title = L"Decay of \hat{\rho}(0,0,k_z)",
    )

    colors = Makie.wong_colors()
    for (i, (label, kz_values, decay)) in enumerate(spectra)
        lines!(ax, kz_values, decay; label = label, color = colors[mod1(i, length(colors))], linewidth = 3)
    end
    axislegend(ax; position = :rt, framevisible = true)

    output = joinpath(@__DIR__, "..", "figs", "screened_density_kz_decay.png")
    save(output, fig)
    return output
end

function save_refined_qz_decay_figure(vs, bounds, cases)
    spectra = [
        begin
            kz_values, decay = normalized_refined_qz_decay(
                vs,
                vs.density,
                bounds,
                EPS_IN,
                EPS_OUT,
                mode;
                upsample_factor = QZ_UPSAMPLE_FACTOR,
                num_points = 450,
                tol = TOL,
            )
            (label, kz_values, decay)
        end
        for (label, mode) in cases
    ]

    fig = Figure(size = (900, 560))
    ax = Axis(
        fig[1, 1],
        xlabel = L"k_z",
        ylabel = L"|\hat{\rho}(0,0,k_z)| / |\hat{\rho}(0,0,0)|",
        yscale = log10,
        title = "Refined 1D Q(z) decay",
    )

    colors = Makie.wong_colors()
    for (i, (label, kz_values, decay)) in enumerate(spectra)
        lines!(ax, kz_values, decay; label = label, color = colors[mod1(i, length(colors))], linewidth = 3)
    end
    axislegend(ax; position = :rt, framevisible = true)

    output = joinpath(@__DIR__, "..", "figs", "screened_density_kz_decay_refined_qz.png")
    save(output, fig)
    return output
end

function save_global_nufft_qz_decay_figure(vs, bounds, cases)
    spectra = [
        begin
            kz_values, decay = normalized_global_nufft_qz_decay(
                vs,
                vs.density,
                bounds,
                EPS_IN,
                EPS_OUT,
                mode;
                upsample_factor = QZ_UPSAMPLE_FACTOR,
                pad_factor = QZ_GLOBAL_PAD_FACTOR,
                num_points = 450,
                tol = TOL,
            )
            (label, kz_values, decay)
        end
        for (label, mode) in cases
    ]

    fig = Figure(size = (900, 560))
    ax = Axis(
        fig[1, 1],
        xlabel = L"k_z",
        ylabel = L"|\hat{\rho}(0,0,k_z)| / |\hat{\rho}(0,0,0)|",
        yscale = log10,
        title = "Global padded NUFFT Q(z) decay",
    )

    colors = Makie.wong_colors()
    for (i, (label, kz_values, decay)) in enumerate(spectra)
        lines!(ax, kz_values, decay; label = label, color = colors[mod1(i, length(colors))], linewidth = 3)
    end
    axislegend(ax; position = :rt, framevisible = true)

    output = joinpath(@__DIR__, "..", "figs", "screened_density_kz_decay_global_nufft_qz.png")
    save(output, fig)
    return output
end

function main()
    datagrid, vs = default_volume_source(; orbital_file = ORBITAL_FILE, tol = TOL)
    bounds = box_bounds(L, L, LZ)
    cases = analysis_cases()

    slice_output = save_slice_figure(datagrid, bounds, cases)
    decay_output = save_kz_decay_figure(vs, bounds, cases)
    refined_decay_output = save_refined_qz_decay_figure(vs, bounds, cases)
    global_nufft_decay_output = save_global_nufft_qz_decay_figure(vs, bounds, cases)

    println("Saved slice figure to: $slice_output")
    println("Saved kz-decay figure to: $decay_output")
    println("Saved refined qz-decay figure to: $refined_decay_output")
    println("Saved global nufft qz-decay figure to: $global_nufft_decay_output")
end

main()

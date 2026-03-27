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
const BANDWIDTHS = [0.1, 0.5, 1.0]
const QZ_UPSAMPLE_FACTOR = 8
const QZ_GLOBAL_PAD_FACTOR = 4
const Z_UPSAMPLE_FACTORS = [1, 2, 3]
const FULL_SOURCE_TOL = 0.0

function analysis_cases()
    cases = Tuple{String, BI.AbstractScreeningMode}[("Sharp", BI.SharpScreening())]
    for bandwidth in BANDWIDTHS
        push!(cases, ("SoftMix b = $(bandwidth)", BI.SoftMixInversePermittivity(bandwidth)))
    end
    return cases
end

function analysis_inputs(; verbose::Bool = false)
    if verbose
        println("[screened-density] Loading centered volume source")
        flush(stdout)
    end
    datagrid, vs = default_volume_source(; orbital_file = ORBITAL_FILE, tol = TOL)
    bounds = box_bounds(L, L, LZ)
    cases = analysis_cases()
    if verbose
        println("[screened-density] Loaded source and analysis cases")
        flush(stdout)
    end
    return datagrid, vs, bounds, cases
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
            kz_values, decay = normalized_tkm3d_kz_slab_max_decay(vs, rho; kmax = BI._estimate_tkm3dc_kmax(vs), eps = TOL)
            (label, kz_values, decay)
        end
        for (label, mode) in cases
    ]

    positive_values = Float64[]
    for (_, _, decay) in spectra
        append!(positive_values, filter(>(0.0), decay))
    end
    floor_value = isempty(positive_values) ? 1e-16 : minimum(positive_values)
    clamped_spectra = [
        (label, kz_values, max.(decay, floor_value))
        for (label, kz_values, decay) in spectra
    ]

    fig = Figure(size = (900, 560))
    ax = Axis(
        fig[1, 1],
        xlabel = L"k_z",
        ylabel = L"\max_{k_x,k_y} |\hat{\rho}(k_x,k_y,k_z)| \,/\, \max_{k_x,k_y} |\hat{\rho}(k_x,k_y,0)|",
        yscale = log10,
        title = L"Decay of \max_{k_x,k_y} |\hat{\rho}(k_x,k_y,k_z)|",
    )

    colors = Makie.wong_colors()
    for (i, (label, kz_values, decay)) in enumerate(clamped_spectra)
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

function save_z_upsampled_kz_decay_figure(datagrid, bounds, cases; verbose::Bool = false)
    if verbose
        println("[screened-density] Building z-upsampled kz-decay figure")
        flush(stdout)
    end
    spectra = [
        begin
            if verbose
                println("[screened-density] Case: $label")
                flush(stdout)
            end
            curve_data = [
                (factor == 1 ? "Original" : "z x $(factor)", kz_values, decay)
                for (factor, kz_values, decay) in z_upsampled_tkm3d_decay_curves(
                    datagrid,
                    bounds,
                    EPS_IN,
                    EPS_OUT,
                    mode;
                    upsample_factors = Z_UPSAMPLE_FACTORS,
                    source_tol = FULL_SOURCE_TOL,
                    tol = TOL,
                    verbose = verbose,
                    label = label,
                )
            ]
            (label, curve_data)
        end
        for (label, mode) in cases
    ]

    n_panels = length(spectra)
    ncols = 3
    nrows = cld(n_panels, ncols)
    fig = Figure(size = (480 * ncols, 360 * nrows + 80))
    colors = Makie.wong_colors()
    first_axis = nothing

    for (index, (panel_label, curve_data)) in enumerate(spectra)
        row = cld(index, ncols)
        col = mod1(index, ncols)
        ax = Axis(
            fig[row, col],
            title = panel_label,
            xlabel = L"k_z",
            ylabel = col == 1 ? L"\max_{k_x,k_y} |\hat{\rho}(k_x,k_y,k_z)| \,/\, \max_{k_x,k_y} |\hat{\rho}(k_x,k_y,0)|" : "",
            yscale = log10,
        )
        isnothing(first_axis) && (first_axis = ax)

        positive_values = Float64[]
        for (_, _, decay) in curve_data
            append!(positive_values, filter(>(0.0), decay))
        end
        floor_value = isempty(positive_values) ? 1e-16 : minimum(positive_values)

        for (curve_index, (curve_label, kz_values, decay)) in enumerate(curve_data)
            lines!(
                ax,
                kz_values,
                max.(decay, floor_value);
                label = curve_label,
                color = colors[mod1(curve_index, length(colors))],
                linewidth = 3,
            )
        end
    end

    axislegend(first_axis; position = :lb, framevisible = true)
    Label(
        fig[0, 1:ncols],
        "Decay after global z-only NUFFT upsampling of the original density",
        fontsize = 22,
    )

    output = joinpath(@__DIR__, "..", "figs", "screened_density_kz_decay_z_upsampled.png")
    save(output, fig)
    if verbose
        println("[screened-density] Saved z-upsampled kz-decay figure to: $output")
        flush(stdout)
    end
    return output
end

function main_all()
    datagrid, vs, bounds, cases = analysis_inputs(; verbose = true)
    println("[screened-density] Rendering slice figure")
    flush(stdout)
    slice_output = save_slice_figure(datagrid, bounds, cases)
    println("[screened-density] Rendering slab-max kz-decay figure")
    flush(stdout)
    decay_output = save_kz_decay_figure(vs, bounds, cases)
    println("[screened-density] Rendering refined qz-decay figure")
    flush(stdout)
    refined_decay_output = save_refined_qz_decay_figure(vs, bounds, cases)
    println("[screened-density] Rendering global nufft qz-decay figure")
    flush(stdout)
    global_nufft_decay_output = save_global_nufft_qz_decay_figure(vs, bounds, cases)
    z_upsampled_decay_output = save_z_upsampled_kz_decay_figure(datagrid, bounds, cases; verbose = true)

    println("Saved slice figure to: $slice_output")
    println("Saved kz-decay figure to: $decay_output")
    println("Saved refined qz-decay figure to: $refined_decay_output")
    println("Saved global nufft qz-decay figure to: $global_nufft_decay_output")
    println("Saved z-upsampled kz-decay figure to: $z_upsampled_decay_output")
end

function main_z_upsampled()
    datagrid, _, bounds, cases = analysis_inputs(; verbose = true)
    output = save_z_upsampled_kz_decay_figure(datagrid, bounds, cases; verbose = true)
    println("Saved z-upsampled kz-decay figure to: $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_all()
end

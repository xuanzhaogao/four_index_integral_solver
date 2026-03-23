using FFTW
using FINUFFT
using CairoMakie
using HCubature
using SpecialFunctions

const DEFAULT_SCREEN_NYQST_SIGMA = 1.0
const DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS = 8.0
const DEFAULT_SCREEN_NYQST_NPOINTS = 256
const DEFAULT_SCREEN_NYQST_NUFFT_EPS = 1.0e-12
const DEFAULT_ERF_SCREEN_BANDWIDTHS = [0.01, 0.1, 0.5]
const DEFAULT_ERF_SCREEN_FIG_PATH = joinpath(dirname(@__DIR__), "figs", "screened_gauss_erf_decay.svg")
const DEFAULT_SCREEN_LINE_CONV_TARGET_HALF_WIDTH_SIGMAS = 4.0
const DEFAULT_SCREEN_LINE_CONV_TARGET_NPOINTS = 256
const DEFAULT_SCREEN_LINE_CONV_EXCLUSION_SIGMAS = 0.05
const DEFAULT_SCREEN_LINE_CONV_RTOL = 1.0e-6
const DEFAULT_SCREEN_LINE_CONV_ATOL = 1.0e-10
const DEFAULT_SCREEN_LINE_CONV_FIG_PATH = joinpath(dirname(@__DIR__), "figs", "screened_gauss_line_convolution.svg")

function screen_function(x::Float64)
    if x < 0.0
        return 0.0
    elseif x > 0.0
        return 1.0
    end
    return 0.5
end

gaussian_screen_base(x::Float64, sigma::Float64) = exp(-0.5 * (x / sigma)^2) / sqrt(2.0 * pi * sigma^2)
gaussian_screen_base_ft(k::Float64, sigma::Float64) = exp(-0.5 * (sigma * k)^2)

screened_gaussian_1d(x::Float64, sigma::Float64) = screen_function(x) * gaussian_screen_base(x, sigma)
erf_screen_function(x::Float64, bandwidth::Float64) = 0.5 * (1.0 + erf(x / (sqrt(2.0) * bandwidth)))
erf_screened_gaussian_1d(x::Float64, sigma::Float64, bandwidth::Float64) =
    erf_screen_function(x, bandwidth) * gaussian_screen_base(x, sigma)

function screened_gaussian_1d_ft(k::Float64, sigma::Float64)
    t = sigma * k / sqrt(2.0)
    # This avoids the overflow in exp(-t^2) * erfc(im * t) for large |t|.
    return 0.5 * exp(-t^2) - im * dawson(t) / sqrt(pi)
end

function screen_centered_spatial_axis(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
    npoints::Int = DEFAULT_SCREEN_NYQST_NPOINTS,
)
    if !iseven(npoints)
        throw(ArgumentError("npoints must be even"))
    end

    half_width = half_width_sigmas * sigma
    dx = 2.0 * half_width / npoints
    axis = dx .* collect((-npoints ÷ 2):(npoints ÷ 2 - 1))
    return axis, dx
end

function screen_truncated_interval(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
)
    half_width = half_width_sigmas * sigma
    return -half_width, half_width
end

function screen_centered_frequency_axis(npoints::Int, dx::Float64)
    dk = 2.0 * pi / (npoints * dx)
    axis = dk .* collect((-npoints ÷ 2):(npoints ÷ 2 - 1))
    return axis, dk
end

function screen_refined_k_axis(k::AbstractVector{<:Real}; refinement::Int = 8)
    if refinement < 1
        throw(ArgumentError("refinement must be at least 1"))
    end
    return collect(range(first(k), last(k); length = refinement * (length(k) - 1) + 1))
end

screen_extended_mode_axis(mode_extent::Int) = collect(-mode_extent:mode_extent)

screen_finufft_nodes_from_x(x::AbstractVector{<:Real}, dk::Float64) = dk .* collect(Float64.(x))

screen_continuous_fft(values::AbstractVector{<:Real}, dx::Float64) =
    dx .* fftshift(fft(ifftshift(collect(values))))

function screen_continuous_nufft_type1(
    x::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    dx::Float64;
    mode_extent::Int = length(x),
    eps::Float64 = DEFAULT_SCREEN_NYQST_NUFFT_EPS,
)
    npoints = length(x)
    dk = 2.0 * pi / (npoints * dx)
    theta = screen_finufft_nodes_from_x(x, dk)
    weights = ComplexF64.(dx .* values)
    modes = screen_extended_mode_axis(mode_extent)
    fk = vec(nufft1d1(theta, weights, -1, eps, length(modes)))
    k = dk .* modes
    return k, fk
end

function _validated_bandwidths(bandwidths::AbstractVector{<:Real})
    isempty(bandwidths) && throw(ArgumentError("bandwidths must be nonempty"))
    values = Float64.(bandwidths)
    all(>(0.0), values) || throw(ArgumentError("bandwidths must be positive"))
    return values
end

function _validated_target_npoints(target_npoints::Int)
    target_npoints >= 2 || throw(ArgumentError("target_npoints must be at least 2"))
    return target_npoints
end

function _validated_interior_breakpoints(
    xmin::Real,
    xmax::Real,
    breakpoints::AbstractVector{<:Real},
)
    values = Float64[]
    for breakpoint in breakpoints
        xmin < breakpoint < xmax || continue
        push!(values, Float64(breakpoint))
    end
    sort!(unique!(values))
    return values
end

function screen_convolution_target_axis(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    half_width_sigmas::Float64 = DEFAULT_SCREEN_LINE_CONV_TARGET_HALF_WIDTH_SIGMAS,
    npoints::Int = DEFAULT_SCREEN_LINE_CONV_TARGET_NPOINTS,
)
    validated_npoints = _validated_target_npoints(npoints)
    half_width = half_width_sigmas * sigma
    return collect(range(-half_width, half_width; length = validated_npoints))
end

function _positive_frequency_magnitude_curve(
    k::AbstractVector{<:Real},
    fhat::AbstractVector{<:Complex},
)
    positive = findall(>=(0.0), k)
    return (
        k_positive = Float64.(k[positive]),
        magnitude_positive = max.(abs.(fhat[positive]), eps(Float64)),
    )
end

function _named_decay_curve(
    label::AbstractString,
    fx::AbstractVector{<:Real},
    fhat::AbstractVector{<:Complex},
    k::AbstractVector{<:Real},
)
    positive_curve = _positive_frequency_magnitude_curve(k, fhat)
    return (
        label = String(label),
        fx = Float64.(fx),
        fhat = ComplexF64.(fhat),
        k_positive = positive_curve.k_positive,
        magnitude_positive = positive_curve.magnitude_positive,
    )
end

function regularized_line_convolution(
    source::Function,
    x0::Real;
    xmin::Real,
    xmax::Real,
    delta::Real,
    breakpoints::AbstractVector{<:Real} = Float64[],
    rtol::Real = DEFAULT_SCREEN_LINE_CONV_RTOL,
    atol::Real = DEFAULT_SCREEN_LINE_CONV_ATOL,
)
    xmin < xmax || throw(ArgumentError("xmin must be less than xmax"))
    delta > 0 || throw(ArgumentError("delta must be positive"))
    rtol >= 0 || throw(ArgumentError("rtol must be nonnegative"))
    atol >= 0 || throw(ArgumentError("atol must be nonnegative"))

    x0_value = Float64(x0)
    delta_value = Float64(delta)
    edges = Float64[Float64(xmin), Float64(xmax)]
    append!(edges, _validated_interior_breakpoints(xmin, xmax, breakpoints))

    left_exclusion = x0_value - delta_value
    right_exclusion = x0_value + delta_value
    if xmin < left_exclusion < xmax
        push!(edges, left_exclusion)
    end
    if xmin < right_exclusion < xmax
        push!(edges, right_exclusion)
    end
    sort!(unique!(edges))

    value = 0.0
    error = 0.0
    for i in 1:(length(edges) - 1)
        a = edges[i]
        b = edges[i + 1]
        a < b || continue

        midpoint = 0.5 * (a + b)
        if abs(midpoint - x0_value) < delta_value
            continue
        end

        segment_value, segment_error = hcubature(
            y -> source(y[1]) / abs(y[1] - x0_value),
            [a],
            [b];
            rtol = Float64(rtol),
            atol = Float64(atol),
        )
        value += segment_value
        error += segment_error
    end

    return (value = value, error = error)
end

function _named_convolution_curve(
    label::AbstractString,
    source::Function,
    targets::AbstractVector{<:Real};
    xmin::Real,
    xmax::Real,
    delta::Real,
    breakpoints::AbstractVector{<:Real} = Float64[],
    rtol::Real = DEFAULT_SCREEN_LINE_CONV_RTOL,
    atol::Real = DEFAULT_SCREEN_LINE_CONV_ATOL,
)
    values = Vector{Float64}(undef, length(targets))
    errors = Vector{Float64}(undef, length(targets))
    for (idx, target) in enumerate(targets)
        result = regularized_line_convolution(
            source,
            target;
            xmin,
            xmax,
            delta,
            breakpoints = breakpoints,
            rtol = rtol,
            atol = atol,
        )
        values[idx] = result.value
        errors[idx] = result.error
    end

    return (
        label = String(label),
        targets = Float64.(targets),
        potential = values,
        errors = errors,
    )
end

function run_gauss_screen_nyqst_study(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
    npoints::Int = DEFAULT_SCREEN_NYQST_NPOINTS,
    nufft_eps::Float64 = DEFAULT_SCREEN_NYQST_NUFFT_EPS,
)
    x, dx = screen_centered_spatial_axis(; sigma, half_width_sigmas, npoints)
    fx = screened_gaussian_1d.(x, Ref(sigma))
    k, dk = screen_centered_frequency_axis(npoints, dx)
    fhat_numerical = screen_continuous_fft(fx, dx)
    fhat_exact = screened_gaussian_1d_ft.(k, Ref(sigma))
    k_nufft, fhat_nufft = screen_continuous_nufft_type1(x, fx, dx; mode_extent = npoints, eps = nufft_eps)
    fhat_exact_nufft = screened_gaussian_1d_ft.(k_nufft, Ref(sigma))

    return (
        x = x,
        dx = dx,
        fx = fx,
        k = k,
        dk = dk,
        fhat_numerical = fhat_numerical,
        fhat_exact = fhat_exact,
        k_nufft = k_nufft,
        fhat_nufft = fhat_nufft,
        fhat_exact_nufft = fhat_exact_nufft,
    )
end

function run_erf_screen_decay_study(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    bandwidths::AbstractVector{<:Real} = DEFAULT_ERF_SCREEN_BANDWIDTHS,
    half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
    npoints::Int = DEFAULT_SCREEN_NYQST_NPOINTS,
)
    validated_bandwidths = _validated_bandwidths(bandwidths)
    x, dx = screen_centered_spatial_axis(; sigma, half_width_sigmas, npoints)
    k, dk = screen_centered_frequency_axis(npoints, dx)
    nyquist_limit = pi / dx
    positive_k = Float64.(k[k .>= 0.0])
    unscreened_fx = gaussian_screen_base.(x, Ref(sigma))
    unscreened_fhat = ComplexF64.(gaussian_screen_base_ft.(k, Ref(sigma)))
    unscreened_curve = _named_decay_curve("unscreened", unscreened_fx, unscreened_fhat, k)
    step_fx = screened_gaussian_1d.(x, Ref(sigma))
    step_fhat = screen_continuous_fft(step_fx, dx)
    step_curve = _named_decay_curve("step", step_fx, step_fhat, k)

    curves = map(validated_bandwidths) do bandwidth
        fx = erf_screened_gaussian_1d.(x, Ref(sigma), Ref(bandwidth))
        fhat = screen_continuous_fft(fx, dx)
        curve = _named_decay_curve("bandwidth = $(bandwidth)", fx, fhat, k)
        return merge((bandwidth = bandwidth,), curve)
    end

    return (
        x = x,
        dx = dx,
        k = k,
        dk = dk,
        nyquist_limit = nyquist_limit,
        positive_k = positive_k,
        bandwidths = validated_bandwidths,
        unscreened_curve = unscreened_curve,
        step_curve = step_curve,
        curves = curves,
    )
end

function run_screened_line_convolution_study(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    bandwidths::AbstractVector{<:Real} = DEFAULT_ERF_SCREEN_BANDWIDTHS,
    source_half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
    target_half_width_sigmas::Float64 = DEFAULT_SCREEN_LINE_CONV_TARGET_HALF_WIDTH_SIGMAS,
    target_npoints::Int = DEFAULT_SCREEN_LINE_CONV_TARGET_NPOINTS,
    exclusion_sigmas::Float64 = DEFAULT_SCREEN_LINE_CONV_EXCLUSION_SIGMAS,
    rtol::Float64 = DEFAULT_SCREEN_LINE_CONV_RTOL,
    atol::Float64 = DEFAULT_SCREEN_LINE_CONV_ATOL,
)
    validated_bandwidths = _validated_bandwidths(bandwidths)
    exclusion_sigmas > 0 || throw(ArgumentError("exclusion_sigmas must be positive"))
    xmin, xmax = screen_truncated_interval(; sigma, half_width_sigmas = source_half_width_sigmas)
    targets = screen_convolution_target_axis(; sigma, half_width_sigmas = target_half_width_sigmas, npoints = target_npoints)
    delta = exclusion_sigmas * sigma

    step_curve = _named_convolution_curve(
        "step",
        x -> screened_gaussian_1d(x, sigma),
        targets;
        xmin,
        xmax,
        delta,
        breakpoints = [0.0],
        rtol = rtol,
        atol = atol,
    )

    curves = map(validated_bandwidths) do bandwidth
        curve = _named_convolution_curve(
            "bandwidth = $(bandwidth)",
            x -> erf_screened_gaussian_1d(x, sigma, bandwidth),
            targets;
            xmin,
            xmax,
            delta,
            rtol = rtol,
            atol = atol,
        )
        return merge(
            (bandwidth = bandwidth, difference_from_step = curve.potential .- step_curve.potential),
            curve,
        )
    end

    return (
        sigma = sigma,
        xmin = xmin,
        xmax = xmax,
        delta = delta,
        targets = Float64.(targets),
        bandwidths = validated_bandwidths,
        source_half_width_sigmas = source_half_width_sigmas,
        target_half_width_sigmas = target_half_width_sigmas,
        step_curve = step_curve,
        curves = curves,
        rtol = rtol,
        atol = atol,
    )
end

function build_erf_screen_decay_figure(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    bandwidths::AbstractVector{<:Real} = DEFAULT_ERF_SCREEN_BANDWIDTHS,
    half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
    npoints::Int = DEFAULT_SCREEN_NYQST_NPOINTS,
)
    study = run_erf_screen_decay_study(; sigma, bandwidths, half_width_sigmas, npoints)

    fig = Figure(size = (900, 500), fontsize = 18)
    ax = Axis(
        fig[1, 1];
        xlabel = "k / (π / dx)",
        ylabel = "|fhat(k)|",
        yscale = log10,
        title = "Gaussian decay with screening",
    )

    lines!(
        ax,
        study.unscreened_curve.k_positive ./ study.nyquist_limit,
        study.unscreened_curve.magnitude_positive;
        color = :gray35,
        linewidth = 3,
        label = study.unscreened_curve.label,
    )

    lines!(
        ax,
        study.step_curve.k_positive ./ study.nyquist_limit,
        study.step_curve.magnitude_positive;
        color = :black,
        linestyle = :dash,
        linewidth = 3,
        label = study.step_curve.label,
    )

    palette = Makie.wong_colors()
    for (idx, curve) in enumerate(study.curves)
        lines!(
            ax,
            curve.k_positive ./ study.nyquist_limit,
            curve.magnitude_positive;
            color = palette[mod1(idx, length(palette))],
            linewidth = 3,
            label = curve.label,
        )
    end

    xlims!(ax, 0.0, 1.0)
    # axislegend(ax; position = :rt)

    Legend(fig[1, 2], ax, orientation = :vertical)

    return fig, study
end

function build_screened_line_convolution_figure(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    bandwidths::AbstractVector{<:Real} = DEFAULT_ERF_SCREEN_BANDWIDTHS,
    source_half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
    target_half_width_sigmas::Float64 = DEFAULT_SCREEN_LINE_CONV_TARGET_HALF_WIDTH_SIGMAS,
    target_npoints::Int = DEFAULT_SCREEN_LINE_CONV_TARGET_NPOINTS,
    exclusion_sigmas::Float64 = DEFAULT_SCREEN_LINE_CONV_EXCLUSION_SIGMAS,
    rtol::Float64 = DEFAULT_SCREEN_LINE_CONV_RTOL,
    atol::Float64 = DEFAULT_SCREEN_LINE_CONV_ATOL,
)
    study = run_screened_line_convolution_study(;
        sigma,
        bandwidths,
        source_half_width_sigmas,
        target_half_width_sigmas,
        target_npoints,
        exclusion_sigmas,
        rtol,
        atol,
    )

    normalized_targets = study.targets ./ sigma
    fig = Figure(size = (900, 800), fontsize = 18)
    ax_potential = Axis(
        fig[1, 1];
        xlabel = "x / σ",
        ylabel = "u_δ(x)",
        title = "Regularized same-line 1/r convolution",
    )
    ax_difference = Axis(
        fig[2, 1];
        xlabel = "x / σ",
        ylabel = "|u_δ^erf(x) - u_δ^step(x)",
        yscale = log10,
    )

    lines!(
        ax_potential,
        normalized_targets,
        study.step_curve.potential;
        color = :black,
        linestyle = :dash,
        linewidth = 3,
        label = study.step_curve.label,
    )
    hlines!(ax_difference, [0.0]; color = :gray55, linestyle = :dot, linewidth = 2)
    vlines!(ax_potential, [0.0]; color = :gray70, linestyle = :dot, linewidth = 1.5)
    vlines!(ax_difference, [0.0]; color = :gray70, linestyle = :dot, linewidth = 1.5)

    palette = Makie.wong_colors()
    for (idx, curve) in enumerate(study.curves)
        color = palette[mod1(idx, length(palette))]
        lines!(
            ax_potential,
            normalized_targets,
            curve.potential;
            color = color,
            linewidth = 3,
            label = curve.label,
        )
        lines!(
            ax_difference,
            normalized_targets,
            abs.(curve.difference_from_step);
            color = color,
            linewidth = 3,
        )
    end

    axislegend(ax_potential; position = :rb)
    return fig, study
end

function main(;
    sigma::Float64 = DEFAULT_SCREEN_NYQST_SIGMA,
    bandwidths::AbstractVector{<:Real} = DEFAULT_ERF_SCREEN_BANDWIDTHS,
    half_width_sigmas::Float64 = DEFAULT_SCREEN_NYQST_HALF_WIDTH_SIGMAS,
    npoints::Int = DEFAULT_SCREEN_NYQST_NPOINTS,
    decay_output_path::AbstractString = DEFAULT_ERF_SCREEN_FIG_PATH,
    convolution_output_path::AbstractString = DEFAULT_SCREEN_LINE_CONV_FIG_PATH,
    convolution_target_half_width_sigmas::Float64 = DEFAULT_SCREEN_LINE_CONV_TARGET_HALF_WIDTH_SIGMAS,
    convolution_target_npoints::Int = DEFAULT_SCREEN_LINE_CONV_TARGET_NPOINTS,
    convolution_exclusion_sigmas::Float64 = DEFAULT_SCREEN_LINE_CONV_EXCLUSION_SIGMAS,
    convolution_rtol::Float64 = DEFAULT_SCREEN_LINE_CONV_RTOL,
    convolution_atol::Float64 = DEFAULT_SCREEN_LINE_CONV_ATOL,
)
    decay_fig, decay_study = build_erf_screen_decay_figure(; sigma, bandwidths, half_width_sigmas, npoints)
    convolution_fig, convolution_study = build_screened_line_convolution_figure(;
        sigma,
        bandwidths,
        source_half_width_sigmas = half_width_sigmas,
        target_half_width_sigmas = convolution_target_half_width_sigmas,
        target_npoints = convolution_target_npoints,
        exclusion_sigmas = convolution_exclusion_sigmas,
        rtol = convolution_rtol,
        atol = convolution_atol,
    )

    mkpath(dirname(decay_output_path))
    mkpath(dirname(convolution_output_path))
    save(decay_output_path, decay_fig)
    save(convolution_output_path, convolution_fig)

    println("saved erf-screened Gaussian decay figure to $(decay_output_path)")
    println("saved screened Gaussian line-convolution figure to $(convolution_output_path)")
    println("bandwidths = $(collect(decay_study.bandwidths))")
    println("convolution exclusion radius = $(convolution_study.delta)")
    return (
        decay = (fig = decay_fig, study = decay_study),
        convolution = (fig = convolution_fig, study = convolution_study),
    )
end

main()

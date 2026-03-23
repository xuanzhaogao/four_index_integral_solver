using FFTW
using FINUFFT
using CairoMakie, LaTeXStrings
using HCubature
using SpecialFunctions

function step_function(x)
    return x < 0.5 ? 0.5 : 1.0
end

function soft_mix(x, bandwidth)
    return 0.25 * (erf((x - 0.5) / bandwidth)) + 0.75
end

function gauss(x, sigma)
    return exp(-0.5 * (x / sigma)^2) / (sigma * sqrt(2π))
end

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax_1 = Axis(fig[1, 1]; xlabel = L"x / \sigma", ylabel = L"\rho(x)")
    ax_2 = Axis(fig[1, 2]; xlabel = L"k \sigma", ylabel = L"\hat{\rho}(k)", yscale = log10)
    x = range(-5.0, 5.0; length = 512)
    N = length(x)
    dx = step(x)
    k = 2π .* (-(N ÷ 2):((N - 1) ÷ 2)) ./ (N * dx)

    bandwidths = [0.05, 0.1, 0.2]
    u_step = step_function.(x) .* gauss.(x, 1.0)
    lines!(ax_1, x, u_step; label = "step")
    u_step_ft = abs.(fftshift(fft(u_step))) .* dx
    lines!(ax_2, k[k .> 0], u_step_ft[k .> 0]; label = "step FT")

    for bandwidth in bandwidths
        u_erf = soft_mix.(x, bandwidth) .* gauss.(x, 1.0)
        lines!(ax_1, x, u_erf; label = "erf, bw = $bandwidth")
        u_erf_ft = abs.(fftshift(fft(u_erf))) .* dx
        lines!(ax_2, k[k .> 0], u_erf_ft[k .> 0]; label = "erf FT, bw = $bandwidth")
    end

    xlims!(ax_1, -5.0, 5.0)
    xlims!(ax_2, 0.0, maximum(k))

    axislegend(ax_1; position = :lt)

    fig

end

save(joinpath(dirname(@__DIR__), "figs", "screened_gauss_ft.svg"), fig)
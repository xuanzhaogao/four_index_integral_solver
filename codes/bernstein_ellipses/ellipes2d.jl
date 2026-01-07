using LinearAlgebra
using FastGaussQuadrature
using CairoMakie

f = x -> sin(x)

function integrate(f, xs, ws, xt, yt)
    res = 0.0
    for (x, w) in zip(xs, ws)
        r = sqrt((x - xt)^2 + (yt)^2)
        res += f(x) * w / r
    end
    return res
end

function bernstein_rho_from_pole(z::Complex{T}) where T <: Real
    s  = sqrt(z*z - 1)          # complex sqrt (principal branch)
    w1 = z + s
    w2 = z - s
    return max(abs(w1), abs(w2))
end

xa, xb = -2.0, 2.0
ya, yb = -1.5, 1.5

xt_s = range(xa, xb, length = 400)
yt_s = range(ya, yb, length = 300)

xs_ref, ws_ref = gausslegendre(128)
res_ref = zeros(length(xt_s), length(yt_s))
for (i, xt) in enumerate(xt_s)
    for (j, yt) in enumerate(yt_s)
        res_ref[i, j] = integrate(f, xs_ref, ws_ref, xt, yt)
    end
end

# check if n = 16, how error changes

n = 16
xs, ys = gausslegendre(n)
res = zeros(length(xt_s), length(yt_s))
for (i, xt) in enumerate(xt_s)
    for (j, yt) in enumerate(yt_s)
        res[i, j] = integrate(f, xs, ys, xt, yt)
    end
end

error_predicted = zeros(length(xt_s), length(yt_s))
for (i, xt) in enumerate(xt_s)
    for (j, yt) in enumerate(yt_s)
        rho = bernstein_rho_from_pole(ComplexF64(xt, yt))
        error_predicted[i, j] = rho ^ (- 2 * n)
    end
end

begin
    fig_2 = Figure()
    ax_2 = Axis(fig_2[1, 1], aspect=DataAspect(), title = "f = sin(x), n = 16")
    hm_2 = contourf!(ax_2, xt_s, yt_s, log10.(abs.(res .- res_ref)), levels = -19.0:2.0:1.0)
    Colorbar(fig_2[1, 2], hm_2)
    hm_3 = contour!(ax_2, xt_s, yt_s, log10.(abs.(error_predicted)), color = :red, levels = -19.0:2.0:1.0)

    xlims!(ax_2, xa, xb)
    ylims!(ax_2, ya, yb)

    save("ellipes2d_n16.png", fig_2, dpi = 300)

    fig_2
end

# then I want to see how error decay as I increase n for a fixed xt, yt
xt, yt = 0.1, 0.2

xs_ref, ws_ref = gausslegendre(256)
f_ref = integrate(f, xs_ref, ws_ref, xt, yt)

n_s = 4:4:128
res = zeros(length(n_s))
for (i, n) in enumerate(n_s)
    xs, ys = gausslegendre(n)
    res[i] = integrate(f, xs, ys, xt, yt)
end

# plot the error decay
rho = bernstein_rho_from_pole(ComplexF64(xt, yt))


begin
    fig_3 = Figure(size = (500, 400), fontsize = 15)
    ax_3 = Axis(fig_3[1, 1], yscale = log10, xlabel = "n", ylabel = "error", title = "f = sin(x), xt = $(xt), yt = $(yt)")
    scatter!(ax_3, n_s, abs.(res .- f_ref))
    lines!(ax_3, n_s, rho .^ (- 2 * n_s))
    ylims!(ax_3, 1e-18, 1e-1)
    fig_3
end

# check the consistency

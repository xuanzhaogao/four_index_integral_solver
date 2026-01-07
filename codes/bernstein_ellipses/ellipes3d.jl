using LinearAlgebra
using FastGaussQuadrature
using CairoMakie

# f = (x, y) -> sin(x) * cos(2.0 * y)
f = (x, y) -> exp(-x^2 - y^2)

function integrate2d(f, xs, ys, ws_x, ws_y, xt, yt, zt)
    res = 0.0
    for (x, w_x) in zip(xs, ws_x)
        for (y, w_y) in zip(ys, ws_y)
            r2 = (x - xt)^2 + (y - yt)^2 + (zt)^2
            res += f(x, y) * w_x * w_y / r2
        end
    end
    return res
end

function bernstein_rho_from_pole(z::Complex{T}) where T <: Real
    s  = sqrt(z*z - 1)          # complex sqrt (principal branch)
    w1 = z + s
    w2 = z - s
    return max(abs(w1), abs(w2))
end

# check the error on y = 0 plane first
xa, xb = -2.0, 2.0
za, zb = -1.5, 1.5

xt_s = range(xa, xb, length = 200)
zt_s = range(za, zb, length = 150)

# for yt in [2.8]
for yt in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5]

begin
    n_ref = 128
    xs_ref, ws_x = gausslegendre(n_ref)
    ys_ref, ws_y = gausslegendre(n_ref)
    res_ref = zeros(length(xt_s), length(zt_s))
    for i in 1:length(xt_s)
        xt = xt_s[i]
        for j in 1:length(zt_s)
            zt = zt_s[j]
            res_ref[i, j] = integrate2d(f, xs_ref, ys_ref, ws_x, ws_y, xt, yt, zt)
        end
    end

    n = 16
    xs, ws_x2 = gausslegendre(n)
    ys, ws_y2 = gausslegendre(n)
    res = zeros(length(xt_s), length(zt_s))
    for (i, xt) in enumerate(xt_s)
        for (j, zt) in enumerate(zt_s)
            res[i, j] = integrate2d(f, xs, ys, ws_x2, ws_y2, xt, yt, zt)
        end
    end

    error_predicted = zeros(length(xt_s), length(zt_s))
    for (i, xt) in enumerate(xt_s)
        for (j, zt) in enumerate(zt_s)
            # a rough estimatation of rho_x
            # rho_min = Inf
            # for yy in -1.0:0.01:1.0
            #     rho_1 = bernstein_rho_from_pole(ComplexF64(xt, sqrt(yy^2 + zt^2)))
            #     rho_min = min(rho_min, rho_1)
            # end
            rho_x_min = Inf
            rho_y_min = Inf
            for t in -1.0:0.1:1.0
                rho_x = bernstein_rho_from_pole(ComplexF64(xt, sqrt((yt - t)^2 + zt^2)))
                rho_y = bernstein_rho_from_pole(ComplexF64(yt, sqrt((xt - t)^2 + zt^2)))
                rho_x_min = min(rho_x_min, rho_x)
                rho_y_min = min(rho_y_min, rho_y)
            end
            error_predicted[i, j] = max(rho_x_min ^ (- 2 * n), rho_y_min ^ (- 2 * n))
        end
    end
end

begin
    levels_log = -21.0:2.0:1.0
    log_true = log10.(abs.(res .- res_ref) .+ eps())
    log_pred = log10.(abs.(error_predicted) .+ eps())

    fig_3d = Figure(size = (600, 420))

    ax_err = Axis(
        fig_3d[1, 1],
        aspect = DataAspect(),
        xlabel = "x",
        ylabel = "z",
        title = "error (log10|res-res_ref|), n = 16, yt = $(yt) + predicted contours",
    )
    hm_err = contourf!(ax_err, xt_s, zt_s, log_true; levels = levels_log, colormap = :jet)

    # Make each predicted-error contour line use the SAME color mapping as the error plot
    contour!(
        ax_err,
        xt_s,
        zt_s,
        log_pred;
        levels = levels_log,
        # colormap = :jet,
        color = :black,
        colorrange = (first(levels_log), last(levels_log)),
        linewidth = 2.0,
        # linestyle = :dash,
    )
    Colorbar(fig_3d[1, 2], hm_err)

    xlims!(ax_err, xa, xb)
    ylims!(ax_err, za, zb)

    save("3d_figs/ellipes3d_overlay_n16_yt$(yt).png", fig_3d, dpi = 300)
end

fig_3d
end
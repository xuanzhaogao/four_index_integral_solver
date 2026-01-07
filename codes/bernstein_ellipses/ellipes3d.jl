using LinearAlgebra
using FastGaussQuadrature
using CairoMakie

f = (x, y) -> sin(x) * cos(2.0 * y)

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

# check the error on y = 0 plane first
xa, xb = -2.0, 2.0
za, zb = -1.5, 1.5

xt_s = range(xa, xb, length = 200)
zt_s = range(za, zb, length = 150)

yt = 0.0

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

begin
    fig_3d = Figure()
    ax_3d = Axis(fig_3d[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "z", title = "f = sin(x) * cos(2y), n = 16, yt = 0.0")
    hm_3d = contourf!(ax_3d, xt_s, zt_s, log10.(abs.(res .- res_ref)))
    Colorbar(fig_3d[1, 2], hm_3d)

    save("ellipes3d_n16.png", fig_3d, dpi = 300)

    fig_3d
end
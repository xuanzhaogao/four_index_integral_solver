using BoundaryIntegral
using BoundaryIntegral.FastGaussQuadrature
import BoundaryIntegral as BI
using Roots
using CairoMakie

eps_m = 1000.0
ang = pi / 2

root_even(n) = fzero(g -> BI.theta_shooting_even(ang, eps_m, g), n)
root_odd(n) = fzero(g -> BI.theta_shooting_odd(ang, eps_m, g), n)

d_even = [root_even(n) for n in 1:5]
d_odd = [root_odd(n) for n in 1:5]

ds = sort(vcat(d_even..., d_odd...)) .- 1.0
f = x -> sum([(cos((di + 1) * x) + sin((di + 1) * x)) * x^(di) for di in ds])

begin
    fig = Figure()
    ax = Axis(fig[1, 1])
    xs = 0.001:0.001:0.1
    lines!(ax, xs, f.(xs))

    lines!(ax, xs, xs .^ (-ds[1]) .* f.(xs))

    fig
end
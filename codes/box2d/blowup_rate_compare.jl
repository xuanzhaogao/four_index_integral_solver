include("utils.jl")
using Roots

function root_even(al, e, n)
    return fzero(g -> BI.theta_shooting_even(al, e, g), n)
end

function root_odd(al, e, n)
    return fzero(g -> BI.theta_shooting_odd(al, e, g), n)
end

gg = 0:0.01:3.5    # also used below

# higher powers w/o rootfinding...
angs = range(0, 2pi, length=200)
eps = [10.0, 1.0]     # diel wedge (ang is for diel)
#eps = [1.0, Inf]   # PEC corner (ang is for vacuum)
dd = [BI.theta_ODE_det([a], eps, g) for a in angs, g in gg]

d_even_1 = [root_even(a, eps[1], 1.0) for a in angs]
d_odd_1 = [root_odd(a, eps[1], 1.0) for a in angs]
d_even_2 = [root_even(a, eps[1], 2.0) for a in angs]
d_odd_2 = [root_odd(a, eps[1], 2.0) for a in angs]
d_even_3 = [root_even(a, eps[1], 3.0) for a in angs]
d_odd_3 = [root_odd(a, eps[1], 3.0) for a in angs]

begin
    fig, ax, p = heatmap(angs, gg, -log.(abs.(dd) .+ 1e-10))    # log highlights the zeros
    p.colormap = :jet;
    p.colorrange = (-10, 10);
    ax.xticks = ([0, pi / 2, pi, 3pi / 2, 2pi], ["0", "π/2", "π", "3π/2", "2π"])

    lines!(ax, angs, d_even_1, color = :red, label = "even parity")
    lines!(ax, angs, d_odd_1, color = :blue, label = "odd parity")
    lines!(ax, angs, d_even_2, color = :red, label = "even parity")
    lines!(ax, angs, d_odd_2, color = :blue, label = "odd parity")
    lines!(ax, angs, d_even_3, color = :red, label = "even parity")
    lines!(ax, angs, d_odd_3, color = :blue, label = "odd parity")
end

save("figs/blowup_rate_compare.png", fig, px_per_unit = 2)
fig
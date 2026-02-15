using LinearAlgebra, SpecialFunctions, ForwardDiff
using CairoMakie

f = x -> erf(x) / (x)
grad_f = x -> ForwardDiff.derivative(f, x)

z_0 = 0.01

xs = range(0.0, 10.0, length = 1000)

ys = [grad_f(sqrt(z_0^2 + x^2)) * z_0 / sqrt(x^2 + z_0^2) for x in xs]

lines(xs, ys)
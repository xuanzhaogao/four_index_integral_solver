include(joinpath(@__DIR__, "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using Printf

res = solve_system(system1(); eps = 1e-3, p = 4, r = 2)
d = 0.01
x0 = [0.25, 0.15, 0.5]
X = Matrix{Float64}(undef, 3, 8)
for (j, t) in enumerate((d, 2d, 3d, 4d, -d, -2d, -3d, -4d))
    X[:, j] .= x0 .+ [0.0, 0.0, t]
end
u_inc = Harness.eval_incident(res, X)
u_sc = Harness.eval_scatter(res, X)
for j in 1:8
    @printf("z=%+.3f  u_inc=%+.6e  u_sc=%+.6e\n", X[3, j] - 0.5, u_inc[j], u_sc[j])
end

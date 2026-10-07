# Physics validation: (a) System I near-conductor limit vs image-charge answer;
# (b) raw transmission derivatives across the System I top face (clear truth).

include(joinpath(@__DIR__, "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, Printf

# (a) conductor limit: cube eps1 = 1e6; image-charge estimate for half-space
sys = system1(eps1 = 1e6)
res = solve_system(sys; eps = 1e-3, p = 4, r = 2)
V = eval_V(res)
r12 = 0.2
rimg = norm([0.2, 0.0, 0.6] .- [0.0, 0.0, 0.4])
V_img = (1 / (4pi)) * (1 / r12 - 1 / rimg)
V_img_flip = (1 / (4pi)) * (1 / r12 + 1 / rimg)
@printf("conductor-limit V = %.6f ; image estimate %.6f (sign-flip would be %.6f)\n", V, V_img, V_img_flip)

# (b) System I (eps1=10) transmission across top face at (0.25, 0.15, 0.5)
res10 = solve_system(system1(); eps = 1e-3, p = 4, r = 2)
d = 0.01
x0 = [0.25, 0.15, 0.5]
X = Matrix{Float64}(undef, 3, 8)
for (j, t) in enumerate((d, 2d, 3d, 4d, -d, -2d, -3d, -4d))
    X[:, j] .= x0 .+ [0.0, 0.0, t]
end
phi = eval_phi(res10, X)
deriv0(f1, f2, f3, dd) = (-2.5f1 + 4f2 - 1.5f3) / dd
dn_p = deriv0(phi[1], phi[2], phi[3], d)
dn_m = -deriv0(phi[5], phi[6], phi[7], d)
# wider-stencil variants to gauge FD/eval noise
dn_p2 = deriv0(phi[2], phi[3], phi[4], d) # nodes 2d,3d,4d -> derivative at d, crude
@printf("top face: dn+ = %+.6f  dn- = %+.6f   eps-*dn- = %+.6f vs eps+*dn+ = %+.6f\n",
        dn_p, dn_m, 10.0 * dn_m, 1.0 * dn_p)
@printf("phi at +/-d: %.6f %.6f ; jump (extrap) = %.2e\n", phi[1], phi[5],
        abs((3phi[1] - 3phi[2] + phi[3]) - (3phi[5] - 3phi[6] + phi[7])))

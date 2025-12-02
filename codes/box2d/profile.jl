using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov, IterativeSolvers
using LinearAlgebra, Random
using BenchmarkTools
using Profile, ProfileSVG

rects = [BI.square(-1.25, -0.25), BI.square(-0.25, 0.25), BI.square(0.25, -0.75), BI.square(-0.75, -1.25)]
eps_boxes = [2.0, 3.0, 4.0, 5.0]
r_src = (0.1, 0.2)
eps_src = 1.0

n_panel = 8
n_quad = 16
n_adapt = 20

mbox = BI.dielectric_mbox2d(eps_boxes, rects, n_panel, n_quad, n_adapt)
n_points = BI.num_points(mbox)

lhs_fmm2d = BI.Lhs_dielectric_mbox2d_fmm2d(mbox, 1e-4)
lhs_direct = BI.Lhs_dielectric_mbox2d(mbox)

rhs = BI.Rhs_dielectric_mbox2d(mbox, r_src, eps_src)

@btime $(lhs_fmm2d) * $(ones(n_points))
@btime $(lhs_direct) * $(ones(n_points))

@time Krylov.gmres(lhs_fmm2d, rhs, atol = 1e-4);
@time Krylov.gmres(lhs_direct, rhs, atol = 1e-4);

function profile_fmm2d(n)
    for i in 1:n
        Krylov.gmres(lhs_fmm2d, rhs, atol = 1e-4)
    end
end

function profile_direct(n)
    for i in 1:n
        Krylov.gmres(lhs_direct, rhs, atol = 1e-4)
    end
end

ProfileSVG.@profview profile_fmm2d(5)
ProfileSVG.@profview profile_direct(5)
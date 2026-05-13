# Quick debug: solve on the unrefined cube with a small p, plot σ at Nyström
# nodes on the +z face, compare against the dense ground-truth operator.

using LinearAlgebra
using FastGaussQuadrature
using BoundaryIntegral
const BI = BoundaryIntegral

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0

p = 8
l_ec = 0.25  # one edge-refinement level
ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)

iface = single_dielectric_box3d(Lx, Ly, Lz, p, l_ec, eps_d, eps_0)
npan = length(iface.panels)
N = npan * p^2
println("panels = $npan,  N = $N")

# Plain FMM and DENSE for comparison
Lhs_fmm = lhs_dielectric_box3d_fmm3d(iface, 1e-12)
Lhs_dense = lhs_dielectric_box3d(iface)
rhs = rhs_dielectric_box3d(iface, ps, eps_src)

sigma_fmm = BI.solve_gmres(Lhs_fmm, rhs, 1e-10, 1e-10)
sigma_lu  = BI.solve_lu(Lhs_dense, rhs)

@info "Residuals" norm(Lhs_fmm*sigma_fmm - rhs) / norm(rhs)  norm(Lhs_dense*sigma_lu - rhs) / norm(rhs)
@info "FMM vs LU sigma diff" norm(sigma_fmm - sigma_lu) / norm(sigma_lu)

# Stats on σ on the +z face
top_idx = Int[]
top_pts = NTuple{3, Float64}[]
top_sigma = Float64[]
offset = 0
for (i, panel) in enumerate(iface.panels)
    is_top = all(abs(c[3] - Lz/2) < 1e-12 for c in panel.corners)
    if is_top
        for k in 1:length(panel.points)
            push!(top_idx, offset + k)
            push!(top_pts, panel.points[k])
            push!(top_sigma, sigma_fmm[offset + k])
        end
    end
    offset += length(panel.points)
end

println("\nTop-face σ stats:")
println("  count = ", length(top_sigma))
println("  min   = ", minimum(top_sigma))
println("  max   = ", maximum(top_sigma))
println("  |min| = ", minimum(abs.(top_sigma)))
println("  |max| = ", maximum(abs.(top_sigma)))
println("  mean  = ", sum(top_sigma)/length(top_sigma))

# Look at σ values at points near the edge x=1 (top face)
println("\nσ near +x edge of +z face (y close to 0):")
for (k, p) in enumerate(top_pts)
    if abs(p[2]) < 0.1 && p[1] > 0
        println("  pt = ($(round(p[1]; digits=4)), $(round(p[2]; digits=4))), σ = $(round(top_sigma[k]; digits=6))")
    end
end

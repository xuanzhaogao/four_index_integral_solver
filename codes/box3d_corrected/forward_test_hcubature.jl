include("utils.jl")

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using HCubature

# Forward test: direct HCubature DT*sigma on full surfaces vs DT_corrected*sigma
Lx = 20.0
Ly = 20.0
Lz = 0.5

p = 4
r = 2

eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^r * 1.01

# Build interface only; no solve needed for prescribed sigma
# (use non-adaptive mesh to keep mapping simple)
tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, p, l_ec, eps_in, eps_out)

# Surface density on all faces
# sigma_func = (x, y, z) -> sin(x) * exp(-y^2) * cos(z)
sigma_func = (x, y, z) -> 1.0

panel_sizes = [length(p.points) for p in tbox.panels]
offsets = cumsum([0; panel_sizes[1:end-1]])

panel_indices(i) = (offsets[i] + 1):(offsets[i] + panel_sizes[i])

sigma = zeros(Float64, BI.num_points(tbox))
for (i, pt) in enumerate(BI.eachpoint(tbox))
    x, y, z = pt.panel_point.point
    sigma[i] = sigma_func(x, y, z)
end

fmm_tol = 1e-10
up_tol = 1e-8
range_factor = 6.0
hcubature_atol = 1e-6

DT_direct = BI.laplace3d_DT_fmm3d(tbox, fmm_tol)

DT_corrected = BI.laplace3d_DT_fmm3d_corrected(
    tbox,
    fmm_tol,
    up_tol,
    max_order,
    include_edges_src = false,
    include_edges_trg = false,
    range_factor = range_factor,
)

rd = DT_direct * sigma
rc = DT_corrected * sigma

# Direct HCubature evaluation on the full surfaces (face-wise, not panel-wise)
faces = [
    (( Lx / 2,  Ly / 2,  Lz / 2), (-Lx / 2,  Ly / 2,  Lz / 2), (-Lx / 2, -Ly / 2,  Lz / 2), ( Lx / 2, -Ly / 2,  Lz / 2)), # z=+Lz/2
    # (( Lx / 2,  Ly / 2, -Lz / 2), ( Lx / 2, -Ly / 2, -Lz / 2), (-Lx / 2, -Ly / 2, -Lz / 2), (-Lx / 2,  Ly / 2, -Lz / 2)), # z=-Lz/2
    (( Lx / 2, -Ly / 2, -Lz / 2), ( Lx / 2,  Ly / 2, -Lz / 2), ( Lx / 2,  Ly / 2,  Lz / 2), ( Lx / 2, -Ly / 2,  Lz / 2)), # x=+Lx/2
    ((-Lx / 2, -Ly / 2, -Lz / 2), (-Lx / 2, -Ly / 2,  Lz / 2), (-Lx / 2,  Ly / 2,  Lz / 2), (-Lx / 2,  Ly / 2, -Lz / 2)), # x=-Lx/2
    ((-Lx / 2,  Ly / 2, -Lz / 2), (-Lx / 2,  Ly / 2,  Lz / 2), ( Lx / 2,  Ly / 2,  Lz / 2), ( Lx / 2,  Ly / 2, -Lz / 2)), # y=+Ly/2
    ((-Lx / 2, -Ly / 2, -Lz / 2), ( Lx / 2, -Ly / 2, -Lz / 2), ( Lx / 2, -Ly / 2,  Lz / 2), (-Lx / 2, -Ly / 2,  Lz / 2)), # y=-Ly/2
]

function face_integral(trg_point::NTuple{3, Float64}, trg_normal::NTuple{3, Float64}, face)
    a, b, c, d = face
    cc = (a .+ b .+ c .+ d) ./ 4
    bma = b .- a
    dma = d .- a
    Lx_face = norm(b .- a)
    Ly_face = norm(d .- a)
    scale = Lx_face * Ly_face / 4
    function integrand(x)
        u, v = x
        p = cc .+ bma .* (u / 2) .+ dma .* (v / 2)
        return BI.laplace3d_grad(p, trg_point, trg_normal) * sigma_func(p[1], p[2], p[3]) * scale
    end
    res, _ = hcubature(integrand, [-1.0, -1.0], [1.0, 1.0]; atol = hcubature_atol)
    return res
end


n_trg = 40
llx = Lx / 2 - 0.1
lly = Ly / 2 - 0.1
xs = range(-llx, llx, length = n_trg)
ys = range(-lly, lly, length = n_trg)
z_trg = -Lz / 2
trg_normal = (0.0, 0.0, -1.0)

res_hcub = zeros(Float64, n_trg, n_trg)
res_direct = zeros(Float64, n_trg, n_trg)
res_corrected = zeros(Float64, n_trg, n_trg)

res_direct_approx = BI.interface_approx(tbox, rd)
res_corrected_approx = BI.interface_approx(tbox, rc)

for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
    @show i, j
    trg_point = (x, y, z_trg)
    acc = 0.0
    for face in faces
        acc += face_integral(trg_point, trg_normal, face)
    end
    res_hcub[i, j] = acc
    res_direct[i, j] = res_direct_approx(trg_point)
    res_corrected[i, j] = res_corrected_approx(trg_point)
end

maximum(abs.(res_hcub .- res_corrected))
maximum(abs.(res_hcub .- res_direct))


begin
    fig = Figure(size = (500, 500), fontsize = 20)
    ax = Axis(fig[1, 1], aspect = DataAspect())
    hm1 = heatmap!(ax, xs, ys, log10.(abs.(res_corrected .- res_hcub)); colormap = :viridis)
    Colorbar(fig[1, 2], hm1)

    save("figs/forward_test_hcubature.png", fig)

    fig
end
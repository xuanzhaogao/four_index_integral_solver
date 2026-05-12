# Diagnostic: for a single corner panel pair, compare
#   A_raw[k,j]      = ω_j K(s_j, t_k)               (GGQ-raw)
#   A_corr[k,j]     = HCubature on Müntz Lagrange ψ_j
#   row sum agreement on a known density σ in the m-mode Müntz span

include(joinpath(@__DIR__, "..", "src", "gram.jl"))
include(joinpath(@__DIR__, "ggq_corner_box2d.jl"))

using HCubature, Printf

function diag()
    α1 = -0.329682
    γshift = [0.0, 0.6594, 2.0, 2.6594, 4.0, 4.6594, 6.0, 6.6594]   # m=8
    m  = length(γshift)

    Cbig, _, _ = generalized_power_basis(BigFloat(α1), BigFloat.(γshift))
    Cf  = Float64.(Cbig)

    # GGQ rule for these gammas
    γs_full = [0.6703178814, 1.3296821186, 2.6703178814, 3.3296821186,
               4.6703178814, 5.3296821186, 6.6703178814, 7.3296821186,
               8.6703178814, 9.3296821186, 10.6703178814, 11.3296821186,
               12.6703178814, 13.3296821186, 14.6703178814, 15.3296821186]
    x01, w01, _, _ = build_ggq(γs_full, m; nsteps=40)
    @printf("Built GGQ m=%d, nodes=%s\n", m, round.(x01, digits=4))

    # Set up perpendicular panels at vertex (-0.5, 0.5).
    L_corner = 0.05
    vertex = (-0.5, 0.5)
    # Source: top-edge horizontal panel; goes vertex -> (vertex.x + L, vertex.y)
    src_dir   = (1.0, 0.0)
    src_normal = (0.0, 1.0)
    src_corners = [vertex, (vertex[1]+L_corner, vertex[2])]
    s_nodes  = L_corner .* x01
    src_pts  = [(vertex[1] + s, vertex[2]) for s in s_nodes]

    Qmat = Float64.(eval_Q_matrix(BigFloat.(x01), BigFloat.(γshift), Cbig))
    Minv = inv(Matrix(Qmat'))

    # FlatPanel for src
    src_panel = BI.FlatPanel(src_normal, src_corners, true,
                             m, collect(2 .* x01 .- 1), copy(w01),
                             src_pts,
                             [L_corner * w01[j] / x01[j]^α1 for j in 1:m])

    src_info = CornerPanelInfo(1, vertex, :a, L_corner,
                               src_corners[1], src_corners[2], src_dir,
                               s_nodes, copy(x01), copy(w01),
                               γshift, Cf, Minv)

    # Target: a node on the left-edge perpendicular panel (vertex going down)
    tgt_dir = (0.0, -1.0)
    target_x = x01   # use the same node positions on perpendicular partner
    targets = [(vertex[1], vertex[2] + L_corner * x * tgt_dir[2]) for x in target_x]

    # Test density: σ(s) = s^α × Q_k(s/L) for k = 1, ..., m  (the m basis modes)
    # σ values at GGQ nodes: σ_j = s_j^α × Q_k(x_j) = s_j^α × Qmat[k, j]
    println("\nMode-by-mode test: integral I_k = ∫_0^L K(s, t) σ(s) ds for σ = s^α Q_k(s/L)")
    println("Using HCubature on direct integrand (high accuracy reference).")
    println("Then compare to (a) GGQ-raw: Σ ω_j K(s_j, t) σ_j, (b) corrected matrix sum")

    # pick one target (k=1, the closest to vertex)
    target = targets[2]
    @printf("\nTarget: %s (distance from vertex = %.4f)\n", target,
            norm(target .- vertex))

    # Reference integral I_k via HCubature on s ∈ [0, L]
    function ref_integral_mode(k)
        integrand = function (yv)
            y = yv[1]
            y == 0.0 && return 0.0
            s = L_corner * y
            src_pt = (vertex[1] + s, vertex[2])
            K = dt_kernel(src_pt, target, src_normal)
            # Q_k(y) = Σ_l C[k, l] y^{γshift[l]}
            Qk = 0.0
            for l in 1:m
                γ = γshift[l]
                Qk += Cf[k, l] * (γ == 0.0 ? 1.0 : y^γ)
            end
            return K * y^α1 * Qk * L_corner
        end
        val, _ = hquadrature(integrand, 0.0, 1.0; rtol=1e-12, atol=1e-16, maxevals=100_000)
        return val
    end

    # GGQ-raw: σ_j = s_j^α Q_k(x_j); sum Σ ω_j K(s_j, t) σ_j
    function ggq_raw_mode(k)
        s = 0.0
        for j in 1:m
            σ_j = (s_nodes[j])^α1 * Qmat[k, j]
            ω_j = L_corner * w01[j] / x01[j]^α1
            K_j = dt_kernel(src_pts[j], target, src_normal)
            s += ω_j * K_j * σ_j
        end
        return s
    end

    # Corrected: same σ_j vector × corrected entries
    function corrected_mode(k)
        s = 0.0
        for j in 1:m
            σ_j = (s_nodes[j])^α1 * Qmat[k, j]
            A_kj = corrected_entry(target, src_panel, src_info, j, α1; rtol=1e-11)
            s += A_kj * σ_j
        end
        return s
    end

    @printf("%-3s  %-18s  %-18s  %-18s  %-12s %-12s\n",
            "k", "Reference", "GGQ-raw", "Corrected",
            "raw - ref", "corr - ref")
    for k in 1:m
        Iref  = ref_integral_mode(k)
        Iraw  = ggq_raw_mode(k)
        Icorr = corrected_mode(k)
        @printf("%-3d  %+.10e  %+.10e  %+.10e  %.2e   %.2e\n",
                k, Iref, Iraw, Icorr, abs(Iraw - Iref), abs(Icorr - Iref))
    end
end

diag()

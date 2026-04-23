"""
Benchmark: GL vs GL+GJ convergence on 1×1 dielectric box.

Geometry: 1×1 box, ε_in=1, ε_out=200, point source at (0.001, 0.001).
Accuracy: Gauss's law  |dot(ws, x) - (1/ε_out - 1/ε_in)|.

Two methods:
  GL     : Gauss-Legendre on all panels (use_singular=false)
  GL+GJ  : Gauss-Jacobi on corner panels, GL on the rest (use_singular=true)

Varied parameters:
  p     = n_quad ∈ {4, 6, 8, 10, 12, 16, 20}
  level = h-refinement level; l_panel and l_corner scale proportionally:
    level 1 (coarse)   : l_panel=0.20, l_corner=0.050
    level 2 (medium)   : l_panel=0.10, l_corner=0.025
    level 3 (fine)     : l_panel=0.05, l_corner=0.0125
    level 4 (very fine): l_panel=0.025, l_corner=0.00625
"""

import BoundaryIntegral as BI
using LinearAlgebra, Printf

function main()
    # ── problem setup ────────────────────────────────────────────────────────
    eps_in  = 1.0
    eps_out = 200.0
    ps      = BI.PointSource((0.001, 0.001), 1.0)
    exact_flux = 1.0/eps_out - 1.0/eps_in
    gamma  = BI.corner_singularity_power(pi/2, eps_in, eps_out)

    println("ε_in=$eps_in, ε_out=$eps_out, γ=$(round(gamma, sigdigits=4))")
    println("exact_flux = $exact_flux")
    println()

    # ── sweep parameters ─────────────────────────────────────────────────────
    p_vals = [4, 6, 8, 10, 12, 16, 20]

    levels = [
        (1, 0.20, 0.050),
        (2, 0.10, 0.025),
        (3, 0.05, 0.0125),
        (4, 0.025, 0.00625),
    ]

    # ── helper ───────────────────────────────────────────────────────────────
    function run_case(nq, l_panel, l_corner, use_singular)
        box = BI.single_dielectric_box2d(1.0, 1.0, nq, l_panel, l_corner,
                                         eps_in, eps_out, Float64;
                                         use_singular=use_singular)
        lhs = BI.lhs_dielectric_box2d(box)
        rhs = BI.rhs_dielectric_box2d(box, ps, eps_in)
        x   = BI.solve_lu(lhs, rhs)
        err  = abs(dot(BI.all_weights(box), x) - exact_flux)
        npts = BI.num_points(box)
        return err, npts
    end

    # ── per-method error tables ───────────────────────────────────────────────
    for method in [:GL, :GLGJ]
        use_singular = (method == :GLGJ)
        name = (method == :GL) ? "Pure GL" : "GL + GJ (corner panels)"
        println("="^70)
        println("Method: $name")
        println("="^70)

        hdr = @sprintf("%-6s", "p")
        for (lv, lp, lc) in levels
            hdr *= @sprintf("  %-20s", " lv$lv(lp=$(lp))")
        end
        println(hdr)
        println("-"^(6 + length(levels)*22))

        for nq in p_vals
            row = @sprintf("%-6d", nq)
            for (lv, lp, lc) in levels
                err, npts = run_case(nq, lp, lc, use_singular)
                row *= @sprintf("  %-8.2e(n=%5d)  ", err, npts)
            end
            println(row)
            flush(stdout)
        end
        println()
    end

    # ── ratio table: GL_err / GJ_err ─────────────────────────────────────────
    println("="^70)
    println("Ratio table: GL_err / GJ_err  (>1 means GJ is better)")
    println("="^70)
    ratio_hdr = @sprintf("%-6s", "p")
    for (lv, lp, lc) in levels
        ratio_hdr *= @sprintf("  %-14s", "lv$lv(lp=$(lp))")
    end
    println(ratio_hdr)
    println("-"^(6 + length(levels)*16))

    for nq in p_vals
        ratio_row = @sprintf("%-6d", nq)
        for (lv, lp, lc) in levels
            err_gl, _ = run_case(nq, lp, lc, false)
            err_gj, _ = run_case(nq, lp, lc, true)
            ratio_row *= @sprintf("  %-14.2f", err_gl / max(err_gj, 1e-20))
        end
        println(ratio_row)
        flush(stdout)
    end
    println()
end

main()

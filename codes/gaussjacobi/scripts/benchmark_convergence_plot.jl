"""
Run GL vs GL+GJ benchmark, save results to CSV, and generate convergence plots.
"""

import BoundaryIntegral as BI
using LinearAlgebra, Printf, CSV, DataFrames, Plots

function main()
    # ── problem setup ────────────────────────────────────────────────────────
    eps_in  = 1.0
    eps_out = 200.0
    ps      = BI.PointSource((0.001, 0.001), 1.0)
    exact_flux = 1.0/eps_out - 1.0/eps_in

    p_vals = [4, 6, 8, 10, 12, 16, 20]

    levels = [
        (1, 0.20, 0.050),
        (2, 0.10, 0.025),
        (3, 0.05, 0.0125),
        (4, 0.025, 0.00625),
    ]

    function run_case(nq, l_panel, l_corner, use_singular)
        box = BI.single_dielectric_box2d(1.0, 1.0, nq, l_panel, l_corner,
                                         eps_in, eps_out, Float64;
                                         use_singular=use_singular)
        lhs = BI.lhs_dielectric_box2d(box)
        rhs = BI.rhs_dielectric_box2d(box, ps, eps_in)
        x   = BI.solve_lu(lhs, rhs)
        abs(dot(BI.all_weights(box), x) - exact_flux), BI.num_points(box)
    end

    # ── collect results ───────────────────────────────────────────────────────
    rows = []
    for (lv, lp, lc) in levels
        for nq in p_vals
            for (method, use_singular) in [("GL", false), ("GLGJ", true)]
                err, npts = run_case(nq, lp, lc, use_singular)
                @printf("method=%-5s  lv=%d  p=%-2d  n=%5d  err=%.3e\n",
                        method, lv, nq, npts, err)
                flush(stdout)
                push!(rows, (method=method, level=lv, l_panel=lp, l_corner=lc,
                             p=nq, n_pts=npts, error=err))
            end
        end
    end

    df = DataFrame(rows)
    csv_path = joinpath(@__DIR__, "benchmark_convergence.csv")
    CSV.write(csv_path, df)
    println("\nSaved: $csv_path")

    # ── plots ─────────────────────────────────────────────────────────────────
    markers = [:circle, :square, :diamond, :utriangle]
    colors_gl   = [:steelblue, :dodgerblue, :royalblue, :navy]
    colors_glgj = [:tomato, :orangered, :firebrick, :darkred]

    # --- Figure 1: error vs p, one line per refinement level, both methods ---
    fig1 = plot(;
        xlabel = "p  (quadrature order)",
        ylabel = "Gauss's law error",
        title  = "Convergence vs p  (ε_in=1, ε_out=200)",
        yscale = :log10,
        legend = :outertopright,
        size   = (800, 500),
        framestyle = :box,
    )
    for (i, (lv, lp, lc)) in enumerate(levels)
        sub_gl   = filter(r -> r.method == "GL"   && r.level == lv, df)
        sub_glgj = filter(r -> r.method == "GLGJ" && r.level == lv, df)
        plot!(fig1, sub_gl.p,   sub_gl.error;
              label  = "GL lv$lv (lp=$lp)",
              color  = colors_gl[i], marker = markers[i], linestyle = :solid)
        plot!(fig1, sub_glgj.p, sub_glgj.error;
              label  = "GL+GJ lv$lv (lp=$lp)",
              color  = colors_glgj[i], marker = markers[i], linestyle = :dash)
    end
    savefig(fig1, joinpath(@__DIR__, "benchmark_error_vs_p.pdf"))
    savefig(fig1, joinpath(@__DIR__, "benchmark_error_vs_p.png"))
    println("Saved: benchmark_error_vs_p.pdf / .png")

    # --- Figure 2: error vs n_pts (work-accuracy), both methods ---
    fig2 = plot(;
        xlabel = "n_pts  (total DOF)",
        ylabel = "Gauss's law error",
        title  = "Work-accuracy  (ε_in=1, ε_out=200)",
        xscale = :log10,
        yscale = :log10,
        legend = :outertopright,
        size   = (800, 500),
        framestyle = :box,
    )
    for (i, (lv, lp, lc)) in enumerate(levels)
        sub_gl   = filter(r -> r.method == "GL"   && r.level == lv, df)
        sub_glgj = filter(r -> r.method == "GLGJ" && r.level == lv, df)
        plot!(fig2, sub_gl.n_pts,   sub_gl.error;
              label  = "GL lv$lv",
              color  = colors_gl[i], marker = markers[i], linestyle = :solid)
        plot!(fig2, sub_glgj.n_pts, sub_glgj.error;
              label  = "GL+GJ lv$lv",
              color  = colors_glgj[i], marker = markers[i], linestyle = :dash)
    end
    savefig(fig2, joinpath(@__DIR__, "benchmark_work_accuracy.pdf"))
    savefig(fig2, joinpath(@__DIR__, "benchmark_work_accuracy.png"))
    println("Saved: benchmark_work_accuracy.pdf / .png")

    # --- Figure 3: for fixed p, error vs refinement level ---
    p_select = [4, 8, 12, 20]
    colors_p = [:green, :darkorange, :purple, :black]
    fig3 = plot(;
        xlabel = "refinement level",
        ylabel = "Gauss's law error",
        title  = "h-convergence at fixed p  (ε_in=1, ε_out=200)",
        yscale = :log10,
        xticks = (1:4, ["lv1\n(lp=0.20)", "lv2\n(lp=0.10)", "lv3\n(lp=0.05)", "lv4\n(lp=0.025)"]),
        legend = :outertopright,
        size   = (800, 500),
        framestyle = :box,
    )
    for (i, p) in enumerate(p_select)
        sub_gl   = filter(r -> r.method == "GL"   && r.p == p, df)
        sub_glgj = filter(r -> r.method == "GLGJ" && r.p == p, df)
        sort!(sub_gl,   :level)
        sort!(sub_glgj, :level)
        plot!(fig3, sub_gl.level,   sub_gl.error;
              label  = "GL p=$p",
              color  = colors_p[i], marker = :circle, linestyle = :solid)
        plot!(fig3, sub_glgj.level, sub_glgj.error;
              label  = "GL+GJ p=$p",
              color  = colors_p[i], marker = :square, linestyle = :dash)
    end
    savefig(fig3, joinpath(@__DIR__, "benchmark_h_convergence.pdf"))
    savefig(fig3, joinpath(@__DIR__, "benchmark_h_convergence.png"))
    println("Saved: benchmark_h_convergence.pdf / .png")

    return df
end

main()

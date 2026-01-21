using Random
using LinearAlgebra
using BoundaryIntegral
using FastGaussQuadrature
using LegendrePolynomials
using CSV
using DataFrames

function make_panel(n_quad, a, b, c, d)
    ns, ws = gausslegendre(n_quad)
    normal = (0.0, 0.0, 1.0)
    return BoundaryIntegral.rect_panel3d_discretize(a, b, c, d, ns, ws, normal; is_edge=false)
end

function combo_error(panel, trg; atol=1e-8, max_order=64, n_trials=5, rng=MersenneTwister(0))
    n_quad_up = BoundaryIntegral.check_quad_order3d(panel, trg, atol, max_order)
    n_ref = min(max_order * 2, max(n_quad_up * 4, n_quad_up))

    M_up = BoundaryIntegral.int_laplace3d_grad(panel.n_quad, n_quad_up, panel, trg)

    errs = Float64[]
    for _ in 1:n_trials
        coeffs = randn(rng, panel.n_quad + 1, panel.n_quad + 1)
        val_up = sum(M_up .* coeffs)
        val_ref = combo_integral(panel, trg, n_ref, coeffs)
        push!(errs, abs(val_up - val_ref))
    end

    return n_quad_up, n_ref, errs
end

function combo_integral(panel, trg, n_quad_up, coeffs)
    ns, ws = gausslegendre(n_quad_up)
    a, b, c, d = panel.corners
    cc = (a .+ b .+ c .+ d) ./ 4
    Lx = norm(b .- a)
    Ly = norm(c .- a)

    acc = 0.0
    for k in 1:n_quad_up
        x = ns[k]
        for l in 1:n_quad_up
            y = ns[l]
            fxy = 0.0
            for i in 0:panel.n_quad
                pix = Pl(x, i)
                for j in 0:panel.n_quad
                    fxy += coeffs[i + 1, j + 1] * pix * Pl(y, j)
                end
            end
            p = cc .+ (b .- a) .* (x / 2) .+ (d .- a) .* (y / 2)
            acc += ws[k] * ws[l] * fxy * BoundaryIntegral.laplace3d_grad(p, trg, panel.normal) * Lx * Ly / 4
        end
    end
    return acc
end

function run_suite()
    n_quads = [2, 4, 6]
    panels = [
        ((-0.5, -0.5, 0.0), (0.5, -0.5, 0.0), (0.5, 0.5, 0.0), (-0.5, 0.5, 0.0)), # square
        ((-0.6, -0.5, 0.0), (0.6, -0.5, 0.0), (0.6, 0.5, 0.0), (-0.6, 0.5, 0.0)), # 1.2:1
        ((-0.5, -0.6, 0.0), (0.5, -0.6, 0.0), (0.5, 0.6, 0.0), (-0.5, 0.6, 0.0)), # 1:1.2
        ((-0.7071, -0.5, 0.0), (0.7071, -0.5, 0.0), (0.7071, 0.5, 0.0), (-0.7071, 0.5, 0.0)), # sqrt(2):1
    ]
    trgs = [
        (0.1, 0.1, 0.01),  # near
        (0.2, -0.15, 0.1), # mid
        (0.3, 0.25, 0.5),  # far
    ]
    rows = Vector{NamedTuple}()
    for atol in [1e-4, 1e-6, 1e-8]
        max_order = 1024
        println("atol=", atol)
        for (pidx, (a, b, c, d)) in enumerate(panels)
            println("panel=", pidx, " corners=", (a, b, c, d))
            for n_quad in n_quads
                panel = make_panel(n_quad, a, b, c, d)
                println("  n_quad=", n_quad)
                for trg in trgs
                    n_quad_up, n_ref, errs = combo_error(panel, trg; atol=atol, max_order=max_order)
                    println("    trg=", trg, " n_quad_up=", n_quad_up, " n_ref=", n_ref, " errs=", errs)
                    @assert length(errs) == 5
                    push!(
                        rows,
                        (
                            atol = atol,
                            panel_idx = pidx,
                            a_x = a[1], a_y = a[2], a_z = a[3],
                            b_x = b[1], b_y = b[2], b_z = b[3],
                            c_x = c[1], c_y = c[2], c_z = c[3],
                            d_x = d[1], d_y = d[2], d_z = d[3],
                            n_quad = n_quad,
                            trg_x = trg[1], trg_y = trg[2], trg_z = trg[3],
                            n_quad_up = n_quad_up,
                            n_ref = n_ref,
                            err1 = errs[1], err2 = errs[2], err3 = errs[3], err4 = errs[4], err5 = errs[5],
                        ),
                    )
                end
            end
        end
    end
    df = DataFrame(rows)
    out_path = joinpath(@__DIR__, "data", "quad_order_check.csv")
    CSV.write(out_path, df)

    df.max_err = map(i -> maximum((
        df.err1[i],
        df.err2[i],
        df.err3[i],
        df.err4[i],
        df.err5[i],
    )), 1:nrow(df))

    summary = combine(
        groupby(df, [:atol, :panel_idx, :n_quad, :trg_x, :trg_y, :trg_z]),
        :max_err => maximum => :max_err,
        :n_quad_up => maximum => :n_quad_up,
        :n_ref => maximum => :n_ref,
    )
    println("\nMax error summary (per atol/panel/n_quad/target):")
    show(summary, allrows = true, allcols = true)
    println()
end

run_suite()

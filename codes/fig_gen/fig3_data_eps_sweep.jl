#=
Epsilon sweep for Figure 3 panel (b).

For contrasts alpha = (eps_d - 1) / (eps_d + 1) in {10/33, 20/33, 30/33},
solve the corrected dielectric-cube BIE on the same edge-local refinement
ladder. The accuracy proxy is the absolute weighted surface integral of the
surface density minus the expected interior-source value, 1 - 1 / eps_d.

Both N and l_min are serialized so the plot can use either quantity as its
x-axis without rerunning the sweep.
=#

using LinearAlgebra
using Serialization
using Printf
using BoundaryIntegral
const BI = BoundaryIntegral

include("fig3_eps_utils.jl")

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_0 = 1.0
const ps = BI.PointSource((0.1, 0.2, 0.3), 1.0)
const p_quad = 4
const k_list = collect(1:8)
const contrast_specs = [(num = 10, den = 33), (num = 20, den = 33), (num = 30, den = 33)]

const fmm_tol = 1e-6
const gmres_tol = 1e-6
const up_tol = 1e-6
const max_order = 64

const jls_path = joinpath(@__DIR__, "fig3_data_eps_sweep.jls")
const prog_path = joinpath(@__DIR__, "fig3_data_eps_sweep_progress.txt")

function solve_and_charge(eps_d::Float64, l_ec::Float64)
    iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
    Lhs = lhs_dielectric_box3d_fmm3d_corrected(iface, fmm_tol, up_tol, max_order;
                                               correct_edges = true)
    rhs = rhs_dielectric_box3d(iface, ps, eps_d)
    sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
    gres = norm(Lhs * sigma - rhs) / norm(rhs)
    q_total = charge_integral(BI.all_weights(iface), sigma)
    q_theory = charge_theory_interior_source(eps_d; eps_0 = eps_0, source_charge = ps.charge)
    return (
        npanels = length(iface.panels),
        N = length(iface.panels) * p_quad^2,
        l_min = panel_min_edge_length(iface),
        charge = q_total,
        charge_theory = q_theory,
        charge_error = q_total - q_theory,
        charge_error_abs = abs(q_total - q_theory),
        gmres_res = gres,
    )
end

function checkpoint(epsilon_results)
    out = (
        geometry = (Lx = Lx, Ly = Ly, Lz = Lz),
        eps_0 = eps_0,
        eps_src = "eps_d",
        source = (point = ps.point, charge = ps.charge),
        p_quad = p_quad,
        k_list = k_list,
        contrast_specs = contrast_specs,
        fmm_tol = fmm_tol,
        gmres_tol = gmres_tol,
        up_tol = up_tol,
        max_order = max_order,
        correct_edges = true,
        x_axes = (:N, :l_min),
        metric = "abs(sum(weights .* sigma) - (1 - 1 / eps_d))",
        epsilon_results = epsilon_results,
    )
    open(io -> serialize(io, out), jls_path, "w")
end

open(prog_path, "w") do io
    println(io, "# fig3 panel(b) epsilon sweep: metric=abs(sum(weights .* sigma) - (1 - 1 / eps_d))")
    println(io, "# source=$(ps.point) charge=$(ps.charge) eps_src=eps_d p_quad=$p_quad fmm_tol=$fmm_tol gmres_tol=$gmres_tol up_tol=$up_tol max_order=$max_order")
    println(io, "# alpha        eps_d        k   l_ec        l_min       npanels   N          charge            charge_theory     abs_error         gmres_res     dt[s]")
    flush(io)
end

epsilon_results = NamedTuple[]
for spec in contrast_specs
    alpha = spec.num / spec.den
    eps_d = Float64(epsilon_from_contrast(alpha; eps_0 = eps_0))
    @info "===== epsilon sweep =====" alpha eps_d
    levels = NamedTuple[]
    for k in k_list
        l_ec = 1.01 / 2.0^k
        @info "  level" alpha eps_d k l_ec
        t0 = time()
        res = solve_and_charge(eps_d, l_ec)
        dt = time() - t0
        level = (
            k = k,
            l_ec = l_ec,
            l_min = res.l_min,
            npanels = res.npanels,
            N = res.N,
            charge = res.charge,
            charge_theory = res.charge_theory,
            charge_error = res.charge_error,
            charge_error_abs = res.charge_error_abs,
            gmres_res = res.gmres_res,
            dt = dt,
        )
        push!(levels, level)
        @info "  done" alpha eps_d k res.N res.l_min res.charge res.charge_theory res.charge_error_abs res.gmres_res dt

        checkpoint(vcat(epsilon_results, [(
            contrast_num = spec.num,
            contrast_den = spec.den,
            alpha = alpha,
            eps_d = eps_d,
            levels = levels,
        )]))
        open(prog_path, "a") do io
            @printf(io, "  %-10.6g   %-10.6g   %d   %-9.4g   %-9.4g   %-7d   %-9d   %-16.8e   %-16.8e   %-16.8e   %-11.3g   %.1f\n",
                    alpha, eps_d, k, l_ec, res.l_min, res.npanels, res.N,
                    res.charge, res.charge_theory, res.charge_error_abs, res.gmres_res, dt)
            flush(io)
        end
    end
    push!(epsilon_results, (
        contrast_num = spec.num,
        contrast_den = spec.den,
        alpha = alpha,
        eps_d = eps_d,
        levels = levels,
    ))
end

checkpoint(epsilon_results)
@info "Saved" jls_path bytes = stat(jls_path).size

println("\n=== Figure 3 epsilon sweep ===")
for r in epsilon_results
    println("\n  alpha = $(r.contrast_num)/$(r.contrast_den), eps_d = $(round(r.eps_d; sigdigits = 6))")
    println("    k    l_min       N          charge           theory           abs_error        gmres_res")
    for l in r.levels
        println("    $(l.k)    $(round(l.l_min; sigdigits = 4))   $(l.N)    $(round(l.charge; sigdigits = 4))    $(round(l.charge_theory; sigdigits = 4))    $(round(l.charge_error_abs; sigdigits = 4))    $(round(l.gmres_res; sigdigits = 3))")
    end
end

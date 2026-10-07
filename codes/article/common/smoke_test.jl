# Smoke test for Harness.jl: tiny end-to-end runs validating API assumptions.
# Run: julia --project=. common/smoke_test.jl

include(joinpath(@__DIR__, "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, Printf

const H = Harness

println("== 1. FMM3D tolerance floor probe (via BI corrected operator) ==")
let sys = system1()
    res0 = solve_system(sys; eps = 1e-3, p = 4, r = 2)
    x = randn(BI.num_points(res0.interface))
    Dref = BI.laplace3d_DT(res0.interface) * x   # dense direct reference
    for tol in (1e-9, 1e-12, 1e-13, 1e-14)
        out = try
            D = BI.laplace3d_DT_fmm3d(res0.interface, tol)
            err = norm(D * x .- Dref, Inf) / norm(Dref, Inf)
            @sprintf("ok  relerr=%.2e", err)
        catch e
            "FAIL ($(sprint(showerror, e)[1:min(end,80)]))"
        end
        @printf("  fmm tol=%.0e : %s\n", tol, out)
    end
end

println("== 2. System I tiny solve (eps=1e-3, p=4, r=2) ==")
sys1 = system1()
t = @elapsed res = solve_system(sys1; eps = 1e-3, p = 4, r = 2)
@printf("  N=%d niter=%d residual=%.2e near=%d adaptive=%d p_up_max=%d corr_mb=%.2f t=%.1fs\n",
        res.N, res.niter, res.residual, res.n_near_pairs, res.n_adaptive_pairs,
        res.p_up_max, res.corr_bytes / 1e6, t)
println("  stats fields: ", propertynames(res.stats))
println("  history length: ", length(res.gmres_history))

println("== 3. V + phi evaluation ==")
tdict = Dict{String, Float64}()
V = eval_V(res; t_out = tdict)
@printf("  V = %.8e   (t_inc=%.2fs t_sc=%.2fs)\n", V, tdict["eval_incident"], tdict["eval_scatter"])
zt = zone_targets(sys1)
phi_near = eval_phi(res, zt.near)
phi_supp = eval_phi(res, zt.supp)
phi_far = eval_phi(res, zt.far)
@printf("  |phi| near/supp/far = %.3e / %.3e / %.3e\n",
        norm(phi_near, Inf), norm(phi_supp, Inf), norm(phi_far, Inf))

println("== 4. Vacuum V sanity (analytic check) ==")
# For two unit Gaussians sigma s, dist a: V = erf(a/(2s)) / a with kernel 1/(4pi r)...
# code kernel: laplace3d_pot = 1/(4pi r). erf-formula: V = erf(a/(2s))/(4pi a)
using SpecialFunctions
Vvac = vacuum_V(sys1, 1e-6)
a = norm(collect(sys1.src_center) .- collect(sys1.tgt_center))
Vana = erf(a / (2 * sys1.src_sigma)) / (4pi * a)
@printf("  V_vac = %.10e  analytic = %.10e  rel = %.2e\n", Vvac, Vana, abs(Vvac - Vana) / abs(Vana))

println("== 5. System II tiny solve + shared-face census (6.1.7) ==")
sys2 = system2()
res2 = solve_system(sys2; eps = 1e-3, p = 4, r = 2)
@printf("  N=%d niter=%d residual=%.2e\n", res2.N, res2.niter, res2.residual)
# census: panels lying in the x=0 plane between the two substrate boxes
let n_shared = 0, gammas = Float64[]
    for (i, p) in enumerate(res2.interface.panels)
        c = (p.corners[1] .+ p.corners[2] .+ p.corners[3] .+ p.corners[4]) ./ 4
        if abs(c[1]) < 1e-12 && abs(p.normal[1]) > 0.999 && c[3] < 0.5
            n_shared += 1
            ei, eo = res2.interface.eps_in[i], res2.interface.eps_out[i]
            push!(gammas, (ei - eo) / (ei + eo))
        end
    end
    println("  shared-face panels: ", n_shared, "  gamma values: ", unique(round.(gammas; digits = 6)))
    println("  expected gamma12 = +/-", round((4.0 - 12.0) / (4.0 + 12.0); digits = 6))
end

println("== 6. Variant B (no near correction) + uniform/greedy refinement ==")
resB = solve_system(sys1; eps = 1e-3, p = 4, r = 2, near_correction = false)
@printf("  B: N=%d niter=%d\n", resB.N, resB.niter)
iso_uniform = uniform_refine(res.interface, 1)
@printf("  uniform_refine x1: %d -> %d points\n", res.N, BI.num_points(iso_uniform))
iso_greedy = refine_to_dof(res.interface, round(Int, 1.5 * res.N))
@printf("  refine_to_dof 1.5x: %d points (target %d)\n", BI.num_points(iso_greedy), round(Int, 1.5 * res.N))

println("== 7. h0-controlled scattered evaluation ==")
X = zt.supp[:, 1:20]
u_default = H.eval_scatter(res, X)
u_h0, st = eval_scatter_with_h0(res, X, Inf)
@printf("  no-refine vs default evaluator: maxdiff = %.2e (n_hcub=%d)\n",
        norm(u_default .- u_h0, Inf), st.n_hcub)
u_h02, st2 = eval_scatter_with_h0(res, X, 0.05)
@printf("  h0=0.05: n_refined=%d n_hcub=%d  maxdiff vs default = %.2e\n",
        st2.n_refined, st2.n_hcub, norm(u_default .- u_h02, Inf))

println("== 8. Bernstein rho_min ==")
slab = slab_system()
res_slab = solve_system(slab; eps = 1e-3, p = 4, r = 2)
rho = bernstein_rho_min(res_slab.interface, res_slab.nb)
@printf("  slab: N=%d near=%d rho_min=%.4f -> stagnation(p=8) ~ %.2e\n",
        res_slab.N, res_slab.n_near_pairs, rho, rho^(-16))

println("== 9. CSV row ==")
append_csv_row(joinpath(@__DIR__, "smoke_test_row.csv"), run_cols(res; V = V, note = "smoke"))
println(read(joinpath(@__DIR__, "smoke_test_row.csv"), String))
rm(joinpath(@__DIR__, "smoke_test_row.csv"))

println("SMOKE TEST DONE")

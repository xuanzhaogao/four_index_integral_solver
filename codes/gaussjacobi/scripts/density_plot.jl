"""
Plot density σ(s) near a corner by interpolating each solution
onto a common fine uniform mesh, then comparing.

For each panel, we use barycentric interpolation from the panel's
quadrature nodes to evaluate the density at arbitrary points.
For GJ panels, the unknown is σ_bare; we interpolate σ_bare and
then multiply by (1+t)^α to get the physical σ.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

import BoundaryIntegral as BI
using LinearAlgebra, Printf, CairoMakie

# ── problem setup ────────────────────────────────────────────────────────────
eps_in  = 1.0
eps_out = 200.0
ps      = BI.PointSource((0.33, 0.44), 1.0)
exact_flux = 1.0/eps_out - 1.0/eps_in
gamma = BI.corner_singularity_power(pi/2, eps_in, eps_out)
alpha = gamma - 1
@printf("γ = %.6f,  α = %.6f\n", gamma, alpha)

corner = (-0.5, 0.5)       # top-left corner
edge_dir = (1.0, 0.0)      # edge going right along y=0.5

# ── helpers ──────────────────────────────────────────────────────────────────
function solve_box(box)
    lhs = BI.lhs_dielectric_box2d(box)
    rhs = BI.rhs_dielectric_box2d(box, ps, eps_in)
    return BI.solve_lu(lhs, rhs)
end

"""
Find panels on a given edge from `corner_pt` in direction `edge_dir`.
"""
function find_edge_panels(box, x, corner_pt, edge_dir; tol=1e-8)
    result = []
    offset = 0
    for p in box.panels
        np = BI.num_points(p)
        x_panel = x[offset+1:offset+np]

        panel_vec = p.corners[2] .- p.corners[1]
        panel_len = norm(panel_vec)
        panel_unit = panel_vec ./ panel_len
        along_edge = abs(abs(dot(panel_unit, edge_dir)) - 1.0) < 0.1

        if along_edge
            on_edge_a = abs(p.corners[1][2] - corner_pt[2]) < tol
            on_edge_b = abs(p.corners[2][2] - corner_pt[2]) < tol
            if on_edge_a && on_edge_b
                s_a = dot(p.corners[1] .- corner_pt, edge_dir)
                s_b = dot(p.corners[2] .- corner_pt, edge_dir)
                push!(result, (panel=p, x_panel=copy(x_panel),
                               s_start=min(s_a, s_b), s_end=max(s_a, s_b)))
            end
        end
        offset += np
    end
    sort!(result, by=r -> r.s_start)
    return result
end

"""
Interpolate density at arc-length s_eval.
For GJ panels: interpolate σ_bare, then multiply by (1+t)^α.
Only evaluates within the convex hull of each panel's nodes.
"""
function interpolate_density(edge_panels, s_eval, corner_pt, edge_dir)
    sigma_eval = fill(NaN, length(s_eval))

    for ep in edge_panels
        p = ep.panel
        np = BI.num_points(p)
        s0, s1 = ep.s_start, ep.s_end

        # determine orientation: does t=-1 map to s_start or s_end?
        s_c1 = dot(p.corners[1] .- corner_pt, edge_dir)
        s_c2 = dot(p.corners[2] .- corner_pt, edge_dir)
        flipped = s_c1 > s_c2

        # node positions in s-space (for clamping eval range)
        s_nodes = Float64[]
        for j in 1:np
            push!(s_nodes, dot(p.points[j] .- corner_pt, edge_dir))
        end
        s_node_min = minimum(s_nodes)
        s_node_max = maximum(s_nodes)

        for (k, s) in enumerate(s_eval)
            # only evaluate within the node range (avoid extrapolation)
            (s_node_min - 1e-14 <= s <= s_node_max + 1e-14) || continue

            # map s to t ∈ [-1, 1]
            t = 2 * (s - s0) / (s1 - s0) - 1
            if flipped
                t = -t
            end

            r = zeros(np)
            BI.barycentric_row!(r, p.gl_xs, p.bary_weights, t)
            val = dot(r, ep.x_panel)

            if p.is_singular
                sigma_eval[k] = val * (1 + t)^p.singular_exponent
            else
                sigma_eval[k] = val
            end
        end
    end
    return sigma_eval
end

"""
Quad node positions in s-space.
"""
function quad_node_s(edge_panels, corner_pt, edge_dir)
    s_all = Float64[]
    for ep in edge_panels
        for pt in ep.panel.points
            push!(s_all, dot(pt .- corner_pt, edge_dir))
        end
    end
    return sort(s_all)
end

# ── build discretizations ────────────────────────────────────────────────────
@printf("Building reference GL (n=32, lp=0.01, lc=0.0002)...\n")
box_ref = BI.single_dielectric_box2d(1.0, 1.0, 32, 0.01, 0.0002, eps_in, eps_out; use_singular=false)
x_ref = solve_box(box_ref)
@printf("  N=%d, flux_err=%.3e\n", BI.num_points(box_ref),
        abs(dot(BI.all_weights(box_ref), x_ref) - exact_flux))

n_quad = 16; l_corner = 0.05
@printf("Building GL (n=%d, lc=%.3f)...\n", n_quad, l_corner)
box_gl = BI.single_dielectric_box2d(1.0, 1.0, n_quad, 0.2, l_corner, eps_in, eps_out; use_singular=false)
x_gl = solve_box(box_gl)
@printf("  N=%d, flux_err=%.3e\n", BI.num_points(box_gl),
        abs(dot(BI.all_weights(box_gl), x_gl) - exact_flux))

@printf("Building GL+GJ (n=%d, lc=%.3f)...\n", n_quad, l_corner)
box_gj = BI.single_dielectric_box2d(1.0, 1.0, n_quad, 0.2, l_corner, eps_in, eps_out; use_singular=true)
x_gj = solve_box(box_gj)
@printf("  N=%d, flux_err=%.3e\n", BI.num_points(box_gj),
        abs(dot(BI.all_weights(box_gj), x_gj) - exact_flux))

ep_ref = find_edge_panels(box_ref, x_ref, corner, edge_dir)
ep_gl  = find_edge_panels(box_gl,  x_gl,  corner, edge_dir)
ep_gj  = find_edge_panels(box_gj,  x_gj,  corner, edge_dir)
@printf("Edge panels: ref=%d, GL=%d, GJ=%d\n", length(ep_ref), length(ep_gl), length(ep_gj))

# ── fine uniform evaluation mesh ─────────────────────────────────────────────
n_eval = 2000
s_eval = collect(range(1e-5, 1.0, length=n_eval))

sig_ref = interpolate_density(ep_ref, s_eval, corner, edge_dir)
sig_gl  = interpolate_density(ep_gl,  s_eval, corner, edge_dir)
sig_gj  = interpolate_density(ep_gj,  s_eval, corner, edge_dir)

valid = .!isnan.(sig_ref) .& .!isnan.(sig_gl) .& .!isnan.(sig_gj)
err_gl = abs.(sig_gl .- sig_ref)
err_gj = abs.(sig_gj .- sig_ref)

@printf("Valid points: %d / %d\n", sum(valid), n_eval)
@printf("Max |σ_GL - σ_ref| = %.3e\n", maximum(err_gl[valid]))
@printf("Max |σ_GJ - σ_ref| = %.3e\n", maximum(err_gj[valid]))

# ── print error table along the full edge ────────────────────────────────────
println("\n" * "="^85)
println("Density error along full edge: corner $corner → $((corner[1]+1.0, corner[2]))")
println("="^85)
@printf("%-14s  %-14s  %-14s  %-14s  %-14s\n",
        "s", "σ_ref", "err_GL", "err_GJ", "GL/GJ")
println("-"^85)
# print at selected s values spanning the edge
s_print = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]
for sp in s_print
    idx = argmin(abs.(s_eval .- sp))
    if valid[idx]
        ratio = err_gl[idx] > 0 && err_gj[idx] > 0 ? err_gl[idx] / err_gj[idx] : NaN
        @printf("%.4e      %+.6e    %.4e      %.4e      %.1f\n",
                s_eval[idx], sig_ref[idx], err_gl[idx], err_gj[idx], ratio)
    end
end
# per-panel summary
println("\n--- Per-panel max error ---")
@printf("%-20s  %-10s  %-14s  %-14s  %-10s\n",
        "panel range", "type", "max err GL", "max err GJ", "GL/GJ")
println("-"^75)
for (i, ep) in enumerate(ep_gl)
    s0, s1 = ep.s_start, ep.s_end
    mask_p = valid .& (s_eval .>= s0) .& (s_eval .<= s1)
    if any(mask_p)
        me_gl = maximum(err_gl[mask_p])
        me_gj = maximum(err_gj[mask_p])
        label = ep.panel.is_singular ? "GJ" : "GL"
        # for GL box all panels are GL; for GJ box the corner ones are GJ
        ep_gj_panel = ep_gj[i]
        label_gj = ep_gj_panel.panel.is_singular ? "GJ" : "GL"
        ratio = me_gl > 0 && me_gj > 0 ? me_gl / me_gj : NaN
        @printf("[%.4f, %.4f]     %-10s  %.4e      %.4e      %.1f\n",
                s0, s1, label_gj, me_gl, me_gj, ratio)
    end
end
println()

# quad node positions
s_nodes_gl = quad_node_s(ep_gl, corner, edge_dir)
s_nodes_gj = quad_node_s(ep_gj, corner, edge_dir)

# power-law fit from ref near corner
mask_fit = .!isnan.(sig_ref) .& (s_eval .< 0.01) .& (abs.(sig_ref) .> 0)
ls = log.(s_eval[mask_fit]); lsig = log.(abs.(sig_ref[mask_fit]))
nf = sum(mask_fit)
b_fit = (nf * dot(ls, lsig) - sum(ls)*sum(lsig)) / (nf*dot(ls,ls) - sum(ls)^2)
C_fit = exp((sum(lsig) - b_fit*sum(ls)) / nf)
@printf("Power-law fit (ref, s<0.01): α=%.4f, C=%.4e\n", b_fit, C_fit)

# evaluate σ at GL and GJ node positions from reference (for node markers)
sig_ref_at_gl = interpolate_density(ep_ref, s_nodes_gl, corner, edge_dir)
sig_ref_at_gj = interpolate_density(ep_ref, s_nodes_gj, corner, edge_dir)
sig_gl_at_gl  = interpolate_density(ep_gl,  s_nodes_gl, corner, edge_dir)
sig_gj_at_gj  = interpolate_density(ep_gj,  s_nodes_gj, corner, edge_dir)

dir = joinpath(@__DIR__, "..", "output", "density")

# ── PLOT 1: density σ(s) log-log ────────────────────────────────────────────
fig = Figure(size=(900, 550))
ax = Axis(fig[1, 1],
    xlabel = "s  (arc-length from corner)",
    ylabel = "|σ(s)|",
    xscale = log10, yscale = log10,
    title  = "Density near corner $corner, top edge\nε_in=$eps_in, ε_out=$eps_out, γ=$(round(gamma,digits=4))")

# power-law
s_pw = collect(range(1e-5, 0.3, length=300))
lines!(ax, s_pw, C_fit .* s_pw .^ b_fit, color=:gray70, linestyle=:dash, linewidth=1.5,
       label="C·s^$(round(b_fit,digits=2))")

# interpolated curves
lines!(ax, s_eval[valid], abs.(sig_ref[valid]), color=:black, linewidth=2,
       label="ref GL (n=32, lc=0.0002)")
lines!(ax, s_eval[valid], abs.(sig_gl[valid]), color=(:blue, 0.8), linewidth=1.5,
       label="GL (n=$n_quad, lc=$l_corner)")
lines!(ax, s_eval[valid], abs.(sig_gj[valid]), color=(:red, 0.8), linewidth=1.5,
       label="GJ (n=$n_quad, lc=$l_corner)")

# node markers
scatter!(ax, s_nodes_gl, abs.(sig_gl_at_gl), color=:blue, markersize=6,
         marker=:circle)
scatter!(ax, s_nodes_gj, abs.(sig_gj_at_gj), color=:red, markersize=6,
         marker=:utriangle)

axislegend(ax, position=:rt)
save(joinpath(dir, "density_near_corner.png"), fig, px_per_unit=3)
save(joinpath(dir, "density_near_corner.pdf"), fig)
println("Saved density_near_corner.{png,pdf}")

# ── PLOT 2: pointwise error ─────────────────────────────────────────────────
fig2 = Figure(size=(900, 550))
ax2 = Axis(fig2[1, 1],
    xlabel = "s  (arc-length from corner)",
    ylabel = "|σ − σ_ref|",
    xscale = log10, yscale = log10,
    title  = "Density error vs reference")

lines!(ax2, s_eval[valid], err_gl[valid] .+ 1e-20, color=:blue, linewidth=1.5,
       label="GL (n=$n_quad, lc=$l_corner)")
lines!(ax2, s_eval[valid], err_gj[valid] .+ 1e-20, color=:red, linewidth=1.5,
       label="GJ (n=$n_quad, lc=$l_corner)")

# vertical lines at panel boundaries
for ep in ep_gl
    vlines!(ax2, [ep.s_end], color=(:blue, 0.2), linewidth=0.5)
end
for ep in ep_gj
    vlines!(ax2, [ep.s_end], color=(:red, 0.2), linewidth=0.5, linestyle=:dash)
end

axislegend(ax2, position=:rt)
save(joinpath(dir, "density_error.png"), fig2, px_per_unit=3)
save(joinpath(dir, "density_error.pdf"), fig2)
println("Saved density_error.{png,pdf}")

# ── PLOT 3: corner zoom with node positions ─────────────────────────────────
fig3 = Figure(size=(900, 550))
ax3 = Axis(fig3[1, 1],
    xlabel = "s  (arc-length from corner)",
    ylabel = "|σ(s)|",
    xscale = log10, yscale = log10,
    title  = "Corner panel zoom (s < $(2*l_corner))")

mask_z = valid .& (s_eval .< 2 * l_corner)
lines!(ax3, s_eval[mask_z], abs.(sig_ref[mask_z]), color=:black, linewidth=2, label="ref")
lines!(ax3, s_eval[mask_z], abs.(sig_gl[mask_z]),  color=:blue, linewidth=1.5, label="GL")
lines!(ax3, s_eval[mask_z], abs.(sig_gj[mask_z]),  color=:red,  linewidth=1.5, label="GJ")
lines!(ax3, s_pw, C_fit .* s_pw .^ b_fit, color=:gray70, linestyle=:dash, linewidth=1,
       label="C·s^α")

# nodes in zoom range
m_gl = s_nodes_gl .< 2 * l_corner
m_gj = s_nodes_gj .< 2 * l_corner
scatter!(ax3, s_nodes_gl[m_gl], abs.(sig_gl_at_gl[m_gl]),
         color=:blue, markersize=8, marker=:circle, label="GL nodes")
scatter!(ax3, s_nodes_gj[m_gj], abs.(sig_gj_at_gj[m_gj]),
         color=:red, markersize=8, marker=:utriangle, label="GJ nodes")

axislegend(ax3, position=:rb)
save(joinpath(dir, "density_corner_zoom.png"), fig3, px_per_unit=3)
save(joinpath(dir, "density_corner_zoom.pdf"), fig3)
println("Saved density_corner_zoom.{png,pdf}")

# ── PLOT 4: error along full edge, s ∈ [0, 1], linear x-axis ────────────────
fig4 = Figure(size=(1000, 550))
ax4 = Axis(fig4[1, 1],
    xlabel = "s  (arc-length from left corner)",
    ylabel = "|σ − σ_ref|",
    yscale = log10,
    title  = "Density error along full edge\nn=$n_quad, lc=$l_corner, ε_in=$eps_in, ε_out=$eps_out")

lines!(ax4, s_eval[valid], err_gl[valid] .+ 1e-20, color=:blue, linewidth=1.5,
       label="GL")
lines!(ax4, s_eval[valid], err_gj[valid] .+ 1e-20, color=:red, linewidth=1.5,
       label="GL+GJ")

# mark panel boundaries
for ep in ep_gl
    vlines!(ax4, [ep.s_end], color=(:gray, 0.3), linewidth=0.5)
end

axislegend(ax4, position=:ct)
save(joinpath(dir, "density_error_full_edge.png"), fig4, px_per_unit=3)
save(joinpath(dir, "density_error_full_edge.pdf"), fig4)
println("Saved density_error_full_edge.{png,pdf}")

# ── PLOT 5: per-panel max error bar chart ────────────────────────────────────
n_panels_edge = length(ep_gl)
panel_labels = String[]
max_err_gl_panels = Float64[]
max_err_gj_panels = Float64[]
panel_types = String[]

for (i, ep) in enumerate(ep_gl)
    s0, s1 = ep.s_start, ep.s_end
    mask_p = valid .& (s_eval .>= s0) .& (s_eval .<= s1)
    any(mask_p) || continue
    push!(max_err_gl_panels, maximum(err_gl[mask_p]))
    push!(max_err_gj_panels, maximum(err_gj[mask_p]))
    is_gj = ep_gj[i].panel.is_singular
    push!(panel_types, is_gj ? "GJ" : "GL")
    push!(panel_labels, @sprintf("[%.2f,%.2f]%s", s0, s1, is_gj ? "*" : ""))
end

fig4 = Figure(size=(1000, 500))
ax4 = Axis(fig4[1, 1],
    xlabel = "panel  (* = GJ corner panel)",
    ylabel = "max |σ − σ_ref| on panel",
    yscale = log10,
    xticks = (1:length(panel_labels), panel_labels),
    xticklabelrotation = pi/4,
    title  = "Per-panel density error: GL vs GL+GJ  (n=$n_quad, lc=$l_corner)")

xs_bar = collect(1:length(panel_labels))
w = 0.35
barplot!(ax4, xs_bar .- w/2, max_err_gl_panels, width=w, color=:blue, label="GL")
barplot!(ax4, xs_bar .+ w/2, max_err_gj_panels, width=w, color=:red,  label="GL+GJ")

axislegend(ax4, position=:rt)
save(joinpath(dir, "density_panel_error.png"), fig4, px_per_unit=3)
save(joinpath(dir, "density_panel_error.pdf"), fig4)
println("Saved density_panel_error.{png,pdf}")

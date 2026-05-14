#=
Figure 7 plot script. Loads fig7_data.jls and generates fig7_convergence.png.

Layout (matches prompts_fig7.md):
  (a) Relative integral error E_V(p, r) = |V_{p,r} - V_ref| / |V_ref|
      versus edge-refinement depth r, one curve per Gauss-Legendre order p.
  (b) GMRES iteration count N_iter(p, r) versus r for the same (p, r) set.
      Shared color/marker scheme with panel (a).

The reference (p_ref, r_ref) is included in the records and is filtered out
of panel (a) (its E_V is zero by construction); it is plotted in panel (b)
since the iteration count is still meaningful.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf

const datapath = joinpath(@__DIR__, "fig7_data.jls")
const data = open(deserialize, datapath, "r")

const p_list = data.sweep.p_list
const p_ref  = data.sweep.p_ref
const V_ref  = data.V_ref

# Per-p color + marker assignment, reused across both panels.
const COLORS = Dict(2 => :black, 4 => :royalblue, 6 => :crimson)
const MARKERS = Dict(2 => :circle, 4 => :rect, 6 => :diamond)

# Pull rows for a given p from the flat results vector. Sort by r for clean
# line connections. Optionally drop the reference row (E_V = 0).
function rows_for(results, p::Int; drop_ref::Bool = false)
    rs = [r for r in results if r.p == p && (!drop_ref || !get(r, :is_ref, false))]
    sort!(rs; by = x -> x.r)
    return rs
end

begin
    fig = Figure(size = (1000, 450), fontsize = 20)

    # ----- Panel (a): E_V vs r --------------------------------------------
    ax_a = Axis(fig[1, 1];
                yscale = log10,
                xlabel = L"r", ylabel = L"\mathcal{E}_V")

    for p in p_list
        rs = rows_for(data.results, p; drop_ref = true)
        isempty(rs) && continue
        xs   = Float64[r.r for r in rs]
        errs = Float64[abs(r.V - V_ref) / abs(V_ref) for r in rs]
        keep = isfinite.(errs) .& (errs .> 0)
        scatterlines!(ax_a, xs[keep], errs[keep];
                      color = COLORS[p],
                      marker = MARKERS[p], markersize = 12, linewidth = 2,
                      label = L"p = %$p")
    end
    axislegend(ax_a; position = :rt, framevisible = false)

    # ----- Panel (b): GMRES iteration count vs r --------------------------
    ax_b = Axis(fig[1, 2];
                xlabel = L"r", ylabel = L"N_{\mathrm{iter}}")

    for p in p_list
        rs = rows_for(data.results, p; drop_ref = false)
        # Deduplicate the reference (it may also appear as the (p_ref, r_ref) row
        # under is_ref = true).
        seen = Set{Int}()
        keep_rs = filter(r -> begin
            key = r.r
            key in seen ? false : (push!(seen, key); true)
        end, rs)
        sort!(keep_rs; by = x -> x.r)
        xs  = Float64[r.r for r in keep_rs]
        its = Float64[r.n_iter for r in keep_rs]
        scatterlines!(ax_b, xs, its;
                      color = COLORS[p],
                      marker = MARKERS[p], markersize = 12, linewidth = 2,
                      label = L"p = %$p")
    end
    axislegend(ax_b; position = :rb, framevisible = false)

    colgap!(fig.layout, 1, 30)

    outpath = joinpath(@__DIR__, "fig7_convergence.png")
    save(outpath, fig; px_per_unit = 4)
    @info "Saved" outpath

    fig
end

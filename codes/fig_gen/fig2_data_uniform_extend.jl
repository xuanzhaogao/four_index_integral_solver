#=
fig2_data_uniform_extend.jl

Extends the *uniform* refinement sweeps of Figure 2 so that each one reaches at
least the largest N of the corresponding adaptive sweep. The original
fig2_data.jl stopped at k = 4, which left the uniform curves spanning
E_f ∈ [1, 3e-3] while the adaptive curves spanned [4e-4, 1e-11] — the two
families barely overlapped, so the figure could not actually be read as a
comparison at fixed accuracy.

This runs incrementally: the (expensive) adaptive sweeps and the panelization
records are read from the existing fig2_data.jls and passed through untouched;
only the new uniform levels are computed.

It also adds an `Ef_abs` field to every record — the absolute error of Eq. (Ef)
of the article. The original `Ef` field (relative, i.e. divided by
fnorm = ‖f‖_{L²(Γ)}) is preserved unchanged.

Uniform N is quantized as N = 6·4^k·p², so an exact match to the adaptive max N
is impossible; we take the smallest k with N_uniform(k) >= max N_adaptive.

Usage:
  julia --project=codes/fig_gen codes/fig_gen/fig2_data_uniform_extend.jl
=#

using Serialization
using Printf

include(joinpath(@__DIR__, "fig2_common.jl"))

const datapath = joinpath(@__DIR__, "fig2_data.jls")
const data     = open(deserialize, datapath, "r")
const fnorm    = data.fnorm

@info "Loaded existing data" datapath fnorm

# Smallest k with 6*4^k*p^2 >= N_target
function k_for_N(N_target::Integer, p::Integer)
    k = 0
    while 6 * 4^k * p^2 < N_target
        k += 1
    end
    return k
end

# ---------------------------------------------------------------------------
# Consistency check: recompute the last existing uniform level with the
# refactored (absolute) E_f and confirm it reproduces the stored relative value.
# ---------------------------------------------------------------------------
let sw = first(data.sweeps)
    r = sw.uniform[end]
    ifc = uniform_box_interface(r.k, sw.p)
    E_abs = E_f_abs_on_interface(ifc, sw.p, 3 * sw.p)
    E_rel = E_abs / fnorm
    reldiff = abs(E_rel - r.Ef) / r.Ef
    @printf("check p=%d k=%d: stored Ef=%.8e  recomputed=%.8e  reldiff=%.2e\n",
            sw.p, r.k, r.Ef, E_rel, reldiff)
    reldiff < 1e-10 || error("refactored E_f does not reproduce stored value (reldiff = $reldiff)")
end

# ---------------------------------------------------------------------------
# Extend each sweep
# ---------------------------------------------------------------------------
new_sweeps = NamedTuple[]

for sw in data.sweeps
    p = sw.p
    q = 3 * p
    N_adapt_max = maximum(r.N for r in sw.adaptive)
    k_have   = maximum(r.k for r in sw.uniform)
    k_target = k_for_N(N_adapt_max, p)

    @info "=== p = $p ===" N_adapt_max k_have k_target

    # pass adaptive records through, adding the absolute error
    adaptive_new = [merge(r, (Ef_abs = r.Ef * fnorm,)) for r in sw.adaptive]

    # existing uniform records: same treatment
    uniform_new = [merge(r, (Ef_abs = r.Ef * fnorm,)) for r in sw.uniform]

    for k in (k_have + 1):k_target
        t0 = time()
        ifc = uniform_box_interface(k, p)
        np  = length(ifc.panels)
        E_abs = E_f_abs_on_interface(ifc, p, q)
        N = np * p^2
        push!(uniform_new, (k = k, N = N, Ef = E_abs / fnorm,
                            npanels = np, Ef_abs = E_abs))
        @printf("  k=%d  panels=%-9d N=%-10d  Ef_abs=%.6e  Ef_rel=%.6e  (%.1fs)\n",
                k, np, N, E_abs, E_abs / fnorm, time() - t0)
        flush(stdout)
        ifc = nothing
        GC.gc()
    end

    sort!(uniform_new, by = r -> r.k)
    push!(new_sweeps, (p = p, adaptive = adaptive_new, uniform = uniform_new))
end

# ---------------------------------------------------------------------------
# Save (backing up the original once)
# ---------------------------------------------------------------------------
out = merge(data, (sweeps = new_sweeps,))

const backup = datapath * ".bak"
if !isfile(backup)
    cp(datapath, backup)
    @info "Backed up original" backup
end

open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved extended data" datapath bytes=stat(datapath).size

for sw in new_sweeps
    println("\n=== p = $(sw.p) ===")
    println("  ADAPTIVE:  N range $(minimum(r.N for r in sw.adaptive)) .. $(maximum(r.N for r in sw.adaptive))")
    println("  UNIFORM :  N range $(minimum(r.N for r in sw.uniform)) .. $(maximum(r.N for r in sw.uniform))")
    println("  k    N            Ef_abs")
    for r in sw.uniform
        @printf("  %-3d  %-11d  %.6e\n", r.k, r.N, r.Ef_abs)
    end
end

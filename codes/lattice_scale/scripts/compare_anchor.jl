# scripts/compare_anchor.jl — cross-check in-memory vs batched four-index integrals.
#
# Compares two independent computation paths for the demo_2x2 campaign:
#   1. In-memory: four_index_integrals(toml) — runs all solves in a single process,
#      no disk I/O, all phases in RAM.
#   2. Batched (file-based): reads the V_full.tsv produced by the full
#      prepare→solve→consolidate→eval→assemble pipeline.
#
# Run after the full pipeline has completed:
#   julia --project scripts/compare_anchor.jl [--toml campaigns/demo_2x2.toml]
#
# If V_full.tsv is not present, only the in-memory result is printed (no comparison).
using BoundaryIntegral, Printf, LinearAlgebra

function parse_v_tsv(path::AbstractString)
    # Header: i j k l V
    # Returns pair_ids (unique (i,j) pairs, sorted by first appearance as row) and V matrix.
    rows = Tuple{Int,Int}[]
    cols = Tuple{Int,Int}[]
    entries = Dict{Tuple{NTuple{2,Int}, NTuple{2,Int}}, Float64}()
    for (n, line) in enumerate(eachline(path))
        n == 1 && continue   # skip header
        f = split(line, '\t')
        length(f) == 5 || continue
        i, j, k, l = parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]), parse(Int, f[4])
        v = parse(Float64, f[5])
        entries[((i, j), (k, l))] = v
        (i, j) in rows || push!(rows, (i, j))
        (k, l) in cols || push!(cols, (k, l))
    end
    rows == cols || @warn "compare_anchor: row/col pair sets differ"
    pair_ids = rows
    n = length(pair_ids)
    idx = Dict(p => i for (i, p) in enumerate(pair_ids))
    V = fill(NaN, n, n)
    for ((ij, kl), v) in entries
        r = get(idx, ij, 0); c = get(idx, kl, 0)
        r > 0 && c > 0 && (V[r, c] = v)
    end
    return pair_ids, V
end

function main()
    toml = joinpath(@__DIR__, "..", "campaigns", "demo_2x2.toml")
    for (i, arg) in enumerate(ARGS)
        arg == "--toml" && i < length(ARGS) && (toml = ARGS[i+1])
    end
    toml = abspath(toml)

    @info "compare_anchor: in-memory solve" toml
    t0 = time()
    res = four_index_integrals(toml)
    t_mem = time() - t0
    @info "in-memory done" t=round(t_mem; digits=1) n_pairs=length(res.pair_ids) size_V=size(res.V)

    # Symmetry of the in-memory result
    scale_mem = maximum(abs.(res.V))
    asym_mem  = maximum(abs.(res.V .- transpose(res.V))) / scale_mem
    @printf("In-memory max|V| = %.6e   max rel asymmetry = %.3e\n", scale_mem, asym_mem)

    # Print in-memory matrix (small campaigns only)
    if length(res.pair_ids) <= 30
        println("\nIn-memory V matrix (pair_ids = $(res.pair_ids)):")
        for r in axes(res.V, 1)
            for c in axes(res.V, 2)
                @printf("  %10.4e", res.V[r, c])
            end
            println()
        end
    end

    # Compare against batched V_full.tsv if present
    c = load_campaign(toml)
    tsv_path = joinpath(c.root, "V_full.tsv")
    if !isfile(tsv_path)
        @info "V_full.tsv not found — skipping batched comparison" tsv_path
        @info "Run the full pipeline first: prepare → solve → consolidate → eval → assemble"
        return
    end

    @info "Loading batched V_full.tsv" tsv_path
    bat_pairs, V_bat = parse_v_tsv(tsv_path)

    # Align pair ordering
    if res.pair_ids != bat_pairs
        @warn "Pair id ordering differs between in-memory and batched; realigning."
    end
    mem_idx = Dict(p => i for (i, p) in enumerate(res.pair_ids))
    bat_idx = Dict(p => i for (i, p) in enumerate(bat_pairs))
    common = intersect(res.pair_ids, bat_pairs)
    length(common) == length(res.pair_ids) == length(bat_pairs) ||
        @warn "Pair sets not identical" n_mem=length(res.pair_ids) n_bat=length(bat_pairs) n_common=length(common)

    n = length(common)
    V_m = Matrix{Float64}(undef, n, n)
    V_b = Matrix{Float64}(undef, n, n)
    for (r, pr) in enumerate(common), (cc, pc) in enumerate(common)
        V_m[r, cc] = res.V[mem_idx[pr], mem_idx[pc]]
        V_b[r, cc] = V_bat[bat_idx[pr], bat_idx[pc]]
    end

    scale = maximum(abs.(V_b))
    diff  = maximum(abs.(V_m .- V_b))
    rel   = diff / max(scale, eps())
    asym_bat = maximum(abs.(V_b .- transpose(V_b))) / scale

    println("\n--- Comparison: in-memory vs batched ---")
    @printf("max|V_batched|        = %.6e\n", scale)
    @printf("max|V_mem - V_bat|    = %.6e\n", diff)
    @printf("max rel diff          = %.3e   (expect < ~1e-2 at production tols)\n", rel)
    @printf("batched rel asymmetry = %.3e\n", asym_bat)
    @printf("in-mem rel asymmetry  = %.3e\n", asym_mem)
end

main()

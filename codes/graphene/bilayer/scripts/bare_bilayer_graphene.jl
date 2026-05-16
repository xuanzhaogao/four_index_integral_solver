using CSV
using DataFrames
using Printf
using Statistics

include(joinpath(@__DIR__, "..", "..", "monolayer", "src", "MonolayerOrbitalLoader.jl"))
include(joinpath(@__DIR__, "..", "..", "monolayer", "src", "MonolayerBareIntegrals.jl"))
include(joinpath(@__DIR__, "..", "src", "BilayerBareIntegrals.jl"))
using .MonolayerOrbitalLoader
using .MonolayerBareIntegrals
using .BilayerBareIntegrals

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/bilayer/k_161601_nb_288_c_15"
const ORBITALS = (
    joinpath(REF_DIR, "graphene_00001.xsf"),
    joinpath(REF_DIR, "graphene_00002.xsf"),
    joinpath(REF_DIR, "graphene_00003.xsf"),
    joinpath(REF_DIR, "graphene_00004.xsf"),
)
const SOURCE_TOL = 1e-3
const VOLUME_TOL = 1e-3
const MIRROR_PAD_LEVEL = 0

const OUT_CSV = joinpath(@__DIR__, "..", "data", "bare_bilayer_graphene.csv")

# Bilayer orbital geometry (Wannier centres, Å):
#   1: (0.0,    0.0,      7.5016)  layer 1 A
#   2: (0.0,    1.4232,   7.4968)  layer 1 B
#   3: (0.0,    1.4232,  10.8532)  layer 2 (above B1, dimer site)
#   4: (1.2325, 0.7116,  10.8484)  layer 2 (hollow over layer 1)
# In-plane (same-layer) pairs: (1,2), (3,4)
# Inter-layer pairs:           (1,3), (1,4), (2,3), (2,4)
const IN_PLANE_PAIRS = Set([(1, 2), (3, 4)])
const INTER_LAYER_PAIRS = Set([(1, 3), (1, 4), (2, 3), (2, 4)])

pair_class(i, j) = (i, j) in IN_PLANE_PAIRS ? "in_plane" :
                   (i, j) in INTER_LAYER_PAIRS ? "inter_layer" :
                   i == j ? "onsite" : "unknown"

annotate_rows(rows) = [merge(r, (; pair_class = pair_class(r.i, r.j))) for r in rows]

function print_avg_row(label::AbstractString, val::Real, ref::Real)
    @printf("%-20s %10.4f   ref = %10.4f   rel.err = %+7.3f %%\n",
            label, val, ref, 100 * (val - ref) / ref)
end

function summarize(table::DataFrame)
    println()
    println("=== Per-channel orbital averages (eV) ===")
    onsite_avg = mean(table[table.channel .== "onsite",  :u_ev])
    inter_avg  = mean(table[table.channel .== "nn",      :u_ev])
    sf_avg     = mean(table[table.channel .== "hund_sf", :u_ev])
    ph_avg     = mean(table[table.channel .== "hund_ph", :u_ev])

    ref = BilayerBareIntegrals.BILAYER_REF_AVG_EV
    print_avg_row("intra-orbital",      onsite_avg, ref.intra)
    print_avg_row("inter-orbital",      inter_avg,  ref.inter)
    print_avg_row("Hund's (spin-flip)", sf_avg,     ref.hund_sf)
    print_avg_row("Hund's (pair-hop) ", ph_avg,     ref.hund_ph)

    println()
    println("=== In-plane vs inter-layer breakdown (eV) ===")
    for channel in ("nn", "hund_sf", "hund_ph")
        ip = table[(table.channel .== channel) .& (table.pair_class .== "in_plane"),    :u_ev]
        il = table[(table.channel .== channel) .& (table.pair_class .== "inter_layer"), :u_ev]
        @printf("%-9s  in-plane avg = %.4f (n=%d)   inter-layer avg = %.4f (n=%d)\n",
                channel, mean(ip), length(ip), mean(il), length(il))
    end
end

function main()
    println("--- bilayer bare-Coulomb sweep (4 orbitals, 10 pairs) ---")
    println("REF_DIR = ", REF_DIR)
    @printf("source_tol = %.1e, volume_tol = %.1e, mirror_pad_level = %d\n",
            SOURCE_TOL, VOLUME_TOL, MIRROR_PAD_LEVEL)
    println()

    raw_rows = compute_all_bilayer_bare_channels(;
        orbitals         = ORBITALS,
        source_tol       = SOURCE_TOL,
        volume_tol       = VOLUME_TOL,
        mirror_pad_level = MIRROR_PAD_LEVEL,
    )
    rows = annotate_rows(raw_rows)

    mkpath(dirname(OUT_CSV))
    table = DataFrame(rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    show(table; allrows = true, allcols = true)
    summarize(table)
end

main()

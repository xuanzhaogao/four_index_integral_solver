# Sweep bare 4-channel integrals across three monolayer Wannier datasets
# (k_161601, k_252501, k_323201) and test Malte's Madelung-correction hypothesis.

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
include(joinpath(@__DIR__, "..", "src", "MonolayerBareIntegrals.jl"))
using .MonolayerOrbitalLoader
using .MonolayerBareIntegrals

"""
    parse_coqui_bare(path)

Parse the "bare interactions (orbital-average)" block of a CoQui
`_coqui_thc_crpa.out` file. Returns a Dict mapping the four channel symbols
to their real-part values in eV.

Stops parsing at the first `static screened` line. Throws if any of the four
expected lines are missing.
"""
function parse_coqui_bare(path::AbstractString)
    isfile(path) || throw(ArgumentError("CoQui file not found: $path"))
    result = Dict{Symbol, Float64}()
    in_block = false
    for line in eachline(path)
        if occursin("bare interactions (orbital-average)", line)
            in_block = true
            continue
        end
        in_block || continue
        if occursin("static screened", line)
            break
        end
        m = match(r"-\s*(intra-orbital|inter-orbital|Hund's coupling \(spin-flip\)|Hund's coupling \(pair-hopping\))\s*=\s*\(\s*([-+0-9.eE]+)\s*,", line)
        m === nothing && continue
        label, val_str = m.captures[1], m.captures[2]
        key = label == "intra-orbital"                       ? :onsite   :
              label == "inter-orbital"                       ? :nn       :
              label == "Hund's coupling (spin-flip)"         ? :hund_sf  :
              label == "Hund's coupling (pair-hopping)"      ? :hund_ph  :
              error("unreachable: $label")
        result[key] = parse(Float64, val_str)
    end
    for k in (:onsite, :nn, :hund_sf, :hund_ph)
        haskey(result, k) || error("CoQui parser failed to find $k in $path")
    end
    return result
end

const REF_ROOT = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer"
const SOURCE_TOL = 1e-3
const VOLUME_TOL = 1e-3

# Madelung corrections supplied by Malte (2026-05 email). Defined for the
# onsite channel only — the long-wavelength correction does not apply to the
# short-range nn / Hund's channels.
const MADELUNG_REF_EV = Dict(
    "k_161601_nb_144_c_15" =>  1.041,
    "k_252501_nb_144_c_15" =>  0.1842,
    "k_323201_nb_144_c_15" => -0.290,
)

const DATASETS = (
    "k_161601_nb_144_c_15",
    "k_252501_nb_144_c_15",
    "k_323201_nb_144_c_15",
)

const OUT_CSV = joinpath(@__DIR__, "..", "data", "bare_monolayer_kmesh_sweep.csv")

"""
    sweep_dataset(kmesh)

Run the 4-channel bare integral on a single dataset directory and combine
with parsed CoQui references. Returns a vector of NamedTuples (one per channel)
ready for assembly into the output CSV.
"""
function sweep_dataset(kmesh::AbstractString)
    dir = joinpath(REF_ROOT, kmesh)
    xsf_1 = joinpath(dir, "graphene_00001.xsf")
    xsf_2 = joinpath(dir, "graphene_00002.xsf")
    coqui = joinpath(dir, "_coqui_thc_crpa.out")
    isfile(xsf_1) || error("missing XSF 1 in $dir")
    isfile(xsf_2) || error("missing XSF 2 in $dir")
    isfile(coqui) || error("missing CoQui output in $dir")

    println("--- $kmesh ---")
    # mirror_pad_level fixed at 0: the 2026-04-23 report showed pad>=1 is
    # inconclusive on these XSFs (orbital fills the supercell, so reflections
    # produce coherent replicas rather than a vacuum-padded tail).
    rows = compute_all_bare_channels(;
        orbital_1 = xsf_1, orbital_2 = xsf_2,
        source_tol = SOURCE_TOL, volume_tol = VOLUME_TOL,
        mirror_pad_level = 0,
    )
    coqui_vals = parse_coqui_bare(coqui)
    madelung_onsite = get(MADELUNG_REF_EV, kmesh, missing)

    out = NamedTuple[]
    for r in rows
        coqui_v = coqui_vals[r.channel]
        diff_ev = coqui_v - r.u_ev
        madelung_ref = r.channel === :onsite ? madelung_onsite : missing
        madelung_residual = r.channel === :onsite ? diff_ev - madelung_ref : missing
        push!(out, (
            kmesh = kmesh,
            channel = String(r.channel),
            our_u_ev = r.u_ev,
            coqui_u_ev = coqui_v,
            diff_ev = diff_ev,
            madelung_ref_ev = madelung_ref,
            madelung_residual_ev = madelung_residual,
            our_u_raw = r.u_raw,
            Na = r.Na,
            Nb = r.Nb,
            n_source_points = r.n_source_points,
            n_target_points = r.n_target_points,
            tkm_kmax = r.tkm_kmax,
            source_tol = r.source_tol,
            volume_tol = r.volume_tol,
            mirror_pad_level = r.mirror_pad_level,
            shift_x = r.shift_x,
            shift_y = r.shift_y,
            shift_z = r.shift_z,
        ))
    end
    return out
end

function main()
    all_rows = NamedTuple[]
    for kmesh in DATASETS
        rows = sweep_dataset(kmesh)
        append!(all_rows, rows)
        # Show per-dataset summary table immediately for easy mid-run progress
        df_local = DataFrame(rows)
        show(df_local[:, [:kmesh, :channel, :our_u_ev, :coqui_u_ev, :diff_ev, :madelung_ref_ev, :madelung_residual_ev]];
             allrows = true, allcols = true)
        println()
    end
    mkpath(dirname(OUT_CSV))
    table = DataFrame(all_rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    show(table[:, [:kmesh, :channel, :our_u_ev, :coqui_u_ev, :diff_ev, :madelung_ref_ev, :madelung_residual_ev]];
         allrows = true, allcols = true)
    println()
end

main()

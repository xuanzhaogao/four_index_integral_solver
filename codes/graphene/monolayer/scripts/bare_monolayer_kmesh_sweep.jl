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

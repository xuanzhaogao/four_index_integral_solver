module ScreenedHubbardComparison

using CSV
using DataFrames

export DEFAULT_SCREENED_HUBBARD_CSV,
    GRAPHENE_CRPA_EV,
    PAIR_ORDER,
    load_softmix_hubbard_results,
    pair_curve_data

const DEFAULT_SCREENED_HUBBARD_CSV = normpath(joinpath(@__DIR__, "..", "data", "screened_hubbard_graphene.csv"))
const PAIR_ORDER = ["U_00", "U_01", "U_02", "U_03"]
const GRAPHENE_CRPA_EV = Dict(
    "U_00" => 9.3,
    "U_01" => 5.5,
    "U_02" => 4.1,
    "U_03" => 3.6,
)

function load_softmix_hubbard_results(path::AbstractString = DEFAULT_SCREENED_HUBBARD_CSV)
    df = CSV.read(path, DataFrame)
    return df[df.mode .== "SoftMixInversePermittivity", :]
end

function pair_curve_data(df::DataFrame; value_col::Symbol = :u_total_ev)
    curves = NamedTuple[]
    for pair in PAIR_ORDER
        rows = df[df.pair .== pair, [:bandwidth, value_col]]
        nrow(rows) == 0 && continue
        sort!(rows, :bandwidth)
        push!(
            curves,
            (
                pair = pair,
                bandwidths = collect(rows.bandwidth),
                values = Float64.(rows[!, value_col]),
                crpa_ev = GRAPHENE_CRPA_EV[pair],
            ),
        )
    end
    return curves
end

end

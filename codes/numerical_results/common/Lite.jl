# Harness-free helpers for analysis/plot scripts. Deliberately does NOT
# import BoundaryIntegral: the saved run data are plain NamedTuples of
# scalars/vectors, and BI is under active development (in-flight edits
# repeatedly break its precompile), so post-processing must not depend on it.
module Lite

using Serialization

export load_ref, save_ref, append_csv_row

load_ref(path::AbstractString) = deserialize(path)
save_ref(path::AbstractString, x) = (mkpath(dirname(path)); serialize(path, x); x)

function append_csv_row(path::AbstractString, row::NamedTuple)
    mkpath(dirname(path))
    newfile = !isfile(path)
    open(path, "a") do io
        newfile && println(io, join(string.(keys(row)), ","))
        println(io, join([sprint(print, v) for v in values(row)], ","))
    end
end

end # module

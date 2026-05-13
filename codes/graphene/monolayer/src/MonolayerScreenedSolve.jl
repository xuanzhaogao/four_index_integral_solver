module MonolayerScreenedSolve

export parse_coqui_loc

"""
    parse_coqui_loc(path)

Parse a CoQui `_coqui_crpa_loc.out` file's `i j k l v_ijkl U_ijkl` table and
return a Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}
for the four distinct channels:

- (0,0,0,0) → :onsite
- (0,0,1,1) → :nn
- (0,1,0,1) → :hund_ph
- (0,1,1,0) → :hund_sf
"""
function parse_coqui_loc(path::AbstractString)
    isfile(path) || throw(ArgumentError("CoQui loc file not found: $path"))
    targets = Dict{NTuple{4, Int}, Symbol}(
        (0, 0, 0, 0) => :onsite,
        (0, 0, 1, 1) => :nn,
        (0, 1, 0, 1) => :hund_ph,
        (0, 1, 1, 0) => :hund_sf,
    )
    result = Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}()
    pat = r"^\s*(\d)\s+(\d)\s+(\d)\s+(\d)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
    for line in eachline(path)
        m = match(pat, line)
        m === nothing && continue
        ijkl = (parse(Int, m.captures[1]), parse(Int, m.captures[2]),
                parse(Int, m.captures[3]), parse(Int, m.captures[4]))
        sym = get(targets, ijkl, nothing)
        sym === nothing && continue
        v = parse(Float64, m.captures[5])
        U = parse(Float64, m.captures[6])
        result[sym] = (v_ijkl = v, U_ijkl = U)
    end
    for k in (:onsite, :nn, :hund_ph, :hund_sf)
        haskey(result, k) || error("parse_coqui_loc: missing channel $k in $path")
    end
    return result
end

end # module

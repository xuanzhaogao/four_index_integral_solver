function normalize_tolerances(tolerances)
    values = sort!(collect(Float64.(tolerances)); rev = true)
    isempty(values) && throw(ArgumentError("tolerances must not be empty"))
    all(0.0 .< values .< 1.0) || throw(ArgumentError("tolerances must lie in (0, 1)"))
    return values
end

function annotate_reference_errors(pair::AbstractString, samples)
    isempty(samples) && throw(ArgumentError("samples must not be empty"))

    reference = samples[end]
    raw_denominator = abs(reference.U_raw)
    ev_denominator = abs(reference.U_ev)

    return [merge(
        (pair = pair,),
        row,
        (
            rel_err_raw = iszero(raw_denominator) ? abs(row.U_raw - reference.U_raw) : abs(row.U_raw - reference.U_raw) / raw_denominator,
            rel_err_ev = iszero(ev_denominator) ? abs(row.U_ev - reference.U_ev) : abs(row.U_ev - reference.U_ev) / ev_denominator,
        ),
    ) for row in samples]
end

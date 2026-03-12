function annotate_reference_errors(pair::AbstractString, samples)
    isempty(samples) && throw(ArgumentError("samples must not be empty"))

    reference = samples[end]
    raw_denominator = abs(reference.U_raw)
    ev_denominator = abs(reference.U_ev)

    return [(
        pair = pair,
        tol = row.tol,
        U_raw = row.U_raw,
        U_ev = row.U_ev,
        rel_err_raw = iszero(raw_denominator) ? abs(row.U_raw - reference.U_raw) : abs(row.U_raw - reference.U_raw) / raw_denominator,
        rel_err_ev = iszero(ev_denominator) ? abs(row.U_ev - reference.U_ev) : abs(row.U_ev - reference.U_ev) / ev_denominator,
    ) for row in samples]
end

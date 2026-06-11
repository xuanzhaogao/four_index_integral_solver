# driver.jl — thin campaign driver on BoundaryIntegral.jl campaign API.
#
# Usage:  julia --project -t 8 driver.jl <campaign.toml> <phase> [--only ID | --workers N]
#   phases: prepare | solve | consolidate | eval | assemble | status
#
# Precompiles BoundaryIntegral (incl. the Distributed extension) on the HEAD process
# BEFORE spawning workers — prevents the worker-load hang caused by parallel precompile
# races over GPFS.  Workers reload the campaign from the toml path (cheap, cached).
using Pkg
Pkg.precompile()                          # must be BEFORE loading Distributed/SlurmClusterManager

using BoundaryIntegral, Distributed, SlurmClusterManager

# Julia 1.12 fully buffers stdout/stderr when redirected to a file (sbatch logs),
# so nothing appears until process exit — flush every 2 s to keep the log live.
Timer(_ -> (flush(stdout); flush(stderr)), 2; interval = 2)

function _parse_args(args)
    length(args) >= 2 ||
        error("usage: driver.jl <campaign.toml> <phase> [--only ID] [--workers N]")
    toml  = args[1]
    phase = Symbol(args[2])
    only_id = nothing
    nworkers = 0
    i = 3
    while i <= length(args)
        if args[i] == "--only"
            i + 1 <= length(args) || error("--only requires a value")
            only_id = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--workers"
            i + 1 <= length(args) || error("--workers requires a value")
            nworkers = parse(Int, args[i+1]); i += 2
        else
            error("unknown argument $(args[i])")
        end
    end
    return toml, phase, only_id, nworkers
end

function main()
    toml, phase, only_id, nworkers = _parse_args(ARGS)
    c = load_campaign(toml)

    if phase === :prepare
        prepare(c)
    elseif phase === :consolidate
        consolidate(c)
    elseif phase === :assemble
        assemble_v(c)
    elseif phase === :status
        for ph in (:solve, :eval)
            p = pending_batches(c, ph)
            println("$ph: $(length(p)) pending$(isempty(p) ? "" : "  (ids $(first(p, min(10, length(p))))…)")")
        end
    elseif phase in (:solve, :eval)
        if only_id !== nothing
            (phase === :solve ? solve_batch : eval_batch)(c, only_id)
        else
            run_phase(c, phase; workers = nworkers)
        end
    else
        error("unknown phase $phase")
    end
end

main()

# driver.jl — campaign driver (spec §4).
#
# Usage:  julia --project -t 8 driver.jl <campaign.toml> <phase> [--only ID] [--workers N]
#   phases: prepare | solve | consolidate | eval | assemble | status
#
# Worker topology: ONE worker per allocated Slurm task (= one per node with
# --ntasks-per-node=1); each task uses the whole node's cores via OMP threads.
# Restart semantics: status lives on disk; resubmitting after any crash redoes at most
# the in-flight batches.
using Distributed
using Pkg
Pkg.precompile()                      # precompile BEFORE spawning workers (GPFS cache race)
using CampaignLib
using ClusterManagers

function _parse_args(args)
    length(args) >= 2 || error("usage: driver.jl <campaign.toml> <phase> [--only ID] [--workers N]")
    toml, phase = args[1], Symbol(args[2])
    only_id = nothing; nworkers_local = 0
    i = 3
    while i <= length(args)
        if args[i] == "--only"
            i + 1 <= length(args) || error("--only requires a value")
            only_id = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--workers"
            i + 1 <= length(args) || error("--workers requires a value")
            nworkers_local = parse(Int, args[i+1]); i += 2
        else
            error("unknown arg $(args[i])")
        end
    end
    return toml, phase, only_id, nworkers_local
end

function _setup_workers(nworkers_local::Int)
    proj = dirname(Base.active_project())
    exe = "--project=$proj"
    glue = get(ENV, "JULIA_GLUE_THREADS", "8")           # Threads.@threads glue loops
    if haskey(ENV, "SLURM_JOB_ID") && parse(Int, get(ENV, "SLURM_NTASKS", "1")) > 1
        np = parse(Int, ENV["SLURM_NTASKS"])
        @info "spawning $np Slurm workers (one per task)"
        addprocs(SlurmManager(np); exeflags = `$exe -t $glue`)
    elseif nworkers_local > 0
        @info "spawning $nworkers_local local workers"
        addprocs(nworkers_local; exeflags = `$exe -t $glue`)
    else
        @info "no workers: running tasks inline (pilot mode)"
    end
end

function main()
    toml, phase, only_id, nworkers_local = _parse_args(ARGS)
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
            println("$ph: $(length(p)) pending  $(isempty(p) ? "" : "(ids $(first(p, min(10, length(p))))…)")")
        end
    elseif phase in (:solve, :eval)
        runner = phase === :solve ? solve_batch : eval_batch
        if only_id !== nothing
            runner(c, only_id)                            # inline pilot, full timing visible
            return
        end
        _setup_workers(nworkers_local)
        @everywhere eval(:(using CampaignLib))
        pending = pending_batches(c, phase)
        @info "$(phase): $(length(pending)) pending batches on $(nworkers()) workers"
        results = pmap(pending; retry_delays = [30.0], on_error = e -> e) do id
            try
                runner(load_campaign(toml), id)           # campaign reloaded per worker (cached)
                (id, :ok)
            catch err
                msg = sprint(showerror, err, catch_backtrace())
                cc = load_campaign(toml)
                mkpath(logs_dir(cc))
                write(joinpath(logs_dir(cc), "$(phase)_batch_$(lpad(id, 4, '0')).err"), msg)
                rethrow()
            end
        end
        ok = count(r -> r isa Tuple && r[2] === :ok, results)
        @info "$(phase) finished" ok failed=length(results)-ok
        still = pending_batches(c, phase)
        isempty(still) || @warn "still pending (resubmit to retry)" still
    else
        error("unknown phase $phase")
    end
end

main()

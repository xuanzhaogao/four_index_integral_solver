# A/B for the eval-phase regression: the campaign's per-FMM-application cost went 8.8 s
# (Aug 30, job 6965032) to 61 s (kb31, job 7006478) at an IDENTICAL problem size --
# 1.37M refined sources against the full 7.3M-point target set, K sources per batch.
#
# Already excluded, each by measurement rather than argument:
#   max_order        -- reaches solve_batch_core only; eval never passes it
#   the triangle     -- the regression is present at FULL target sets (kb31 batch 1)
#   target scatter   -- `keep` is sorted, `tgt` is a contiguous materialized copy
#   hcubature near   -- "num of hcub calculations" never logged: the list is empty
#   the BI kernels   -- Sec. 5.2 on current BI runs at 10.4 s, close to Aug 30's 8.8 s
#
# So it is in evaluate_batch_potential as the campaign calls it, between BI 2aac4ac (Aug 5,
# what the Aug 30 run used) and e95f8af. This script times ONE batch's evaluation so the two
# BI checkouts can be compared directly. Run it once per checkout via --project.
#
#   BATCH_ID=85 CAMP=lattice_conv_l3_eps2.4_kb31 julia --project -t 96 scripts/ab_eval_regression.jl
#
# Batch 85 is the cheap end (n_targets_used 443,892, 143 s in production); batch 1 is the full
# target set. Time BOTH: if only the full-target case regresses, the cause scales with target
# count, which points somewhere different than if both do.
using BoundaryIntegral, Serialization, Printf
const BI = BoundaryIntegral
const CEPH = "/mnt/ceph/users/xgao1/four_index"
const CAMP = get(ENV, "CAMP", "lattice_conv_l3_eps2.4_kb31")
const BID  = parse(Int, get(ENV, "BATCH_ID", "85"))

c       = load_campaign(joinpath(@__DIR__, "..", "campaigns", CAMP * ".toml"))
targets = open(deserialize, BI.targets_path(c))
store   = open(deserialize, BI.rho_store_path(c))
br      = BI.load_batch_result(BI.batch_path(c, BID))
dg      = BI.load_templates!(c)[1][2]

println("BI at ", read(`git -C $(pkgdir(BI)) rev-parse --short HEAD`, String) |> strip,
        "   campaign $(CAMP)  batch $(BID)  K = $(length(br.pair_ids))")
println("threads = ", Threads.nthreads(), "  OPENBLAS = ", get(ENV, "OPENBLAS_NUM_THREADS", "unset"))

for tri in (true, false)
    t0 = time()
    _, _, rows, ntgt = BI.eval_batch_core(br, targets, store, dg, c; triangle = tri)
    dt = time() - t0
    @printf("triangle=%-5s  rows %5d  n_targets %8d  t %8.1f s  -> %.1f s per source, %.2f us per source-point\n",
            tri, length(rows), ntgt, dt, dt / length(br.pair_ids),
            1e6 * dt / (length(br.pair_ids) * ntgt))
end

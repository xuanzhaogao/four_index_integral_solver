# Why did the campaign's eval phase go from 39 min (Aug 30, job 6965032) to 110 min (kb31,
# job 7006478)?  Normalized by the work actually done -- sources x target POINTS, read from each
# V file's n_targets_used -- cost per source-point went 1.21 us to 5.91 us, i.e. per source
# against the full target set, 8.23 s (Aug 30 fit over K = 1..17) to 61.4 s (kb31, K = 31).
#
# Excluded by measurement, each rather than by argument:
#   max_order 8->64   reaches solve_batch_core only; eval never passes it, and the refined
#                     source count is unchanged (1.373M Aug 30 vs 1.371M kb31)
#   the triangle      the gap is at FULL target sets, where the triangle does nothing
#   target scatter    `keep` is sorted, `tgt` a contiguous materialized copy
#   hcubature near    "num of hcub calculations" absent from all logs: the list is empty, so
#                     the ONLY functions 471ce6e/e95f8af touched in this path never run
#   the BI version    reflog: HEAD was 2aac4ac from Aug 5 to Sep 8, and 2aac4ac..e95f8af is
#                     four files, none of them executing code in this path
#   FINUFFT cap       neither campaign jobscript sets TKM3D_FINUFFT_NTHREADS; same both eras
#
# What is NOT excluded, because nobody has measured it: the K range 17 < K < 31. The Aug 30
# slope is fitted over K = 1..17 and extrapolated; the kb31/k46 rate rests on two points
# (61.4 s at K = 31, 61.8 s at K = 46). A jump anywhere in that gap would explain everything
# with no regression at all. This sweeps it on ONE batch at the full target set.
#
# Second variable: OPENBLAS_NUM_THREADS was 96 on Aug 30 (run_conv_l3_eps2.4.sbatch) and 1 for
# kb31 (run_all.sbatch). Sec. 5.2 measured 1 as much FASTER for its evaluation, but that was a
# different code path, so it is re-tested here in-process at the top K.
#
#   CAMP=lattice_conv_l3_eps2.4_kb31 BATCH_ID=1 KS=4,8,13,17,24,31 \
#     julia --project -t 96 scripts/ab_eval_regression.jl
using BoundaryIntegral, Serialization, Printf, LinearAlgebra
const BI = BoundaryIntegral

const CAMP = get(ENV, "CAMP", "lattice_conv_l3_eps2.4_kb31")
const BID  = parse(Int, get(ENV, "BATCH_ID", "1"))
const KS   = parse.(Int, split(get(ENV, "KS", "4,8,13,17,24,31"), ','))

c       = load_campaign(joinpath(@__DIR__, "..", "campaigns", CAMP * ".toml"))
targets = open(deserialize, BI.targets_path(c))
br      = BI.load_batch_result(BI.batch_path(c, BID))
dg      = BI.load_templates!(c)[1][2]

pos = BI.grid_positions(dg, br.gidx)
At, Bt, Ct = BI.true_cell_vectors(dg)
lb = ((At[1]/dg.nx, At[2]/dg.nx, At[3]/dg.nx), (Bt[1]/dg.ny, Bt[2]/dg.ny, Bt[3]/dg.ny),
      (Ct[1]/dg.nz, Ct[2]/dg.nz, Ct[3]/dg.nz))
mksrc(K) = [BI.VolumeSource(copy(pos), copy(br.weights), br.densities[:, k]; lattice_basis = lb)
            for k in 1:K]

tgt   = targets.positions                      # FULL target set, as Aug 30's batches used
ntgt  = size(tgt, 2)
Kmax  = length(br.pair_ids)

println("BI ", strip(read(`git -C $(pkgdir(BI)) rev-parse --short HEAD`, String)),
        "  campaign $(CAMP) batch $(BID)  K_avail = $Kmax  n_targets = $ntgt")
println("julia threads = ", Threads.nthreads(), "  BLAS = ", BLAS.get_num_threads(),
        "  OPENBLAS_NUM_THREADS = ", get(ENV, "OPENBLAS_NUM_THREADS", "unset"))
println("\nreference: Aug 30 fit t = 7.1 + 8.23*K  ->  8.23 s per source at the full target set\n")

# Sigma carries one column per source, so it must be sliced with the sources -- otherwise
# evaluate_batch_potential raises "Sigma columns != number of sources".
run1(K) = begin
    srcs = mksrc(K)
    sig  = br.sigma[:, 1:K]
    t0 = time()
    BI.evaluate_batch_potential(br.interface, sig, srcs, tgt;
        lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"], c_pad = c.c_pad,
        screen_boxes = c.boxes, screen_epses = c.epses, screen_eps_out = c.eps_out)
    time() - t0
end

println("K sweep at the full target set:")
@printf("%4s  %10s  %12s  %14s  %s\n", "K", "t (s)", "s/source", "us/src-point", "vs Aug 30")
for K in filter(<=(Kmax), KS)
    dt = run1(K)
    @printf("%4d  %10.1f  %12.2f  %14.3f  %6.2fx\n",
            K, dt, dt/K, 1e6*dt/(K*ntgt), (dt/K)/8.23)
    flush(stdout)
end

# stderr is fully buffered when redirected to a file on Julia 1.12, so the library's own @info
# ("num of sources", "num of hcub calculations") would otherwise appear only at exit.
flush(stderr)

get(ENV, "SKIP_BLAS", "") == "1" && exit(0)

# OpenBLAS: the other thing that differed between the two eras.
K = min(Kmax, maximum(KS))
println("\nOpenBLAS thread sensitivity at K = $K:")
for nb in (1, 96)
    BLAS.set_num_threads(nb)
    dt = run1(K)
    @printf("  BLAS threads %3d ->  %8.1f s   (%.2f s per source)\n", nb, dt, dt/K)
    flush(stdout)
end

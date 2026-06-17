# Diagnostic: is the interface-build time (and its growth with K/n_src) mainly the
# NUFFT steps (ltkm3dc = TKM type-1+type-2 FINUFFT) vs the FMM (lfmm3d) vs other?
#
# Profiles the multi/envelope build at two cutoffs and classifies each profile sample by
# whether its stack passes through the NUFFT/TKM path, the FMM path, or neither.
# (Classifying by the Julia ltkm3dc/lfmm3d wrapper frames is robust even when the C
# FFTW/FINUFFT leaves don't resolve.) Reports build time + NUFFT% for each, so the
# GROWTH (cutoff 1.5 -> 2.5) can be attributed.
#
# Run at the PRODUCTION config (96 threads) to include the FFTW thread pathology.
# Build only, no solve. Run: JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 julia -t 96 --project=<BI> diag_build_profile.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using Printf, Profile

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const TEMPL = [BI.read_xsf(joinpath(REF_DIR, "graphene_00001.xsf"))[2],
               BI.read_xsf(joinpath(REF_DIR, "graphene_00002.xsf"))[2]]
const CENTROID = [ntuple(d -> Float64(BI.density_centroid(dg)[d]), 3) for dg in TEMPL]
const A1 = (2.465, 0.0, 0.0)
const A2 = (-1.2325, 2.1347526, 0.0)
addv(a, b) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
scalev(s, a) = (s*a[1], s*a[2], s*a[3])
distv(a, b) = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2)

const LZ, EPS_IN, EPS_OUT, L = 3.35, 3.5, 1.0, 90.0
const N_QUAD, EDGE, RHS_TOL, MAXDEPTH, SUPPORT_RTOL = 6, 4, 1e-3, 12, 1e-4
const L_EC = LZ / 2.0^EDGE * 1.01
const BOXES = BI.BoxGeom[(center = (0.0, 0.0, CENTROID[1][3]), Lx = L, Ly = L, Lz = LZ)]
const EPSES = Float64[EPS_IN]

function batch(cutoff; nrange = 3)
    center = CENTROID[1]
    sites = Tuple{Int,NTuple{3,Float64}}[(1, center)]
    for t in 1:2, n1 in -nrange:nrange, n2 in -nrange:nrange
        (t == 1 && n1 == 0 && n2 == 0) && continue
        pos = addv(CENTROID[t], addv(scalev(Float64(n1), A1), scalev(Float64(n2), A2)))
        distv(pos, center) <= cutoff && push!(sites, (t, pos))
    end
    insts = Dict(id => BI.OrbitalInstance(id, t, BI.snap_orbital(TEMPL[t], CENTROID[t], pos))
                 for (id, (t, pos)) in enumerate(sites))
    return BI.assemble_lattice_batch(TEMPL, insts, [(1, id) for id in 1:length(sites)];
                                     support_rtol = SUPPORT_RTOL)
end

function build(cutoff)
    b = batch(cutoff)
    env = BI.envelope_volume_source(b)
    kmax = BI._estimate_tkm3dc_kmax(env)
    return BI.multi_dielectric_box3d_rhs_adaptive(N_QUAD, L_EC, BOXES, EPSES, env, RHS_TOL;
        eps_out = EPS_OUT, max_depth = MAXDEPTH, tkm_kmax = kmax)
end

const NUFFT_RE = r"ltkm3dc|nufft|finufft|tkm3d|fftw"
const FMM_RE = r"lfmm3d|fmm3d"

function profile_build(cutoff)
    Profile.init(n = 10^8, delay = 0.01)
    Profile.clear()
    local itf
    t = @elapsed (@profile itf = build(cutoff))
    data = Profile.fetch(include_meta = false)
    lidict = Profile.getdict(data)
    nufft = 0; fmm = 0; other = 0; total = 0
    cur_n = false; cur_f = false; started = false
    for ip in data
        if ip == 0
            if started
                total += 1
                cur_n ? (nufft += 1) : cur_f ? (fmm += 1) : (other += 1)
            end
            cur_n = false; cur_f = false; started = false
        else
            started = true
            frames = get(lidict, ip, nothing)
            if frames !== nothing
                for sf in frames
                    s = lowercase(string(sf.func) * "|" * string(sf.file))
                    occursin(NUFFT_RE, s) && (cur_n = true)
                    occursin(FMM_RE, s) && (cur_f = true)
                end
            end
        end
    end
    return (; n_pts = BI.num_points(itf), t, total,
            nufft_pct = 100nufft/max(total,1), fmm_pct = 100fmm/max(total,1),
            other_pct = 100other/max(total,1))
end

@printf("threads = %d   OMP = %s\n", Threads.nthreads(), get(ENV, "OMP_NUM_THREADS", "?"))
build(0.0)  # warm-up / compile

results = NamedTuple[]
for cutoff in (1.5, 2.5)
    r = profile_build(cutoff)
    nufft_s = r.t * r.nufft_pct / 100
    @printf("cutoff %.1f  build %.1f s  pts %d  | NUFFT %.1f%% (%.1f s)  FMM %.1f%%  other %.1f%%  [%d samples]\n",
            cutoff, r.t, r.n_pts, r.nufft_pct, nufft_s, r.fmm_pct, r.other_pct, r.total)
    push!(results, (; cutoff, r.t, r.nufft_pct, nufft_s))
    GC.gc()
end

if length(results) == 2
    dt = results[2].t - results[1].t
    dn = results[2].nufft_s - results[1].nufft_s
    @printf("\nGROWTH 1.5 -> 2.5:  total +%.1f s   NUFFT +%.1f s  (%.0f%% of the growth is NUFFT)\n",
            dt, dn, 100dn/max(dt, eps()))
end

# ground truth: what IS the 'other' 99.5%? print the top self-time frames of the
# last profiled build (cutoff 2.5, still in the Profile buffer).
println("\n-- top frames of the cutoff 2.5 build (flat, count = samples through frame) --")
Profile.print(format = :flat, sortedby = :count, mincount = 10000)
println("\nDIAG BUILD PROFILE DONE")

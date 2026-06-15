# 3D NUFFT ~50–80× slower at 96 threads than at ≤64 for specific mode-grid dimensions

## Summary

A 3D type-1 (and type-2) transform of ~3.5×10⁵ nonuniform points onto a
**(413, 375, 325)** mode grid runs in ~0.3–0.5 s at 1–64 threads but jumps to
**~27 s at 96 threads** — a ~50–80× slowdown from *adding* threads. A
neighbouring, slightly *larger* grid **(441, 401, 353)** is unaffected
(~0.4 s at 96 threads). The cost is in the FFT execution, not planning
(`finufft_makeplan` is ~0 with the default ESTIMATE planner). The machine is a
dual-socket 2×48-core node, and the cliff appears exactly at 96 = full
subscription of both sockets, which suggests an FFTW threaded-FFT / NUMA
decomposition pathology for these particular transform dimensions.

## Environment

| | |
|---|---|
| FINUFFT.jl | 3.5.2 |
| finufft_jll | 2.5.1+0 |
| FFTW_jll | 3.3.12+0 |
| Julia | 1.12.6 |
| CPU | 2× AMD EPYC 9474F (48-core each; 96 cores total) |
| OS / kernel | Rocky Linux 8.10, kernel 6.6.142 |
| upsampfac | default/auto (= 1.25 here; reproduces at forced `upsampfac=1.25`) |
| tol | 1e-3 (also reproduces at 1e-4) |

## Minimal reproducer (depends only on FINUFFT.jl)

```julia
import FINUFFT as FN
using Random, Printf

N = 356_000; tol = 1e-3
Random.seed!(1)
x = 2π .* rand(N) .- π; y = 2π .* rand(N) .- π; z = 2π .* rand(N) .- π
c = complex.(randn(N))

function bench(ms, nth)                       # min of 2 warm type-1 transforms
    FN.nufft3d1(x, y, z, c, -1, tol, ms...; nthreads = nth)
    minimum(@elapsed(FN.nufft3d1(x, y, z, c, -1, tol, ms...; nthreads = nth)) for _ in 1:2)
end

for ms in ((413, 375, 325), (441, 401, 353))
    for nth in (1, 8, 16, 32, 64, 96)
        @printf("modes %-16s  nthreads %3d : %7.2f s\n", string(ms), nth, bench(ms, nth))
    end
end
```

Run (Julia's own thread count is irrelevant; FINUFFT `nthreads` is set per call):

```
julia --project=. repro.jl     # project with FINUFFT 3.5.2 pinned
```

## Results

type-1 `nufft3d1` wall time (s):

| mode grid | 1 | 8 | 16 | 32 | 64 | **96** |
|---|--:|--:|--:|--:|--:|--:|
| (413, 375, 325) | 2.30 | 0.71 | 0.53 | 0.51 | 0.60 | **27.71** |
| (441, 401, 353) | 3.09 | 3.79 | 2.41 | 3.25 | 3.48 | **3.78** |

A guru-interface type-2 (`finufft_makeplan` + `finufft_exec`) on the same
grids shows the same effect, isolating it to the FFT execution:

| mode grid | metric | 16 thr | 96 thr |
|---|---|--:|--:|
| (413, 375, 325) | plan+setpts | 0.01 | 0.05 |
| (413, 375, 325) | **exec** | **0.40** | **27.33** |
| (441, 401, 353) | exec | 0.49 | 0.39 |

## What this is and isn't

- **Not problem size:** the slow grid (5.0×10⁷ modes) is *smaller* than the
  fast grid (6.2×10⁷). Both are 2,3,5-smooth after FINUFFT's upsampling.
- **Not the spreader / not planning:** `finufft_makeplan` (ESTIMATE) is ~0;
  the time is entirely in `finufft_exec` (the threaded FFT).
- **Not upsampfac:** auto picks 1.25; forcing `upsampfac=1.25` reproduces it,
  and forcing 2.0 is uniformly slow on both grids (separate concern).
- **Not tolerance:** reproduces at both 1e-3 and 1e-4.
- **Thread-count + dimension specific:** the slow grid is fine at 1–64 threads
  and only blows up at 96 (= all 96 cores across 2 sockets). The fast grid is
  fine at every thread count.

## Hypothesis / question

This looks like an FFTW multi-threaded plan that decomposes badly for certain
3D upsampled sizes when using the full core count on a dual-socket NUMA node
(possible cross-socket work split or a degenerate thread tiling of one axis).
Is this a known FFTW threading issue for specific sizes, and/or could FINUFFT
guard against it (e.g. nf-size selection that avoids thread-unfriendly
factorizations, or an `nthreads` heuristic for pathological dimensions)? A
per-call `nthreads ≤ 64` is an effective workaround, but 96 threads being
~50× slower than 64 on the same transform is surprising.

(Discovered in a boundary-integral solver where a volume-potential evaluation
landed on the (413,375,325) grid and ran ~50 s instead of ~1 s at 96 threads;
a neighbouring box that lands on (441,401,353) is unaffected.)

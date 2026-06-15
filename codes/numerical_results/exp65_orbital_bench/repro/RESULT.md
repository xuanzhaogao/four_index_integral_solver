# FINUFFT 96-thread FFT pathology — standalone reproduction

**Stack:** FINUFFT.jl 3.5.2, finufft_jll 2.5.1+0 (FFTW backend), Julia 1.12.6,
AMD EPYC Genoa 96-core (Flatiron worker7001), `tol = 1e-3`, N = 356,000
uniform-random nonuniform points. Reproducer: `repro_finufft_threads.jl`
(depends only on FINUFFT; pinned `Project.toml`/`Manifest.toml` in this dir).

## Observation

A 3D transform on mode grid **(413, 375, 325)** is ~**80x slower at 96 threads
than at 64**, while a nearby grid (441, 401, 353) is unaffected. The cost is in
the FFT execution (`finufft_exec`); plan time is ~0 (FFTW ESTIMATE).

```
grid (modes)                 nthr    type1(s)  t2_plan(s)  t2_exec(s)
BAD  ltkm3dc-box (413,375,325)  1       2.30      0.00       2.02
BAD  ltkm3dc-box (413,375,325)  8       0.71      0.00       0.51
BAD  ltkm3dc-box (413,375,325) 16       0.53      0.01       0.40
BAD  ltkm3dc-box (413,375,325) 32       0.51      0.01       0.35
BAD  ltkm3dc-box (413,375,325) 64       0.60      0.04       0.34
BAD  ltkm3dc-box (413,375,325) 96      27.71      0.05      27.33   <-- 80x cliff
GOOD field-box   (441,401,353)  1       3.09      0.00       2.38
GOOD field-box   (441,401,353)  8       3.79      0.01       0.63
GOOD field-box   (441,401,353) 16       2.41      0.01       0.49
GOOD field-box   (441,401,353) 32       3.25      0.01       0.40
GOOD field-box   (441,401,353) 64       3.48      0.04       1.98
GOOD field-box   (441,401,353) 96       3.78      0.05       0.39   <-- fine
```

## Characterization

- The cliff is **specific to the full core count (96)** for the BAD grid — it
  is fast (0.34-0.6 s) at 1..64 threads, then jumps to ~27 s at 96. The GOOD
  grid is fast at every thread count including 96.
- It is the threaded **FFT execution**, not planning (FFTW ESTIMATE, plan ~0).
- It is **dimension-specific.** Both grids are 2,3,5-smooth after FINUFFT's
  upsampling+rounding (nf = next-2,3,5-smooth of ~1.25x the mode count); the
  BAD grid's upsampled size evidently decomposes badly across 96 FFTW threads.

## Why it surfaced

In the graphene BIE pipeline (exp65) the screened-density volume-potential
call `TKM3D.ltkm3dc` built its per-call Fourier box from the source+target
extent, landing on mode grid (413,375,325); at the production 96 threads each
of its two transforms (type-1 + type-2) took ~27 s, i.e. the ~54 s full call.
The `PrecomputedVolumeField` path uses a margin-padded box -> (441,401,353),
which incidentally lands on a thread-friendly size and avoids the cliff.

## Fix options (NOT a global thread cap)

1. **Report upstream** (FINUFFT / FFTW): a smooth 3D size ~5e7 should not be
   80x slower at 96 vs 64 threads — looks like an FFTW threaded-plan
   decomposition pathology for specific dimensions.
2. **Pad mode dimensions to thread-friendly sizes** at the call site (what the
   field box does for free). This keeps full threading for the well-behaved
   sizes — preferable to capping `nthreads` globally.
3. Per-transform `nthreads` heuristic only for flagged-bad dimensions.

## Reproduce

```
JULIA_NUM_THREADS=1 julia --project=exp65_orbital_bench/repro \
    exp65_orbital_bench/repro/repro_finufft_threads.jl
```

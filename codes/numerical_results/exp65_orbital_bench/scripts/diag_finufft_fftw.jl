# Test the proposed fix: does FFTW_MEASURE/PATIENT avoid the FINUFFT 96-thread FFT
# pathology that FFTW_ESTIMATE (FINUFFT's default) falls into? Separates FFTW PLANNING
# (finufft_makeplan, where MEASURE pays) from EXECUTION (finufft_exec, where the
# pathology shows). Type-1 3D NUFFT on the documented-pathological mode grid.
#
# Run: OMP_NUM_THREADS=96 julia -t 2 --project=<BI> diag_finufft_fftw.jl

using Printf
import TKM3D
const FF = TKM3D.FINUFFT

const ESTIMATE, MEASURE, PATIENT = 64, 0, 32   # FFTW planner flags
const NMODES = [413, 375, 325]                 # REPORT's pathological grid (27 s @96, 0.48 s @16)
const N = 300_000
const TOL = 1e-3

x = (rand(N) .* 2π) .- π
y = (rand(N) .* 2π) .- π
z = (rand(N) .* 2π) .- π
c = randn(ComplexF64, N)

function planexec(nth, fftwflag)
    tplan = @elapsed begin
        plan = FF.finufft_makeplan(1, NMODES, -1, 1, TOL; nthreads = nth, fftw = fftwflag)
        FF.finufft_setpts!(plan, x, y, z)
    end
    local out
    texec = @elapsed out = FF.finufft_exec(plan, c)
    FF.finufft_destroy!(plan)
    return tplan, texec
end

planexec(16, ESTIMATE)  # warm-up / compile
@printf("type-1 NUFFT  modes=%s  Nsrc=%d  tol=%.0e\n\n", Tuple(NMODES), N, TOL)
@printf("%-10s %-9s %-12s %-12s\n", "fftw", "nthreads", "plan (s)", "exec (s)")
for (nm, fl) in (("ESTIMATE", ESTIMATE), ("MEASURE", MEASURE), ("PATIENT", PATIENT))
    for nth in (16, 96)
        tp, te = planexec(nth, fl)
        @printf("%-10s %-9d %-12.2f %-12.2f\n", nm, nth, tp, te)
    end
end
println("\nDIAG FINUFFT DONE")

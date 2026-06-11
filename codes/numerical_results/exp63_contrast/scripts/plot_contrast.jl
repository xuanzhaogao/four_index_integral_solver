# 6.3 figure: (a) V-error vs eps2 against per-case references,
# (b) N_iter vs eps2 with normalized GMRES residual-history inset,
# (c) screening curve V/V_vac vs eps2 (conductor-limit saturation).
# Output: figs/fig63_contrast.pdf

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using CairoMakie, LaTeXStrings

const DATA = joinpath(@__DIR__, "..", "data")
const RAW = joinpath(DATA, "raw")
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

e2s = sort([parse(Float64, match(r"^ratio_test_e2_(.+)\.jls$", f).captures[1])
            for f in readdir(RAW) if occursin(r"^ratio_test_e2_.*\.jls$", f)])
tests = Dict(e2 => load_ref(joinpath(RAW, "ratio_test_e2_$(e2).jls")) for e2 in e2s)
refs = Dict(e2 => load_ref(joinpath(RAW, "ratio_ref_e2_$(e2).jls")) for e2 in e2s)
Vvac = load_ref(joinpath(RAW, "v_vacuum.jls")).Vvac

relerr = [abs(tests[e].V - refs[e].V) / abs(refs[e].V) for e in e2s]

const CTEST = "#0072B2"
const CREF = "#D55E00"
const HIST_E2 = [6.0, 400.0, 1.0e11]                # low / high / conductor limit
const CHIST = Dict(6.0 => "#0072B2", 400.0 => "#D55E00", 1.0e11 => "#009E73")

e2lab(e) = e >= 1e10 ? L"10^{11}" : L"%$(round(Int, e))"

fig = Figure(size = (1180, 330))

ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10,
    xlabel = L"\varepsilon_2", ylabel = L"|V - V_\mathrm{ref}| / |V_\mathrm{ref}|",
    title = "(a) accuracy vs contrast, ε = 10⁻⁴", xticks = (e2s, e2lab.(e2s)))
scatterlines!(ax1, e2s, relerr; color = CTEST, marker = :circle)
hlines!(ax1, [1e-4]; color = :black, linestyle = :dot)
text!(ax1, e2s[2], 1.15e-4; text = L"\varepsilon = 10^{-4}", fontsize = 12)

ax2 = Axis(fig[1, 2]; xscale = log10,
    xlabel = L"\varepsilon_2", ylabel = L"N_\mathrm{iter}",
    title = "(b) GMRES iterations", xticks = (e2s, e2lab.(e2s)))
scatterlines!(ax2, e2s, [tests[e].niter for e in e2s]; color = CTEST,
    marker = :circle, label = "test (p=6, r=4)")
scatterlines!(ax2, e2s, [refs[e].niter for e in e2s]; color = CREF,
    marker = :utriangle, label = "reference (p=8, r=6)")
axislegend(ax2; position = :lt, framevisible = false, labelsize = 11)

# inset: normalized residual histories (test runs), selected contrasts
ax2i = Axis(fig[1, 2]; width = Relative(0.46), height = Relative(0.42),
    halign = 0.92, valign = 0.12, yscale = log10,
    xlabel = "iteration", ylabel = L"\|r_k\|/\|r_0\|",
    xlabelsize = 9, ylabelsize = 9, xticklabelsize = 8, yticklabelsize = 8)
for e in HIST_E2
    h = tests[e].history
    lines!(ax2i, 0:(length(h) - 1), h ./ h[1]; color = CHIST[e], linewidth = 1.4)
end
translate!(ax2i.blockscene, 0, 0, 150)

ax3 = Axis(fig[1, 3]; xscale = log10,
    xlabel = L"\varepsilon_2", ylabel = L"V / V_\mathrm{vac}",
    title = "(c) screening", xticks = (e2s, e2lab.(e2s)))
scatterlines!(ax3, e2s, [tests[e].V for e in e2s] ./ Vvac; color = CTEST, marker = :circle)
hlines!(ax3, [tests[1.0e11].V / Vvac]; color = "#009E73", linestyle = :dash)
text!(ax3, e2s[2], tests[1.0e11].V / Vvac * 1.02;
    text = "conductor limit", fontsize = 11, color = "#009E73")

save(joinpath(FIGS, "fig63_contrast.pdf"), fig)
save(joinpath(FIGS, "fig63_contrast.png"), fig; px_per_unit = 2)
println("wrote figs/fig63_contrast.{pdf,png}")

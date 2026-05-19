# Fig 4: convergence of near-field GL quadrature — notes

Working notes from the d-sweep / p-sweep investigation of the standard
`p × p` Gauss–Legendre near-field error for two parallel unit panels.

## Setup

- Panels: `P = [-1, 1]² × {0}` (target), `Q = [-1, 1]² × {d}` (source).
- Kernel: switched from Laplace double-layer (`z·(4πr³)⁻¹`) to **single-layer
  `1/(4π r)`** — the algebraic-branch structure is the same (`r² = 0` is a
  branch curve), but the single-layer doesn't vanish on the panel surface,
  which makes the LHS Bernstein contour and the RHS rates easier to read.
- Density: `σ_Q(η₁, η₂) = exp(−η₁² − η₂²)`, sampled at the `p × p` GL grid
  and extended off-node by tensor barycentric Lagrange — exactly what the
  production solver consumes.
- Target points: the `p²` GL nodes on `P`.
- Reference integral: HCubature on `K · P_σ` with `rtol = 1e-12`.
- Error metric (fig4 right panel): relative L²,
  `E_std = ‖I_p − I_ref‖₂ / ‖I_ref‖₂`. The plot script also stores absolute
  error `E_abs = ‖I_p − I_ref‖₂` and `‖I_ref‖₂` per `(p, d)`.

## The puzzle

Eyeballing the right panel of `fig4_near_correction.png` (E_std vs `d` for
`p ∈ {4, 6, 8}`), the empirical decay looked like

    E_std  ~  ρ(d)⁻ᵖ ,                ρ(d) = d + √(1 + d²),

i.e. **half** the canonical tensor-product Bernstein rate `ρ⁻²ᵖ`. The
question: is the rate truly half, or is something else going on?

## Diagnostic 1 — p-sweep at fixed d

`fig4_pscan.jl` fixes `d ∈ {0.3, 1.0, 3.0}` and sweeps
`p ∈ {4, 6, 8, 10, 12, 14}`, then linearly fits `log₁₀ E_std` against `p`.

| d   | ρ      | empirical slope | `−2 log₁₀ ρ` (canonical) | `−log₁₀ ρ` (half-rate) |
|-----|--------|-----------------|--------------------------|------------------------|
| 0.3 | 1.344  | **−0.328**       | **−0.257**                | −0.128                  |
| 1.0 | 2.414  | **−0.785**       | **−0.766**                | −0.383                  |
| 3.0 | 6.162  | **−1.306**       | **−1.580**                | −0.790                  |

- `d=1.0` is the cleanest: `−0.785` matches `−2 log₁₀ ρ = −0.766` essentially exactly.
- `d=0.3` is a hair steeper than the canonical line (mild pre-asymptotic).
- `d=3.0` looks shallower only because `E_std` hits the HCubature
  reference / round-off floor (`~10⁻¹⁵`) past `p ≈ 10`; the LS fit on the
  full range is dragged by the saturated tail.

**Conclusion:** the actual convergence rate **is** the canonical Bernstein
`ρ(d)⁻²ᵖ`. The empirical `ρ⁻ᵖ` in the d-sweep is an artefact of reading an
exponent off a curve where both `ρ(d)` and the prefactor change with `d`.

Output: `figs/fig4_pscan.{pdf,png}` — the data hugs the dashed `ρ⁻²ᵖ` family
and is nowhere near the dotted `ρ⁻ᵖ` family.

## Diagnostic 2 — prefactor `C(p, d) = E_abs(p, d) · ρ(d)²ᵖ`

`fig4_prefactor.jl` divides out the canonical exponential. If the canonical
form held exactly with a `d`-independent prefactor, the curves would be
flat in `d`. They are not — two regimes appear:

### Near-singular regime, `d ∈ [0.1, ~1.0]`

`C_abs(p, d)` is essentially flat and nearly p-independent.

| d    | C(p=4) | C(p=6) | C(p=8) |
|------|--------|--------|--------|
| 0.1  | 0.69   | 0.53   | 0.45   |
| 0.5  | 0.18   | 0.18   | 0.17   |
| 1.0  | **0.094** | **0.078** | **0.077** |

At `d = 1` the absolute prefactor collapses onto a single number
`C_abs ≈ 0.08` for all three `p`. **This is the regime where the canonical
Bernstein bound `E_abs ≈ 0.08 · ρ⁻²ᵖ` is tight.**

### Far-field regime, `d ≳ 1.5`

`C_abs` grows rapidly in `d`, more steeply for larger `p`. Far-field-only
(i.e. `d > 1.5`) log-log slopes of `C_abs(d)`:

| p | α (slope of `C_abs` in d) |
|---|--------------------------|
| 4 | 1.7 |
| 6 | 3.4 |
| 8 | 6.2 |

Roughly `α(p) ≈ 0.7 p`. With `ρ ≈ 2d` at large `d`, this gives
`E_abs(p, d) ~ d^(0.7p) · d⁻²ᵖ = d^(−1.3 p)`. The shape on a `log E` vs
`log d` plot at fixed `p` is therefore noticeably shallower than the naive
`d⁻²ᵖ` slope — which is exactly the "looks like `ρ⁻ᵖ`" effect that started
the investigation.

## Resolution of the puzzle

The two slopes measure complementary cuts of the same surface
`E_abs(p, d) = C(p, d) · ρ(d)⁻²ᵖ`:

- **Fixed `d`, vary `p`** ⇒ slope vs `p` is `∂ₚ log C − 2 log ρ`. In both
  regimes `∂ₚ log C ≈ 0`, so the slope is purely `−2 log ρ`. ✓
- **Fixed `p`, vary `log d`** ⇒ slope vs `log d` is
  `d (∂_{log d}) log C − 2p · d log ρ / d log d`. At large `d`,
  `d log ρ / d log d → 1`, but the prefactor term contributes `+0.7p` from
  the far-field `α(p) ≈ 0.7p`, partially cancelling the `−2p` from `ρ⁻²ᵖ`.

So both observations are consistent. There is **no half-rate**: the true
asymptotic rate in `p` is canonical Bernstein. The far-field "softening"
in `d` is because the **Bernstein bound is loose** for smooth integrands —
when `d ≫ panel size`, the kernel `1/r ≈ 1/d` varies slowly across the
panel and the actual GL error decays faster than `ρ⁻²ᵖ`. Multiplying the
actual error by `ρ⁺²ᵖ` then blows the prefactor up.

## Physical picture

```
near-singular (d ≲ L_half)           far-field (d ≳ L_half)
──────────────────────────────       ──────────────────────────────
branch curve close to [-1,1]²         branch curve far from panel
integrand quasi-singular              integrand smooth, slowly varying
Bernstein bound TIGHT                 Bernstein bound LOOSE upper bound
E_abs ≈ 0.08 · ρ⁻²ᵖ                   E_abs ≪ 0.08 · ρ⁻²ᵖ ; C grows
```

The crossover sits near `d ≈ L_half`. Below it, the dashed
`ρ⁻²ᵖ/(ρ²−1)` lines on the figure track the data. Above it, the data
peels off and decays slower than the bound predicts.

## Implications for the main figure

Current `fig4_plot.jl` uses

```julia
factors = [1.1, 2.0, 2.6]
f_temp  = x -> (x + √(1+x²))^(-2p) / ((x+√(1+x²))² − 1) / factors[i]
```

calibrated at `d ≈ 1`. The slope is correct everywhere; the per-p
constants are right where Bernstein is tight; the dashed lines diverge
below the data for `d ≳ 1.5` because Bernstein over-predicts the rate
there.

Three reasonable presentation choices:

1. **Trim x-axis** to `d ∈ [0.1, 3]` so the dashed lines track the data
   throughout the displayed range.
2. **Re-label** the dashed line as `Bernstein upper bound, ρ⁻²ᵖ/(ρ²−1)`;
   the divergence then reads correctly as "bound goes loose".
3. **Single universal prefactor**: replace the formula with
   `E_abs ≈ 0.08 · ρ⁻²ᵖ` (or, in relative-error form,
   `E_std ≈ 0.08 · ρ⁻²ᵖ / ‖I_ref(d)‖`). Cleaner, but only valid in the
   near-field.

Not yet applied — pending choice.

## File index

| file | what it does |
|------|--------------|
| `fig4_data.jl`      | Sweep `(p, d)` for `p ∈ {4,6,8}`, `d ∈ [10⁻¹, 10¹]`. Saves `E_std`, `E_abs`, `‖I_ref‖` per pair. |
| `fig4_plot.jl`      | Renders `figs/fig4_near_correction.{pdf,png}` — main figure (Bernstein ellipses LHS, `E_std` vs `d` RHS). |
| `fig4_pscan.jl`     | Diagnostic: fixed `d ∈ {0.3, 1.0, 3.0}`, sweep `p ∈ {4,…,14}`. Confirms the `ρ⁻²ᵖ` rate at fixed `d`. Saves `figs/fig4_pscan.{pdf,png}`. |
| `fig4_prefactor.jl` | Diagnostic: plots `E_abs` and `E_abs · ρ²ᵖ` vs `d`. Exposes the two-regime structure of the prefactor. Saves `figs/fig4_prefactor.{pdf,png}`. |

## To regenerate

```bash
# Data sweep (a few minutes due to HCubature reference at rtol=1e-12):
ssh worker7002 '/mnt/home/xgao1/.juliaup/bin/julia \
    --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/fig_gen \
    /mnt/home/xgao1/work/four_index_integral_solver/codes/fig_gen/fig4_data.jl'

# Diagnostics + figure:
julia --project=. fig4_pscan.jl     # ~5 min on worker7002 (heavy at p=12,14)
julia --project=. fig4_prefactor.jl # seconds (reads fig4_data.jls)
julia --project=. fig4_plot.jl      # seconds (reads fig4_data.jls)
```

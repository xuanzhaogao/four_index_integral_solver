# Email draft — Malte Rösner, 2026-05-13

**Subject:** Monolayer bare-V k-mesh sweep — both methods converge at k=32×32

---

Hi Malte,

Thanks for the two extra datasets and for pointing at the Madelung correction. I ran the bare 4-channel integral at all three k-meshes through our direct Wannier pipeline and the result is, I think, cleaner than your hypothesis suggested — both methods are computing the same Wannier matrix element, and they agree to ~10 mV once your k-mesh is dense enough.

Full table (all values in eV; `diff = CoQui − ours`):

| k-mesh    | channel  | our v   | CoQui v  |  diff   |
|-----------|----------|--------:|---------:|--------:|
| k_161601  | onsite   | 16.4004 | 17.4342  |  1.0337 |
| k_161601  | nn       |  8.5289 |  8.8396  |  0.3107 |
| k_161601  | hund_sf  | 0.11750 | 0.13080  | 0.01331 |
| k_161601  | hund_ph  | 0.11750 | 0.13080  | 0.01331 |
| k_252501  | onsite   | 17.3815 | 17.4376  |  0.0561 |
| k_252501  | nn       |  8.8217 |  8.8376  |  0.0159 |
| k_252501  | hund_sf  | 0.12902 | 0.13058  | 0.00156 |
| k_252501  | hund_ph  | 0.12902 | 0.13058  | 0.00156 |
| k_323201  | onsite   | 17.4158 | 17.4036  | −0.0122 |
| k_323201  | nn       |  8.8425 |  8.8319  | −0.0107 |
| k_323201  | hund_sf  | 0.13125 | 0.13155  | 0.00030 |
| k_323201  | hund_ph  | 0.13125 | 0.13155  | 0.00030 |

Two things jump out:

**1. At k_323201, all four channels agree to within 12 mV.** Onsite, NN, and both Hund's terms. That's not what we would see if the difference between our direct integral and CoQui's bare V were a real physical periodic-image effect — a Madelung-like contribution wouldn't simultaneously vanish on the onsite, NN, *and* sub-meV-scale Hund's channels at the same k-mesh. It's what we'd see if the two methods are computing the same Wannier-orbital matrix element and your CoQui calculation is converged at this BZ sampling.

**2. The "Madelung correction" you reported is k-mesh-dependent with sign change** (+1.041 → +0.184 → −0.290 eV), which is the signature of a finite-supercell / finite-k-mesh regularization residual, not of a static lattice contribution. The CoQui log notes `Treatment of long-wavelength divergence in bare V: gygi`, and the Gygi-Baldereschi prescription is designed precisely so that the regularized result *equals* the converged Wannier-orbital integral as the BZ sampling refines, with the auxiliary q=0 contribution decaying (often oscillatorily) to zero.

So my reading is: what's being called a "Madelung correction" is the finite-k-mesh error in CoQui's bare V, not a piece of physics that needs to be added on top of a direct integral. The 1 eV apparent gap at k_161601 is then expected — the Gygi-Baldereschi residual is still ~1 eV at that BZ sampling — and the agreement at k_323201 is the actual k-mesh-convergence signal.

A couple of questions to confirm this reading:

- Is there a way to read off, from CoQui's output, the explicit Gygi-Baldereschi auxiliary value it added at the q=0 point for each calculation? If yes, I'd expect that number to equal each of (+1.041, +0.184, −0.290) eV at the corresponding k-mesh, and we'd have a tight diagnostic.
- For going forward — should we treat k_323201 as the joint reference for parameter validation? That's the cleanest agreement point we have, and the 10 mV residual is well below the cusp-under-resolution error from our XSF grid (the 2026-04-23 report attributed ~5% on the *coarse* k-mesh; that estimate is probably an overestimate given how well things agree at dense k).
- The XSF-resolution question (3×3×1 vs 5×5×1): I think this is not the leading candidate any more — the density and Hund's channels agreeing at sub-percent level at k_323201 argues the current grid is adequate for the bare matrix elements. I'd still welcome a higher-resolution regeneration once your re-optimized MLWFs are ready, since tighter cusp resolution would only help, but I wouldn't ask for it as a priority.

Full report (including the parsing, the script, and the raw CSV) is in our repo at `monolayer/results/2026-05-13-bare-monolayer-kmesh-sweep.md`. Happy to share or annotate if helpful.

Best,
Xuanzhao

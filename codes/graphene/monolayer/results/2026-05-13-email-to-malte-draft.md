# Email draft — Malte Rösner, 2026-05-13

**Subject:** Monolayer bare-V k-mesh sweep — both methods converge at k=32×32

---

Hi Malte,

Thanks for the two extra datasets and for pointing at the Madelung correction. I ran the bare 4-channel integral at all three k-meshes through our direct Wannier pipeline, the results are shown in the table below.

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

At k_323201, all four channels agree with each other within 12 mV, showing that both methods are computing the same Wannier matrix element. The previous mismatch at k_161601 is therefore a convergence artifact, not a physical effect from periodic boundary conditions or Madelung contributions.

The XSF-resolution question (3×3×1 vs 5×5×1): I think this is not the leading candidate any more — the density and Hund's channels agreeing at sub-percent level at k_323201 argues the current grid is adequate for the bare matrix elements.

I think we can use this data to do the dielectric solves and see if the same convergence pattern holds for the screened interactions.

Best,
Xuanzhao

# Bare Coulomb Interaction Parameters for Graphene π Orbitals

Computed via real-space volume integral using VASP XSF orbital files for a 5×5×1 graphene supercell.

---

## Formula

The bare Hubbard interaction between (normalized) Wannier functions $w_m$ and $w_n$ is

$$
U_{mn}(\mathbf{R}) = \iint |\tilde{w}_m(\mathbf{r})|^2 \, v_C(\mathbf{r} - \mathbf{r}') \, |\tilde{w}_n(\mathbf{r}')|^2 \, d^3r \, d^3r'
$$

where $v_C(\mathbf{r}) = e^2 / (4\pi\varepsilon_0 |\mathbf{r}|)$ is the bare Coulomb kernel and
$\tilde{w}$ denotes the **normalized** orbital density obtained from the XSF wavefunction $\psi$ via

$$
|\tilde{w}(\mathbf{r})|^2 = \frac{|\psi(\mathbf{r})|^2}{\mathcal{N}}, \qquad \mathcal{N} = \int_{\mathrm{SC}} |\psi(\mathbf{r})|^2 \, d^3r.
$$

---

## Numerical implementation

The original reference run used `lfbc3d` (`FBCPoisson.jl`), which applies the free-boundary Poisson kernel

$$
G(\mathbf{r}) = \frac{1}{4\pi|\mathbf{r}|}
$$

and the current validation run uses `TKM3D.ltkm3dc`, which evaluates the same free-space Laplace kernel in Fourier space. In both cases the raw code output (in Å$^{-1}$) is related to the physical interaction by

$$
U_{mn}^{\rm phys} \;[\text{eV}] = \underbrace{4\pi}_{\text{restore Coulomb}} \times \underbrace{\frac{e^2}{4\pi\varepsilon_0}}_{\text{= 14.3996 eV·Å}} \times \frac{U_{mn}^{\rm raw} \;[\text{Å}^{-1}]}{\mathcal{N}_m \cdot \mathcal{N}_n}
$$

### Normalization

The XSF wavefunction $\psi$ satisfies (verified numerically):

$$
\mathcal{N} = \int_{\rm SC} |\psi(\mathbf{r})|^2 \, d^3r \approx V_{\rm prim} \simeq 84.26 \;\text{Å}^3
$$

where $V_{\rm prim} = |\mathbf{a}_1 \cdot (\mathbf{a}_2 \times \mathbf{c})| \approx 84.3\;\text{Å}^3$ is the primitive-cell volume.
The factor $\mathcal{N} \neq 1$ reflects VASP's internal normalization convention.

---

## Lattice geometry (5×5×1 supercell)

| Vector | Value (Å) |
|---|---|
| $\mathbf{A} = 5\mathbf{a}_1$ | $[12.243,\; 0,\; 0]$ |
| $\mathbf{B} = 5\mathbf{a}_2$ | $[-6.121,\; 10.615,\; 0]$ |
| $\mathbf{C}$ | $[0,\; 0,\; 15.920]$ |

Primitive lattice constant $a = 2.465$ Å, C–C bond $d_{\rm CC} = 1.42$ Å.
Orbital A is centered on sublattice A, orbital B on sublattice B (both from `graphene_0000{1,2}_5x5x1_shifted.xsf`).

---

## Results

### FBCPoisson reference

| Parameter | Pair | Distance (Å) | $U^{\rm raw}$ (Å$^{-1}$) | $U^{\rm phys}$ (eV) | Point-charge limit (eV) |
|---|---|---|---|---|---|
| $U_{00}$ | A–A on-site | 0 | 673.070 | **17.153** | — |
| $U_{01}$ | A–B nearest-neighbor | 1.420 | 346.012 | **8.817** | 10.141 |
| $U_{02}$ | A–A next-nearest | 2.465 | 218.411 | **5.566** | 5.842 |
| $U_{03}$ | A–B 2nd shell | 2.845 | 190.572 | **4.856** | 5.062 |

The point-charge limit is $e^2/(4\pi\varepsilon_0 d)$.
All inter-site values fall below the point-charge limit, as expected for extended Wannier functions.

### TKM3D validation

`compute_bare_hubbard_graphene_tkm3d.jl` reproduces the same four channels using `ltkm3dc` with a geometry-derived shared cutoff `k_{\max} = 39.0625` and tolerance sweep `10^{-2}, 10^{-3}, 10^{-4}`. Taking the `10^{-4}` run as the TKM reference gives:

| Parameter | $U^{\rm raw}$ (Å$^{-1}$) | $U^{\rm phys}$ (eV) | rel. err. at $10^{-3}$ | rel. err. at $10^{-2}$ |
|---|---|---|---|---|
| $U_{00}$ | 673.525 | **17.1645** | $3.89\times 10^{-5}$ | $1.33\times 10^{-3}$ |
| $U_{01}$ | 346.008 | **8.8170** | $3.92\times 10^{-5}$ | $1.25\times 10^{-3}$ |
| $U_{02}$ | 218.395 | **5.5657** | $3.34\times 10^{-5}$ | $9.75\times 10^{-4}$ |
| $U_{03}$ | 190.573 | **4.8562** | $3.12\times 10^{-5}$ | $9.13\times 10^{-4}$ |

These values agree closely with the `lfbc3d` reference. The practical lesson from this dataset is that letting `estimate_kcut3dc` infer an anisotropic Nyquist box from Cartesian coordinate gaps is too expensive for the skewed graphene grid; the validation script therefore uses the shared geometry-derived cutoff instead.

---

## Script

`compute_bare_hubbard_graphene.jl` — `lfbc3d` reference sweep over `N_FFT`.

`compute_bare_hubbard_graphene_tkm3d.jl` — `ltkm3dc` validation sweep over `tol`.

---

## Computational details

- FBC reference: `lfbc3d` with $N_{\rm FFT} = 256$, tolerance $10^{-6}$
- TKM validation: `ltkm3dc` with shared geometry-derived $k_{\max} = 39.0625$ and tolerances $10^{-2}, 10^{-3}, 10^{-4}$
- Input: squared XSF wavefunction (`datagrid.values .*= datagrid.values`)
- Z-shift applied to both orbitals: $-7.920155$ Å (centers graphene at $z \approx 0$)
- Off-site shifts for $U_{02}$/$U_{03}$: `VolumeSource` position shift by $\mathbf{a}_1 = [2.465, 0, 0]$ Å
- $\mathcal{N}_1 = 84.2638\;\text{Å}^3$, $\mathcal{N}_2 = 84.2723\;\text{Å}^3$

---

## Physical context

These are **bare** (unscreened) values.
After constrained-RPA screening the on-site $U_{00}$ is typically reduced to $\sim 9\text{–}10$ eV for graphene π orbitals (cf. Wehling et al., *PRL* **106**, 236805 (2011)).

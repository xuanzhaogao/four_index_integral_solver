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

The `lfbc3d` solver (`FBCPoisson.jl`) uses the free-boundary Poisson kernel

$$
G(\mathbf{r}) = \frac{1}{4\pi|\mathbf{r}|}
$$

so the raw code output (in Å$^{-1}$) is related to the physical interaction by

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

| Parameter | Pair | Distance (Å) | $U^{\rm raw}$ (Å$^{-1}$) | $U^{\rm phys}$ (eV) | Point-charge limit (eV) |
|---|---|---|---|---|---|
| $U_{00}$ | A–A on-site | 0 | 673.070 | **17.153** | — |
| $U_{01}$ | A–B nearest-neighbor | 1.420 | 346.012 | **8.817** | 10.141 |
| $U_{02}$ | A–A next-nearest | 2.465 | 218.411 | **5.566** | 5.842 |
| $U_{03}$ | A–B 2nd shell | 2.845 | 190.572 | **4.856** | 5.062 |

The point-charge limit is $e^2/(4\pi\varepsilon_0 d)$.
All inter-site values fall below the point-charge limit, as expected for extended Wannier functions.

---

## Script

`compute_bare_hubbard_graphene.jl` (same directory) — runs all four integrals and prints results.

---

## Computational details

- Solver: `lfbc3d` with $N_{\rm FFT} = 256$, tolerance $10^{-6}$
- Input: squared XSF wavefunction (`datagrid.values .*= datagrid.values`)
- Z-shift applied to both orbitals: $-7.920155$ Å (centers graphene at $z \approx 0$)
- Off-site shifts for $U_{02}$/$U_{03}$: `VolumeSource` position shift by $\mathbf{a}_1 = [2.465, 0, 0]$ Å
- $\mathcal{N}_1 = 84.2638\;\text{Å}^3$, $\mathcal{N}_2 = 84.2723\;\text{Å}^3$

---

## Physical context

These are **bare** (unscreened) values.
After constrained-RPA screening the on-site $U_{00}$ is typically reduced to $\sim 9\text{–}10$ eV for graphene π orbitals (cf. Wehling et al., *PRL* **106**, 236805 (2011)).

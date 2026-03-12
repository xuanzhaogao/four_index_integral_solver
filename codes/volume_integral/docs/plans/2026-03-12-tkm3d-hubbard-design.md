# TKM3D Hubbard Validation Design

## Goal
Add a dedicated `ltkm3dc` validation workflow in `codes/volume_integral` that reproduces the graphene bare-interaction channels `U_00`, `U_01`, `U_02`, and `U_03`, then studies their convergence as the TKM tolerance is tightened. The existing `FBCPoisson` scripts remain unchanged during this phase.

## Context
`compute_bare_hubbard_graphene.jl` currently loads two graphene XSF orbitals, constructs four `BoundaryIntegral.VolumeSource` objects, evaluates potentials with `lfbc3d`, and writes a CSV under `data/`. The migration target, `TKM3D.ltkm3dc`, already exists in the project dependency set and can estimate its own spectral cutoff from `eps` when `kmax` is omitted.

## Proposed Architecture
Introduce a new TKM-focused entry point instead of editing the current FBC script in place. The new path will share the same orbital-loading, geometry, and normalization conventions, but will evaluate each interaction channel with `ltkm3dc(eps, sources; charges, targets, pgt = 1)`. A small helper file will hold reusable logic for tolerance validation, convergence-row construction, and CSV shaping so the expensive physics script is testable through pure functions.

## Data Flow
1. Load the two XSF wavefunctions and square them into orbital densities.
2. Build `VolumeSource` objects for the on-site, nearest-neighbor, and shifted neighbor cases.
3. Form preweighted source charges as `density .* weights`.
4. For each selected tolerance, evaluate the TKM potential on the target positions and integrate against the target weighted density to obtain `U_raw`.
5. Convert `U_raw` to eV using the same normalization factors already used by the FBC script.
6. Treat the tightest tolerance run as the reference and compute per-pair relative error columns.
7. Write a convergence CSV containing `pair`, `tol`, `U_raw`, `U_ev`, relative errors, and the estimated `kcut` metadata used for the run.

## Error Handling
- Reject empty or invalid tolerance lists.
- Fail fast if the XSF files are missing.
- Check `values.ier` from `ltkm3dc` and abort on any nonzero code.
- Preserve explicit comments around units and normalization to simplify scientific review.

## Verification Strategy
- Add a focused Julia test around the pure helper logic that prepares tolerance sweeps and reference-relative error rows.
- Run the new TKM script end-to-end on a short tolerance list to confirm it produces the expected CSV structure and a zero relative error for the tightest run.
- After the short smoke test passes, run the intended tolerance sweep to inspect convergence behavior for `U_00` through `U_03`.

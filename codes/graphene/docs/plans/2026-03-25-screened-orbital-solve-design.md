# Screened Orbital Solve Design

**Goal:** Compute screened graphene orbital interactions by solving the dielectric-box boundary problem for a centered source density, then evaluating `u_int + u_scatter` on the same four orbital pairs used in the bare `volume_integral` workflow.

**Context:** The existing graphene code in this repo only analyzes screened densities and their Fourier decay. The reference orbital workflow in [`compute_bare_hubbard_graphene_tkm3d.jl`](/mnt/home/xgao1/work/four_index_integral_solver/codes/volume_integral/compute_bare_hubbard_graphene_tkm3d.jl) computes the unscreened `U_00/U_01/U_02/U_03` channels by evaluating the direct volume potential between four centered/shifted orbital densities.

## Requirements

- Use the same orbital channels as the bare workflow:
  - `U_00`: orbital 1 against orbital 1
  - `U_01`: orbital 1 against orbital 2
  - `U_02`: orbital 1 against orbital 1 shifted by `a1 = (2.465, 0, 0)`
  - `U_03`: orbital 1 against orbital 2 shifted by `a1`
- Center the input density in the dielectric box.
- Use the dielectric-box solve sequence:
  - build interface
  - evaluate RHS
  - build LHS
  - solve for surface charge `sigma`
  - evaluate direct volume potential `u_int`
  - evaluate scattered potential `u_scatter`
  - integrate `u_int + u_scatter` over target orbitals
- Sweep `SharpScreening()` and `SoftMixInversePermittivity(b)` over a small bandwidth set.
- Use:
  - `L = 90.0`
  - `LZ = 2.8`
  - `EPS_IN = 2.4`
  - `EPS_OUT = 1.0`

## Approach

Create a small reusable module for the orbital workflow rather than folding solver logic into the plotting module. The new module will:

- load and square the XSF orbitals,
- center the datagrids using the existing centroid-based shift logic,
- construct the four `VolumeSource` targets,
- apply `screened_volume_source(..., mode)` to the source density,
- build the adaptive interface with `single_dielectric_box3d_rhs_adaptive`,
- form the RHS with `Rhs_dielectric_box3d_fmm3d`,
- solve the corrected LHS with GMRES,
- evaluate `u_int` with `BoundaryIntegral.TKM3D.ltkm3dc`,
- evaluate `u_scatter` with `laplace3d_pottrg_near`,
- write a CSV table with raw and eV interaction values.

## Numerical choices

- Reuse the package’s own `kmax` estimate through `BoundaryIntegral._estimate_tkm3dc_kmax(vs)` for the direct volume potential.
- Use the corrected FMM3D LHS already used in the volume-source convergence scripts.
- Keep the source screening explicit by constructing a screened `VolumeSource` for each mode. This keeps the direct field, RHS, and total field consistent with the requested screened-density interpretation.

## Outputs

- New source module for screened orbital solves.
- New script entrypoint that runs the bandwidth sweep and writes a CSV under `data/`.
- Focused tests for helper behavior and direct-potential evaluation on synthetic sources.

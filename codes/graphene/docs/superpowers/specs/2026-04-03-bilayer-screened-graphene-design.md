# Bilayer Screened Graphene Design

**Goal:** Evaluate screened orbital interactions `U_00` through `U_05` for the bilayer-as-slab graphene model by reusing the existing screened-orbital workflow with updated slab geometry and an extended in-plane neighbor set.

**Context:** In this repository, graphene is modeled as a dielectric slab rather than as two explicit atomic sheets. For this task, "bi-layer graphene" means a single dielectric slab of thickness `d = 6.7 Å`, with the same dielectric constant used in the single-layer run and both orbital densities centered at `z = d/4`.

## Requirements

- Keep the same screened-solve pipeline used in [results/2026-03-26-screened-graphene-report.md](/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/results/2026-03-26-screened-graphene-report.md):
  - load and square the XSF orbital densities
  - center the orbitals
  - build the dielectric-box interface
  - solve for surface charge `sigma`
  - evaluate `u_int` and `u_scatter`
  - integrate against the target orbitals
- Change the slab thickness from `2.8 Å` to `6.7 Å`.
- Keep `eps_in` equal to the single-layer case.
- Shift both source and target orbitals so their centroids are located at `z = d/4 = 1.675 Å`.
- Extend the target set from four channels to six:
  - `U_00`: on-site orbital 1 to orbital 1
  - `U_01`: on-site orbital 1 to orbital 2
  - `U_02`: center orbital to the existing `a1`-shifted like-orbital target
  - `U_03`: center orbital to the existing `a1`-shifted cross-orbital target
  - `U_04`: center orbital to the fourth-nearest-neighbor target
  - `U_05`: center orbital to the fifth-nearest-neighbor target
- Follow the same reporting style as the March 26 screened-graphene workflow and produce actual numerical values for `U_00` to `U_05`.

## Approach

Extend the shared orbital helper layer instead of writing a one-off script. The source construction path will continue to apply one common in-plane centering shift derived from orbital 1, but it will also accept a target `z` center so the same helper can place both orbitals at `z = d/4`. The pair-construction layer will be generalized from a hard-coded four-channel tuple to a small explicit table of pair labels and lattice translations. This keeps the screened solver unchanged while allowing the script to request the six bilayer-slab channels in a reproducible way.

The production script will be updated to use `Lz = 6.7`, preserve `eps_in = 2.4`, and emit the six-channel interaction table. The output should include the same raw and eV quantities already reported by the current script, plus the effective `z` placement and the target translation metadata so the bilayer slab run is traceable.

## Neighbor Convention

This repository already fixes the first shifted channel by `a1 = (2.465, 0, 0)`. For this task, `U_04` and `U_05` will be implemented by adding two more explicit in-plane translation vectors to the pair table using the graphene lattice convention already used by the workflow. The implementation should define these vectors in one place so the test suite can verify the intended positions directly.

## Testing

- Extend the pair-target tests to verify the six labels and the exact translation vectors.
- Add a test that verifies the new source-centering helper places both orbitals at the requested `z` center without changing the shared in-plane centering behavior.
- Keep the existing one-shot full-target interaction test unchanged except for any necessary pair-count expectations.

## Output

- Updated solver helpers that support six target channels and configurable `z` placement.
- Updated production script for the `d = 6.7 Å` slab case.
- Verified numerical results for `U_00` to `U_05`.

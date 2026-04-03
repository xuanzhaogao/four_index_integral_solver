# Paper Shell Labeling Design

**Goal:** Align the reported `U_00` through `U_05` labels in the slab-model graphene workflow with the shell ordering used in Rosner et al., while keeping the current dielectric-slab approximation unchanged.

## Context

The current code evaluates six interaction channels for a graphene orbital source in a dielectric slab. The solver geometry is not the same as the bilayer graphene model in [PhysRevB.92.085102.pdf](/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/docs/PhysRevB.92.085102.pdf): the paper uses AB-stacked bilayer graphene with explicit layer separation and a momentum-dependent effective dielectric function, whereas the code uses one slab with both orbitals centered at the same `z`.

The user does not want to replace the slab model. They only want the `U_00` to `U_05` shell labels to mirror the paper’s shell convention as closely as possible under the slab approximation.

## Requirements

- Keep the current screened slab workflow and geometry unchanged.
- Make the shell-label convention explicit in code instead of leaving it implicit in ad hoc translation vectors.
- Define `U_00` through `U_05` as the paper-style shell ordering:
  - `U_00`: on-site
  - `U_01`: nearest neighbor
  - `U_02`: next-nearest neighbor
  - `U_03`: third shell
  - `U_04`: fourth shell representative
  - `U_05`: fifth shell representative
- Document that the slab approximation does **not** reproduce the paper’s full AB-bilayer sublattice splitting for every second shell.
- Preserve the current production run behavior and output values unless the shell ordering itself is wrong.

## Effective Dielectric Interpretation

The paper does not define one scalar effective dielectric constant for freestanding bilayer graphene. Instead it introduces a momentum-dependent effective two-dimensional dielectric function `ε_eff^2D(q)` in Eq. (11), with vacuum surroundings corresponding to `ε2 = ε3 = 1`. The current slab code still uses the repository’s scalar dielectric-box parameters `eps_in` and `eps_out`; these should remain unchanged, but the report should state clearly that this is a slab-model approximation rather than the paper’s `ε_eff^2D(q)` model.

## Approach

Add one explicit layout constant for the paper-style shell representatives and use it as the default pair layout. The constant should encode both the shell label and a short human-readable description so the chosen mapping is obvious from the code and easy to test. The existing displacement vectors already appear to follow the shell ordering by distance, so the implementation should focus on making that convention explicit, testable, and documented rather than changing solver logic.

Add a regression test that checks the representative shell distances and verifies the shell labels increase monotonically in distance from the source centroid. Update the report to state that `U_00` through `U_05` follow the paper’s shell numbering under a slab approximation and that the effective dielectric treatment differs from the paper’s momentum-dependent WFCE dielectric function.

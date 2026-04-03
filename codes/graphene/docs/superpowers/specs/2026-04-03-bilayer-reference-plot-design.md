# Bilayer Reference Plot Design

**Goal:** Generate a one-off comparison figure that plots the current `eps_in = 2.4` bilayer-slab softmix results against the freestanding bilayer graphene reference values from Table I of Rosner et al.

## Context

The repository already has a plotting workflow for comparing screened Hubbard values against cRPA references, but that module is currently tailored to the earlier four-channel monolayer-style comparison. The user asked not to refactor that shared path and instead wants the same visual style reused for a one-off bilayer comparison figure.

## Requirements

- Leave the shared comparison module and its existing plot script untouched.
- Create a dedicated script that reads the current `data/screened_hubbard_graphene.csv`.
- Plot the `eps_in = 2.4` softmix `U_00` through `U_05` values against the freestanding BLG Table I reference values.
- Reuse the same general look as the existing comparison figure:
  - one line per pair over bandwidth
  - scatter markers at sampled bandwidths
  - dashed horizontal reference lines
- Make the title and output filename clearly bilayer-specific.
- Add a focused test for the six-value bilayer reference dictionary used by the new plot.

## Approach

Create a standalone plotting script under `scripts/` that defines its own pair order and BLG reference dictionary, then renders a Makie figure using the same structure as the existing `plot_screened_hubbard_vs_crpa.jl` script. Add the minimal reference-value regression to the comparison tests so the hard-coded Table I values are checked in one place.

The figure should be framed as a slab-approximation comparison to the paper reference, not as an exact WFCE reproduction.

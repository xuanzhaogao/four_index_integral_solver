# Screened Hubbard cRPA Plot Design

Date: 2026-03-26

## Goal

Create a combined comparison plot of the computed screened Hubbard interactions versus softmix bandwidth `b`, overlaid with the graphene cRPA reference values from `data/PhysRevLett.106.236805.pdf`.

## Data

- Computed results: `data/screened_hubbard_graphene.csv`
- cRPA reference values from Table I of `PhysRevLett.106.236805.pdf`
  - `U_00 = 9.3 eV`
  - `U_01 = 5.5 eV`
  - `U_02 = 4.1 eV`
  - `U_03 = 3.6 eV`

## Plot Layout

- One combined figure.
- `x` axis: bandwidth `b`
- `y` axis: interaction energy in `eV`
- Four solid curves with markers for the computed `U_00`, `U_01`, `U_02`, `U_03`
- Four matching horizontal dashed lines for the cRPA references
- Output: `figs/screened_hubbard_vs_crpa.png`

## Implementation

- Add a small helper module to:
  - load the CSV
  - filter/sort the softmix rows
  - expose the cRPA reference constants in one place
- Add a plotting script using `CairoMakie`
- Add lightweight tests for the helper logic so the plotting data path is verified separately from rendering

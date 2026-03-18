# Repository Guidelines

## Project Structure & Module Organization
This repository is a small Julia workflow for volume-integral experiments on graphene orbitals. Top-level scripts such as `compute_bare_hubbard_graphene.jl`, `plot_bare_interaction.jl`, `orbital_orbital.jl`, and `orbital_orbital_fft.jl` are the main entry points. Generated tabular outputs live in `data/`, figures in `figs/`, and narrative analysis notes in `results/`. Package metadata is tracked in `Project.toml` and `Manifest.toml`. Several scripts expect XSF density inputs outside this repo at `../../density_data/`.

## Build, Test, and Development Commands
Use the project environment for every run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. compute_bare_hubbard_graphene.jl
julia --project=. plot_bare_interaction.jl
```

The first command installs pinned dependencies. The second regenerates `data/hubbard_graphene.csv`. The third rebuilds `figs/hubbard_graphene_convergence.svg`. Use `julia --project=. orbital_orbital.jl` or `julia --project=. orbital_orbital_fft.jl` for exploratory solver checks.

## Coding Style & Naming Conventions
Follow existing Julia script style: 4-space indentation, descriptive local names, and `snake_case` for files and variables. Keep research constants explicit (`z_shift`, `e2_4pieps0`) and prefer short comments that capture physics assumptions or units instead of restating code. There is no repo formatter config, so match the surrounding style and keep imports grouped at the top.

## Testing Guidelines
There is no automated `test/` suite yet. Treat script reruns as regression checks: rerun the relevant Julia entry point, confirm the CSV or SVG output updates as expected, and review any changed values in `results/` notes before committing. If you add reusable logic, create `test/runtests.jl` with Julia's `Test` standard library and make the new test runnable with `julia --project=. test/runtests.jl`.

## Commit & Pull Request Guidelines
Recent history uses short imperative commits, sometimes with prefixes such as `feat:`, `docs:`, and `chore:`. Prefer that convention over vague messages like `update`. For pull requests, include: the scientific goal, scripts touched, regenerated artifact paths, dependency changes, and any required external inputs. Attach updated plots when figure output changes, and call out assumptions about files under `../../density_data/`.

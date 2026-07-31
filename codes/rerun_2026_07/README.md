# Rerun campaign, July 2026 — edge correction + diagonal preconditioning

Reruns every article result affected by the two LHS defects described in
`BoundaryIntegral.jl/design/plans/2026-07-31-lhs-conditioning-edges-and-preconditioning.md`.

Nothing here overwrites the data behind the submitted manuscript. Every job
writes to a `_v2`-suffixed location; the published trees stay byte-identical for
comparison.

## What was actually wrong, per section

| Article section | edges as published | precond as published | why it must be rerun |
|---|---|---|---|
| §6.1 — Fig 8, Table 1 | on | **off** | iteration counts only; the *V* errors should reproduce |
| §6.3 — Tables 2, 3 | on | **off** | iteration counts and every timing |
| §6.4 — Fig 10, Table 4 | **off** | **off** | **values**, plus timings |

§6.4 is the only section whose numbers can move. It went through
`solve_dielectric_lattice_batch`, which had no `correct_edges` argument at all,
so the flag could not be set from the campaign TOML.

## The trap that shaped these scripts

`Harness.jl` (§6.1) and `run_multicube.jl` (§6.3) call `Krylov.gmres` and
`Krylov.block_gmres` **directly**. They never go through
`solve_dielectric_box3d_block`, so BI's new `precondition = true` default does
**not** reach them. Rerunning those two experiments unmodified would have
reproduced the published iteration counts exactly and burned eight exclusive
nodes proving nothing.

Both now take a `PRECONDITION` environment variable (default `1`) and apply
`N = Diagonal(BI.dielectric_diagonal_scaling(interface))` as a right
preconditioner, matching what BI does internally. Set `PRECONDITION=0` to
recover the published behaviour.

The §6.4 campaign path needs no such change: `solve_batch_core` passes neither
flag, so it picks up both new defaults automatically.

## Job list

Run `scripts/gen_rerun_campaigns.jl` **before** jobs 05 and 06.

| job | produces | nodes | est. cost | notes |
|---|---|---|---|---|
| `01_fig8_sweep.sbatch` | Fig 8 (b), (c) | 1 | ~2 h | includes the p=8 reference and the no-edge ablation |
| `02_tab1_contrast.sbatch` | Table 1 | 1 | ~2 h | |
| `03_tab2_threads.sbatch` | Table 2 | 1 × array 0–7 | ~1 h each | **needs 04's K=1 record to plot** |
| `04_tab3_ksweep.sbatch` | Table 3 | 1 × array 0–6 | ~4 h total | K=31 peaked at 538 GB |
| `05_lattice_het3x.sbatch` | Table 4, Fig 9 geometry | 10 | was 38.5 node-h | resumable |
| `06_lattice_conv_l3.sbatch` | Fig 10 (a), (b), §6.4 *U* values | 10 | was 83 node-h | resumable, heaviest |
| `07_lec_single.sbatch` | Fig 10 (c) | 1 | ~1 h | levels 1–5 sequentially |

Expect the §6.4 jobs to be **cheaper** than the originals despite the extra
near-pair assembly: the published batches ran at 61–86 GMRES iterations with
neither fix, and the solve stage was 84% of the cost.

## Submission order

The `genoa` partition has 32 nodes, so the two ten-node campaigns plus the two
arrays cannot all run at once. Suggested order:

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07

# 0. retarget the lattice campaigns (writes campaigns/*_v2.toml, refuses to
#    clobber an existing ceph root)
julia --project=../lattice_scale scripts/gen_rerun_campaigns.jl

# 1. the cheap ones first -- they validate the patched drivers before you commit
#    ten nodes to anything
sbatch jobscripts/01_fig8_sweep.sbatch
sbatch jobscripts/02_tab1_contrast.sbatch

# 2. Sec. 6.3. Run 04 before 03 is plotted: plot_scaling.jl reads niter_ref from
#    04's converged K=1 record, and the published 36 is stale once the
#    preconditioner is on.
sbatch jobscripts/04_tab3_ksweep.sbatch
sbatch jobscripts/03_tab2_threads.sbatch

# 3. Sec. 6.4, the expensive ones. 06 is the one that decides whether the
#    published tensor values stand.
sbatch --nodes=10 jobscripts/06_lattice_conv_l3.sbatch
sbatch --nodes=10 jobscripts/05_lattice_het3x.sbatch
sbatch jobscripts/07_lec_single.sbatch
```

Jobs 05 and 06 are resumable — phase status is derived from files on disk, so
resubmitting the identical line continues where a timeout left off.

## After the runs

```bash
BASE=lattice_conv_l3     TAG=_v2 julia --project=../lattice_scale scripts/compare_rerun.jl
BASE=lattice_10x10_het3x TAG=_v2 julia --project=../lattice_scale scripts/compare_rerun.jl
```

This is the number the paper needs: how far the tensor moved, and how much the
iteration count fell. It reports on-site *U* both ways with the shift, the
whole-tensor relative difference, the exchange-symmetry error both ways, and
per-batch iteration counts.

Sanity checks worth applying to its output:

- The preconditioner must not move σ beyond `gmres_rtol`. A materially larger
  shift is a bug in the scaling, not a result.
- On a single-contrast geometry the preconditioner is a scalar multiple of `I`,
  so the iteration count must be *identical*. That is BI's own cheapest guard
  against a sign or indexing error.

## Provenance

Every job sources `jobscripts/_provenance.sh` and logs the BoundaryIntegral and
TKM3D commit hashes plus a list of modified files. This matters right now
because **the two fixes are uncommitted working-tree edits** in
`~/codes/BoundaryIntegral.jl` — five modified files, no commit to cite. Commit
them before submitting, or the logs will be the only record of what ran.

## Files changed outside this directory

The rerun needed six drivers patched. All changes are additive and default to
the new behaviour; the published behaviour is still reachable by environment
variable.

| file | change |
|---|---|
| `numerical_results/common/Harness.jl` | `precondition` kwarg (env `PRECONDITION`, default on), applied as right preconditioner to `Krylov.gmres` |
| `numerical_results/exp66_multicube/scripts/run_multicube.jl` | same for `Krylov.block_gmres`; records `precondition` in the output record |
| `exp61_convergence/scripts/run_fig1_sweep.jl` | `RERUN_TAG` output separation |
| `exp61_convergence/scripts/run_fig1_noedges.jl` | `RERUN_TAG` |
| `exp61_convergence/scripts/analyze_fig1.jl` | `RERUN_TAG` |
| `exp63_contrast/scripts/run_contrast_ratio.jl` | `RERUN_TAG` |
| `exp63_contrast/scripts/analyze_contrast.jl` | `RERUN_TAG` |
| `lattice_scale/scripts/lec_conv_single.jl` | `RERUN_TAG` on the campaign name and output TSV |

The frozen copies in `../../../BIE_ERI_results/` are untouched and still reflect
exactly what produced the submitted figures.

# Full-Target Pair Evaluation Design

Date: 2026-03-26

## Goal

Remove target batching from the screened orbital interaction evaluation so each pair channel evaluates `u_int` and `u_scatter` on the full target orbital in one call.

## Motivation

The previous implementation split each target orbital into `20_000`-point batches. That caused the corrected surface target operator to be rebuilt repeatedly for the same pair channel, which amplified the cost of the `laplace3d_pottrg_fmm3d_corrected_hcubature` path.

## Design

1. Keep the source/interface/GMRES workflow unchanged.
2. Change `evaluate_screened_pair_interaction` to:
   - pass `target_vs.positions` directly to the volume and scattered potential evaluators
   - integrate the returned full-target vectors against `target_vs.weights .* target_vs.density`
3. Remove the `batch_size` plumbing from the orbital solve and script entrypoint.
4. Add a regression test that injects fake volume/scatter evaluators and verifies the full target set is evaluated in exactly one call per evaluator.

## Expected Outcome

- Fewer repeated corrected `pottrg` builds.
- Higher peak memory during per-pair target evaluation.
- Simpler control flow in the pair interaction path.

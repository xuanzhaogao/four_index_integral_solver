# Sourced by every rerun job. Records exactly what code produced the run --
# essential here because the two LHS fixes are (as of writing) UNCOMMITTED
# working-tree edits in the BoundaryIntegral checkout.
provenance () {
  echo "=================== PROVENANCE ==================="
  echo "host        : $(hostname)   $(date -u '+%F %T UTC')"
  echo "slurm job   : ${SLURM_JOB_ID:-none}  array task ${SLURM_ARRAY_TASK_ID:-none}"
  echo "julia       : $($J --version)"
  for R in /mnt/home/xgao1/codes/BoundaryIntegral.jl /mnt/home/xgao1/codes/TKM3D.jl; do
    echo "--- $R"
    echo "    commit  : $(git -C $R rev-parse HEAD)"
    echo "    dirty   : $(git -C $R status --porcelain | wc -l) modified file(s)"
    git -C $R status --porcelain | sed 's/^/      /'
  done
  echo "RERUN_TAG   : ${RERUN_TAG:-<unset>}"
  echo "=================================================="
}

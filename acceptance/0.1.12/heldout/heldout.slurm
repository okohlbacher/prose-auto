#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --time=06:00:00
#SBATCH --job-name=pa0112
#SBATCH --output=pa0112-%x-%A_%a.out
# PREREG 0.1.12 §10.1, arrays A (verify, windows), B (task), B2 (best_psm) and C (scenario) on BinAC2:
#   sbatch --array=1-N%40 heldout.slurm STAGE        STAGE: verify | windows | task | scenario
#   sbatch --array=1-N%16 --mem=1000gb --time=2-00:00:00 heldout.slurm best_psm   (high-memory nodes: ~900 GB peak)
#   a final unspecific best_psm search (> 1 TB) is not run here: B2 leaves best_psm/himem_due, and once the outputs
#   are back on Ceph, `harness.py himem` runs them on a 2.2 TB IBMI node, one at a time, under a lease
# A dry run sets HELDOUT_DRY=<file of KEY=VALUE lines> (PROSE_AUTO_H, _RUNS, _WRONGDB, _FREEZE, HELDOUT_NO_JEV)
# A task that exits 1 has a unit due its rerun (§10.1): resubmit those indices once with twice the memory and time,
#   sbatch --mem=256gb --time=12:00:00 --array=I,J heldout.slurm STAGE   (B2: --mem=2000gb --time=4-00:00:00)
# a unit failing again is final (missing, <unit>.failed). Outputs return to $PROJ/h/ at their paths under H.
# The frozen IBMI tree is mirrored in the workspace (<ws>/ceph -> /ceph inside the container); --cpus-per-task counts
# hyperthreads, so the threads are half of it. Each task checks the container image against FREEZE.md, the harness
# checks every executable and FASTA, and a mismatch fails the task.
set -euo pipefail
STAGE=$1
W=$(ws_find proseauto)
PROJ=/pfs/10/project/bw26g016/prose-auto-0112
T=/ceph/ibmi/abi/oliver/prose-auto-20261029
ACC=$T/prose-auto/acceptance/0.1.12
SIF=$W/ubuntu2404.sif
TH=$(( SLURM_CPUS_PER_TASK / 2 ))
DRY=()
[ -n "${HELDOUT_DRY:-}" ] && mapfile -t DRY < "$HELDOUT_DRY"
FREEZE_MD=$W$ACC/FREEZE.md
for kv in "${DRY[@]}"; do [ "${kv%%=*}" = PROSE_AUTO_FREEZE ] && FREEZE_MD=$W${kv#*=}; done
want=$(grep -o 'BinAC2 container image (SHA-256, measured on BinAC2): [0-9a-f]*' "$FREEZE_MD" | awk '{print $NF}')
have=$(sha256sum "$SIF" | cut -d' ' -f1)
[ "$want" = "$have" ] || { echo "container image $have is not the frozen $want"; exit 3; }
PUB=$TMPDIR/publish.list
rc=0
apptainer exec --bind "$W/ceph:/ceph" --bind "$TMPDIR:$TMPDIR" "$SIF" env -i HOME="$HOME" PATH=/usr/bin:/bin \
  TMPDIR_HELDOUT="$TMPDIR" HELDOUT_PUBLISH="$PUB" HELDOUT_THREADS=$TH OMP_NUM_THREADS=$TH OPENBLAS_NUM_THREADS=$TH \
  MKL_NUM_THREADS=$TH "${DRY[@]}" \
  "$T/env/bin/python" "$ACC/heldout/harness.py" "$STAGE" --index "$SLURM_ARRAY_TASK_ID" || rc=$?
# publish exactly this task's outputs (harness.PUBLISH, relative to H), also after a failure: the attempt counts
OUTH=$T/v0112/heldout
for kv in "${DRY[@]}"; do [ "${kv%%=*}" = PROSE_AUTO_H ] && OUTH=${kv#*=}; done
DEST=$PROJ/h$([ -n "${HELDOUT_DRY:-}" ] && echo _dry)
mkdir -p "$DEST"
[ -s "$PUB" ] && rsync -a -r --files-from="$PUB" "$W$OUTH/" "$DEST/"
echo "TASK_END $STAGE $SLURM_ARRAY_TASK_ID rc=$rc"
exit $rc

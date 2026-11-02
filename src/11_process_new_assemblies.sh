#!/usr/bin/env bash
# Move any oatk output that's landed in data/oatk_assemblies/ into
# data/{mito,plastid}/, then submit the core analysis suite
# (qc_basic_stats -> annotation -> repeats/synteny/phylogeny) as one LSF
# job. Meant for cron (daily - see README.md's "Keeping the dataset
# current" section for the crontab line).
#
# Safe to run even when nothing's finished assembling yet - every step
# here is then a cheap no-op. Also safe to run while some oatk jobs from
# 10_update_from_darwin.sh are still in flight:
# 03_move_assemblies.sh only moves files that already exist for a given
# species, so a species whose job hasn't finished yet (missing a file
# oatk hasn't written out yet) just gets picked up on a later run once
# more of its files appear - `mv` on the same filesystem is a rename, not
# a copy, so it can't corrupt a file oatk still has open for writing. Same
# self-healing pattern already documented in analysis/README.md for
# run_all.sh racing a revision batch, just applied one step earlier.
#
# analysis/run_all.sh is itself incremental (qc_basic_stats 01-03,
# annotation, repeats, synteny all skip a species whose output is already
# newer than its input - see analysis/common/io_utils.py), so running it
# unscoped here every day is correct AND cheap: it only actually
# (re)computes species that are new or whose assembly changed since the
# last run, not the whole ~1250-species dataset from scratch.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

TS="$(date +%Y%m%d_%H%M%S)"
mkdir -p logs
LOG="logs/process_new_assemblies_${TS}.log"

{
  echo "[info] $(date -Iseconds) starting"

  src/03_move_assemblies.sh

  echo "[info] submitting core analysis suite refresh"
  bsub -n 8 -q normal \
    -o "logs/run_all_${TS}.out" -e "logs/run_all_${TS}.err" \
    analysis/run_all.sh --organelle both

  echo "[info] $(date -Iseconds) done (analysis suite runs independently on LSF - check logs/run_all_${TS}.out/.err)"
} >>"${LOG}" 2>&1

echo "log: ${LOG}"

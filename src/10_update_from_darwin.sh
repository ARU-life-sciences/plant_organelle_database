#!/usr/bin/env bash
# Pull a fresh species/read manifest from DTOL's live archive
# (01_make_meta.sh always does a full rescan - by design, so it reflects
# the archive's current state rather than accumulating stale entries) and
# submit oatk for anything not yet assembled.
#
# Meant for cron (weekly - see README.md's "Keeping the dataset current"
# section for the crontab line and why it's split from
# 11_process_new_assemblies.sh). This only diffs the manifest and submits
# bsub jobs, then exits - it does NOT wait for them to finish; that's
# 11_process_new_assemblies.sh's job, run separately and more often.
#
# Usage: src/10_update_from_darwin.sh [coverage] [max_fasta_gb]
#   (both forwarded to 02_submit_oatk_missing_mito.sh - defaults 80, 16)
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

TS="$(date +%Y%m%d_%H%M%S)"
mkdir -p logs
LOG="logs/update_from_darwin_${TS}.log"

{
  echo "[info] $(date -Iseconds) starting DTOL refresh"

  BEFORE="logs/meta_snapshot_${TS}.before.tsv"
  AFTER="logs/meta_snapshot_${TS}.after.tsv"
  [[ -f meta/file.txt ]] && cp meta/file.txt "${BEFORE}"

  src/01_make_meta.sh
  cp meta/file.txt "${AFTER}"

  if [[ -f "${BEFORE}" ]]; then
    new_species="$(comm -13 \
      <(cut -f1 "${BEFORE}" | sort -u) \
      <(cut -f1 "${AFTER}" | sort -u))"
  else
    echo "[info] no prior manifest - treating every species as new"
    new_species="$(cut -f1 "${AFTER}" | sort -u)"
  fi

  n_new="$(grep -c . <<<"${new_species}" || true)"
  echo "[info] ${n_new} species new to the DTOL archive since last refresh"
  if [[ -n "${new_species}" ]]; then
    printf '%s\n' "${new_species}" > "logs/new_species_${TS}.txt"
    echo "[info] new species list -> logs/new_species_${TS}.txt"
  fi

  # 02_submit_oatk_missing_mito.sh is already idempotent (skips any
  # species with both mito and plastid already in data/), so it's safe to
  # run against the whole refreshed manifest every time, not just the
  # diffed new species - it naturally also retries anything that's never
  # successfully assembled at all, new or old (a prior submission that
  # never produced a usable GFA, for instance).
  src/02_submit_oatk_missing_mito.sh "${1:-80}" "${2:-16}"

  echo "[info] $(date -Iseconds) done"
} >>"${LOG}" 2>&1

echo "log: ${LOG}"

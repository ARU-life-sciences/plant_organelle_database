#!/usr/bin/env bash
set -euo pipefail

# Resolve repo root (works even if called from ./src/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

MITO_DIR="${ROOT_DIR}/data/mito"
META_DIR="${ROOT_DIR}/meta"
OUT_TSV="${META_DIR}/dead_ends.tsv"

mkdir -p "${META_DIR}"

if [[ ! -d "${MITO_DIR}" ]]; then
  echo "[err] Missing directory: ${MITO_DIR}" >&2
  exit 1
fi

# Optional sanity checks (don't hard-fail; record errors in TSV)
if ! command -v gfatk >/dev/null 2>&1; then
  echo "[warn] gfatk not found in PATH; will mark all as gfatk_missing" >&2
fi

echo "[info] Root:      ${ROOT_DIR}"
echo "[info] Mito dir:  ${MITO_DIR}"
echo "[info] Output:    ${OUT_TSV}"
echo

shopt -s nullglob

# Write TSV header immediately so you always get a file.
printf "species\tgfa_path\tdead_end_nodes\tstatus\n" > "${OUT_TSV}"

for sp_dir in "${MITO_DIR}"/*; do
  [[ -d "${sp_dir}" ]] || continue
  sp="$(basename "${sp_dir}")"

  gfa_files=( "${sp_dir}"/*.mito.gfa )
  if (( ${#gfa_files[@]} == 0 )); then
    printf "%s\t\tNA\tmissing_gfa\n" "$sp" >> "${OUT_TSV}"
    continue
  fi

  gfa="${gfa_files[0]}"

  if ! command -v gfatk >/dev/null 2>&1; then
    printf "%s\t%s\tNA\tgfatk_missing\n" "$sp" "$gfa" >> "${OUT_TSV}"
    continue
  fi

  # Run gfatk stats; do NOT let a failure kill the whole script.
  if ! stats_out="$(gfatk stats "$gfa" 2>&1)"; then
    printf "%s\t%s\tNA\tgfatk_error\n" "$sp" "$gfa" >> "${OUT_TSV}"
    continue
  fi

  # Parse first "Dead-end nodes:" integer using awk (no rg/sed dependency).
  dead_end="$(
    printf "%s\n" "$stats_out" | awk '
      /Dead-end nodes:/ {
        for (i=1; i<=NF; i++) if ($i ~ /^[0-9]+$/) { print $i; exit }
      }
    '
  )"

  if [[ -z "${dead_end}" ]]; then
    printf "%s\t%s\tNA\tparse_error\n" "$sp" "$gfa" >> "${OUT_TSV}"
    continue
  fi

  printf "%s\t%s\t%s\tok\n" "$sp" "$gfa" "$dead_end" >> "${OUT_TSV}"
done

echo "[done] Wrote ${OUT_TSV}"

#!/usr/bin/env bash
set -euo pipefail

# Resolve repo root even if script is run from ./src/
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

OATK_DIR="${ROOT_DIR}/data/oatk_assemblies"
MITO_DIR="${ROOT_DIR}/data/mito"
PLASTID_DIR="${ROOT_DIR}/data/plastid"

if [[ ! -d "${OATK_DIR}" ]]; then
  echo "[err] Missing directory: ${OATK_DIR}" >&2
  exit 1
fi

mkdir -p "${MITO_DIR}" "${PLASTID_DIR}"

# so globs that don't match expand to nothing
shopt -s nullglob

echo "[info] Root:      ${ROOT_DIR}"
echo "[info] OATK src:  ${OATK_DIR}"
echo "[info] Mito dst:  ${MITO_DIR}"
echo "[info] Pltd dst:  ${PLASTID_DIR}"
echo

moved_any=0

for sp_dir in "${OATK_DIR}"/*; do
  [[ -d "${sp_dir}" ]] || continue

  sp="$(basename "${sp_dir}")"
  echo "[info] Species: ${sp}"

  # 1) Delete utg files
  utg_files=( "${sp_dir}"/*utg* )
  if (( ${#utg_files[@]} > 0 )); then
    echo "  - deleting ${#utg_files[@]} *utg* files"
    rm -f -- "${utg_files[@]}"
  else
    echo "  - no *utg* files"
  fi

  # 2) Move mito files
  mito_files=(
    "${sp_dir}"/*.mito.*
    "${sp_dir}"/*annot_mito*
  )
  if (( ${#mito_files[@]} > 0 )); then
    dest="${MITO_DIR}/${sp}"
    mkdir -p "${dest}"
    echo "  - moving ${#mito_files[@]} mito files -> ${dest}"
    mv -f -- "${mito_files[@]}" "${dest}/"
    moved_any=1
  else
    echo "  - no mito files to move"
  fi

  # 3) Move plastid files
  pltd_files=(
    "${sp_dir}"/*.pltd.*
    "${sp_dir}"/*annot_pltd*
  )
  if (( ${#pltd_files[@]} > 0 )); then
    dest="${PLASTID_DIR}/${sp}"
    mkdir -p "${dest}"
    echo "  - moving ${#pltd_files[@]} plastid files -> ${dest}"
    mv -f -- "${pltd_files[@]}" "${dest}/"
    moved_any=1
  else
    echo "  - no plastid files to move"
  fi

  echo
done

if (( moved_any == 0 )); then
  echo "[warn] Nothing was moved. (No matching files found?)"
else
  echo "[done] Move/delete pass completed."
fi

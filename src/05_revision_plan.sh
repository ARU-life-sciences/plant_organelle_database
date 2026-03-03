#!/usr/bin/env bash
set -euo pipefail

# Build a unified per-species revision plan by combining:
#  - meta/dead_ends.tsv (species, gfa_path, dead_end_nodes, status)
#  - logs/*.out (TERM_MEMLIMIT / job killed)
#
# Output:
#  - meta/revision_plan.tsv
#
# Action priority:
#   1) mem_fix        (memkilled==1)
#   2) deadend_fix    (dead_end_nodes>0 and status==ok)
#   3) rerun_missing  (status==missing_gfa)
#   4) inspect_parse  (status==parse_error)
#   5) none

DEAD_TSV="${1:-meta/dead_ends.tsv}"
LOGDIR="${2:-logs}"
OUT="${3:-meta/revision_plan.tsv}"

if [[ ! -f "${DEAD_TSV}" ]]; then
  echo "[ERROR] Missing dead-ends TSV: ${DEAD_TSV}" >&2
  exit 1
fi
if [[ ! -d "${LOGDIR}" ]]; then
  echo "[ERROR] Missing logs dir: ${LOGDIR}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUT}")"

tmp_mem="$(mktemp)"
trap 'rm -f "${tmp_mem}"' EXIT

# Collect mem-killed species (unique)
rg -H "TERM_MEMLIMIT|job killed" "${LOGDIR}" \
  | sed -E 's#^'"${LOGDIR}"'/##; s/\.out:.*$//' \
  | sort -u > "${tmp_mem}"

# Build revision plan by joining memkill set with dead_ends.tsv
awk -F'\t' -v OFS='\t' -v memfile="${tmp_mem}" '
  BEGIN {
    # load memkilled species into a set
    while ((getline line < memfile) > 0) {
      if (line != "") mem[line] = 1
    }
    close(memfile)
  }

  NR==1 {
    # dead_ends.tsv header expected: species gfa_path dead_end_nodes status
    print "species","gfa_path","dead_end_nodes","dead_status","memkilled","dead_end_gt0","action"
    next
  }

  {
    sp=$1
    gfa=$2
    den=$3
    st=$4

    mk = (sp in mem) ? 1 : 0
    gt0 = (st=="ok" && den+0>0) ? 1 : 0

    # action priority
    action="none"
    if (mk==1) action="mem_fix"
    else if (gt0==1) action="deadend_fix"
    else if (st=="missing_gfa") action="rerun_missing"
    else if (st=="parse_error") action="inspect_parse"
    else action="none"

    print sp,gfa,den,st,mk,gt0,action
  }
' "${DEAD_TSV}" > "${OUT}"

echo "[OK] wrote ${OUT}"
echo "[OK] counts by action:"
cut -f7 "${OUT}" | tail -n +2 | sort | uniq -c

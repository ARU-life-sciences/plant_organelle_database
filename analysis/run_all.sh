#!/usr/bin/env bash
# Run the full analysis suite in dependency order:
#   qc_basic_stats (the QC gate) -> annotation -> {repeats, synteny, phylogeny}
#
# Usage:
#   analysis/run_all.sh [--species-list FILE] [--organelle mito|pltd|both]
#                        [--force] [--qc-status pass[,flag[,fail]]] [--skip-heavy]
#                        [--phylo-threads N]
# --phylo-threads overrides 04_build_tree.py's default (2 - deliberately
# conservative for interactive/login-node use); pass the core count you
# actually requested from LSF when running this as a batch job.
#
# Meant to be run either interactively for a --species-list-restricted
# mini-project (fast: seconds to low minutes), or wrapped in a single LSF
# job for the full dataset, e.g.:
#   bsub -n 8 -q normal -o analysis_run.out -e analysis_run.err \
#     analysis/run_all.sh --organelle both
# No individual step needs more than modest single-node parallelism - this
# is deliberately NOT a ~1250-job LSF array.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common/args.sh"
source "${SCRIPT_DIR}/common/tool_paths.sh"
source "${SCRIPT_DIR}/common/logging.sh"

SKIP_HEAVY=0
PHYLO_THREADS=""
PASSTHROUGH=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-heavy) SKIP_HEAVY=1; shift ;;
    --phylo-threads) PHYLO_THREADS="$2"; shift 2 ;;
    *) PASSTHROUGH+=("$1"); shift ;;
  esac
done
parse_common_args "${PASSTHROUGH[@]}"

# NOTE: `[[ cond ]] && arr+=(...)` as a bare statement is unsafe under `set -e` -
# when cond is false (the common case: no --species-list, --force off), the
# whole compound command's exit status is 1, which would abort this entire
# script. Use explicit if/fi instead.
FLAGS=(--organelle "${ORGANELLE}")
if [[ -n "${SPECIES_LIST}" ]]; then
  FLAGS+=(--species-list "${SPECIES_LIST}")
fi
if [[ "${FORCE}" == 1 ]]; then
  FLAGS+=(--force)
fi
QC_FLAGS=(--qc-status "${QC_STATUS}")

log_info "=== qc_basic_stats ==="
python3 "${SCRIPT_DIR}/qc_basic_stats/src/01_gfa_stats.py" "${FLAGS[@]}"
python3 "${SCRIPT_DIR}/qc_basic_stats/src/02_contig_stats.py" "${FLAGS[@]}"
python3 "${SCRIPT_DIR}/qc_basic_stats/src/03_gene_matrix.py" "${FLAGS[@]}"
python3 "${SCRIPT_DIR}/qc_basic_stats/src/04_core_gene_list.py" --organelle "${ORGANELLE}"
python3 "${SCRIPT_DIR}/qc_basic_stats/src/05_qc_summary.py" --organelle "${ORGANELLE}"

log_info "=== annotation ==="
python3 "${SCRIPT_DIR}/annotation/src/01_gene_table.py" "${FLAGS[@]}" "${QC_FLAGS[@]}"
python3 "${SCRIPT_DIR}/annotation/src/02_gene_fasta_extract.py" "${FLAGS[@]}"

if [[ "${SKIP_HEAVY}" == 1 ]]; then
  log_info "--skip-heavy set, stopping after annotation"
  exit 0
fi

log_info "=== repeats ==="
python3 "${SCRIPT_DIR}/repeats/src/01_find_repeats.py" "${FLAGS[@]}" "${QC_FLAGS[@]}"
python3 "${SCRIPT_DIR}/repeats/src/02_flag_recomb_repeats.py" $([[ -n "${SPECIES_LIST}" ]] && echo --species-list "${SPECIES_LIST}")
python3 "${SCRIPT_DIR}/repeats/src/03_repeats_summary.py"

log_info "=== synteny (MTPT mode) ==="
python3 "${SCRIPT_DIR}/synteny/src/01_pairwise_paf.py" --mode mtpt "${FLAGS[@]}" "${QC_FLAGS[@]}"
python3 "${SCRIPT_DIR}/synteny/src/02_dotplot.py"
python3 "${SCRIPT_DIR}/synteny/src/03_mtpt_summary.py"

log_info "=== phylogeny ==="
python3 "${SCRIPT_DIR}/phylogeny/src/01_select_core_genes.py" --organelle "${ORGANELLE}"
python3 "${SCRIPT_DIR}/phylogeny/src/02_align_genes.py" --organelle "${ORGANELLE}" $([[ "${FORCE}" == 1 ]] && echo --force)
python3 "${SCRIPT_DIR}/phylogeny/src/03_concat_alignment.py" --organelle "${ORGANELLE}" $([[ "${FORCE}" == 1 ]] && echo --force)
BUILD_TREE_ARGS=(--organelle "${ORGANELLE}")
if [[ "${FORCE}" == 1 ]]; then
  BUILD_TREE_ARGS+=(--force)
fi
if [[ -n "${PHYLO_THREADS}" ]]; then
  BUILD_TREE_ARGS+=(--threads "${PHYLO_THREADS}")
fi
python3 "${SCRIPT_DIR}/phylogeny/src/04_build_tree.py" "${BUILD_TREE_ARGS[@]}"

log_info "=== done ==="

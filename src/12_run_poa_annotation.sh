#!/usr/bin/env bash
# Dataset-wide re-run of the retired gene-calling/editing/trans-splicing/
# GFF3-export stages, now owned by the standalone `poa` package (see
# README.md's "Analysis suite" section and analysis/README.md's "Gene
# calling / editing / splicing / GFF3 export: retired, now `poa`"). This
# repo's own qc_basic_stats (core suite, NOT retired) is deliberately left
# alone - these stages read its existing qc_summary.tsv as their QC gate
# (--qc-status default: pass) rather than rebuilding it.
#
# Runs poa's own bundled stage scripts (python3 -c 'import poa; print(poa.PIPELINE_DIR)')
# against this repo as the data root (PLANT_ORGANELLE_DATA_ROOT=.), in the
# order poa's own README gives, skipping poa's bundled qc_basic_stats 01-06
# (this repo's own copy of that module stays authoritative - see above).
#
# Three of these steps rebuild their aggregate output table from scratch
# every run (see analysis/README.md's documented hazard) - never pass
# --species-list to them once full-dataset results exist, so none of the
# invocations below do.
#
# Usage: bq -n 16 -m 24G -t 48h -J poa_annotation -- src/12_run_poa_annotation.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PLANT_ORGANELLE_DATA_ROOT="$(pwd)"
export TMPDIR=/tmp

A="$(python3 -c 'import poa; print(poa.PIPELINE_DIR)')"
GENE_CALLS="analysis/denovo_annotation/results/gene_calls.tsv"
N="${LSB_DJOB_NUMPROC:-16}"
JOBS=$(( N / 2 ))
[[ "${JOBS}" -ge 1 ]] || JOBS=1

log() { echo "[$(date '+%F %T')] $*" >&2; }

log "=== denovo_annotation: nhmmscan ==="
python3 "${A}/denovo_annotation/src/00_run_nhmmscan.py" --organelle both --jobs "${JOBS}" --cpu 2
log "=== denovo_annotation: trnascan ==="
python3 "${A}/denovo_annotation/src/01_run_trnascan.py" --organelle both --jobs "${JOBS}" --thread 2
log "=== denovo_annotation: barrnap (nhmmer rRNA) ==="
python3 "${A}/denovo_annotation/src/02_run_barrnap.py" --organelle both --jobs "${N}" --threads 1
log "=== denovo_annotation: combined GFF ==="
python3 "${A}/denovo_annotation/src/03_build_combined_gff.py" --organelle both
log "=== denovo_annotation: gene_calls.tsv (full rebuild, unscoped) ==="
python3 "${A}/denovo_annotation/src/04_build_gene_calls.py" --organelle both

log "=== editing: RNA-editing-aware ORF correction (orfedit, full rebuild, unscoped) ==="
python3 "${A}/editing/src/01_scan_editing.py" --organelle both --threads "${N}" --gene-calls "${GENE_CALLS}"

log "=== trans_splicing: spliced-gene reconstruction (transsplice, full rebuild, unscoped) ==="
python3 "${A}/trans_splicing/src/01_reconstruct.py" --threads "${N}" --gene-calls "${GENE_CALLS}"

log "=== unitig_coords: unitig map ==="
python3 "${A}/unitig_coords/src/00_build_unitig_map.py" --organelle both --jobs "${N}"
log "=== unitig_coords: linearization QC ==="
python3 "${A}/unitig_coords/src/01_linearization_qc.py" --organelle both

log "=== gff_export: merged hierarchical GFF3 per species ==="
python3 "${A}/gff_export/src/01_build_gff.py" --organelle both --gene-calls "${GENE_CALLS}"

log "=== done ==="

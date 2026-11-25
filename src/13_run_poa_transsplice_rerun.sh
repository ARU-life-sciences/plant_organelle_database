#!/usr/bin/env bash
# Dataset-wide re-run of just `trans_splicing` and `gff_export`, after
# upgrading `transsplice` (0.3.0 -> 0.3.1: an exon on a non-cis junction is
# now placed next to its neighbour / round a circular contig's origin,
# instead of arbitrarily taking either of two identical repeat copies - see
# poa's CHANGELOG.md "Unreleased" entry) and reinstalling `poa` itself (its
# own gff_export code changed too: GFF junctions=cis,trans,... and
# exception=trans-splicing).
#
# Narrower than src/12_run_poa_annotation.sh on purpose: this fix only
# touches trans-splice reconstruction and its GFF rendering. denovo_annotation
# (gene calling), editing (RNA-editing-aware ORF correction) and
# unitig_coords are all upstream of/independent from transsplice, so redoing
# their expensive scans (nhmmscan/tRNAscan-SE/barrnap, ~90 min of this
# dataset's last full run) would be pure waste.
#
# Validated first against Arabidopsis_thaliana: confirmed real, not just
# relabelling - plastid rps12's exon 2 moved from an IR copy 60kb+ from its
# neighbour to the adjacent copy, same score, matching the CHANGELOG's own
# description; nad1/nad2/nad5 gained correct exception=trans-splicing
# labels with coordinates unchanged.
#
# trans_splicing/01_reconstruct.py has no incremental skip (always a full
# recompute over its manifest - see its own --force help text) and rebuilds
# its aggregate table from scratch every run (analysis/README.md's
# documented hazard) - never pass --species-list here once full-dataset
# results exist, so this script doesn't.
#
# Usage: bq -n 16 -m 16G -t 6h -J poa_transsplice_rerun -- src/13_run_poa_transsplice_rerun.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PLANT_ORGANELLE_DATA_ROOT="$(pwd)"
export TMPDIR=/tmp

A="$(python3 -c 'import poa; print(poa.PIPELINE_DIR)')"
GENE_CALLS="analysis/denovo_annotation/results/gene_calls.tsv"
N="${LSB_DJOB_NUMPROC:-16}"

log() { echo "[$(date '+%F %T')] $*" >&2; }

log "=== trans_splicing: spliced-gene reconstruction (transsplice 0.3.1, full rebuild, unscoped) ==="
python3 "${A}/trans_splicing/src/01_reconstruct.py" --threads "${N}" --gene-calls "${GENE_CALLS}"

log "=== gff_export: merged hierarchical GFF3 per species ==="
python3 "${A}/gff_export/src/01_build_gff.py" --organelle both --gene-calls "${GENE_CALLS}"

log "=== done ==="

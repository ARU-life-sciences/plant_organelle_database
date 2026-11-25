#!/usr/bin/env bash
# Re-run just gff_export, after src/13_run_poa_transsplice_rerun.sh's
# trans_splicing stage completed (new transsplice 0.3.1 results already in
# analysis/trans_splicing/results/) but its gff_export stage was killed
# partway through (job 164202, farm-wide IO slowdown - see hub notes/commit
# messages around 2026-11-25) - it had rewritten only ~1205 of 2441
# species' GFFs when killed, leaving an inconsistent partial mix of old and
# new GFFs in the working tree. build_gff() rewrites every species
# unconditionally each run, so re-running it to completion overwrites that
# partial set with a fully consistent one - no revert needed first.
#
# Usage: bq -n 16 -m 16G -t 6h -J poa_gff_export_only -- src/14_run_poa_gff_export_only.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PLANT_ORGANELLE_DATA_ROOT="$(pwd)"
export TMPDIR=/tmp

A="$(python3 -c 'import poa; print(poa.PIPELINE_DIR)')"

echo "[$(date '+%F %T')] === gff_export: merged hierarchical GFF3 per species (retry) ===" >&2
python3 "${A}/gff_export/src/01_build_gff.py" --organelle both --gene-calls analysis/denovo_annotation/results/gene_calls.tsv
echo "[$(date '+%F %T')] === done ===" >&2

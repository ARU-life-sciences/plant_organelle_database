#!/usr/bin/env bash
# Installs `poa` (the plant organelle annotator gene calling/editing/
# trans-splicing/GFF3-export package this repo's own analysis/
# denovo_annotation, editing, trans_splicing, unitig_coords and gff_export
# modules were retired in favor of - see README.md's "Analysis suite"
# section) into an already-activated conda/mamba environment. Its own
# external-tool dependencies (gfatk, orfedit, transsplice, hmmer,
# tRNAscan-SE) are pinned in environment.yml, not built here.
#
# Usage:
#   mamba env create -f environment.yml
#   conda activate plant_organellar_database
#   ./install.sh
#
# `poa` isn't on bioconda/PyPI yet (name collision with an unrelated,
# existing `poa` package - Partial Order Alignment), so it's installed
# from this cluster's sibling checkout rather than a registry. Update
# POA_REPO below (or export it before running this script) if that
# checkout moves, or switch this to a `pip install poa` /
# `mamba install poa` line once it's republished under a new name and
# that lands - see plant_organelle_annotator's own README for status.
set -euo pipefail

POA_REPO="${POA_REPO:-/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator}"

if [[ ! -d "${POA_REPO}" ]]; then
  echo "[err] poa checkout not found at ${POA_REPO} (set POA_REPO to override)." >&2
  exit 1
fi

echo "[info] installing poa from ${POA_REPO}..."
pip install "${POA_REPO}"
poa check

echo "[info] done. Verifying tool_paths.sh resolves everything cleanly:"
bash -c 'source analysis/common/tool_paths.sh'
echo "[info] if no [warn] lines appeared above, every tool this pipeline needs is on PATH."

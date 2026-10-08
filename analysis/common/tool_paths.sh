#!/usr/bin/env bash
# Single source of truth for tool paths used across analysis/.
# Source this from any bash script: `source "$(dirname "${BASH_SOURCE[0]}")/../../common/tool_paths.sh"`
#
# Every tool resolves PATH first (so `mamba env create -f environment.yml
# && conda activate plant_organellar_database` - or any other install that
# puts these on PATH - just works), falling back to this cluster's known
# absolute install location only if PATH lookup fails. The Sanger-specific
# fallback paths are harmless to keep for external installs: `command -v`
# always wins when the tool is genuinely on PATH, so these paths only ever
# matter on this cluster's own un-conda'd interactive/bsub shells.

GFATK="${GFATK:-$(command -v gfatk || echo "${HOME}/.cargo/bin/gfatk")}"
MINIMAP2="${MINIMAP2:-$(command -v minimap2 || echo /software/team301/minimap2/minimap2)}"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools || echo /software/team301/samtools/samtools)}"
SEQKIT="${SEQKIT:-$(command -v seqkit || echo /software/team301/seqkit)}"
BLASTN="${BLASTN:-$(command -v blastn || echo /software/team301/ncbi-blast-2.16.0+/bin/blastn)}"
MAKEBLASTDB="${MAKEBLASTDB:-$(command -v makeblastdb || echo /software/team301/ncbi-blast-2.16.0+/bin/makeblastdb)}"
CDHIT_EST="${CDHIT_EST:-$(command -v cd-hit-est || echo /software/team301/cdhit/cd-hit-est)}"
# The top-level /software/team301/mafft wrapper is broken (hardcoded to look for
# helper binaries under /usr/local/libexec/mafft, which doesn't exist here) -
# only relevant to that Sanger-specific fallback path; a conda-installed mafft
# on PATH doesn't have this problem and needs no MAFFT_BINARIES override.
if command -v mafft >/dev/null 2>&1; then
  MAFFT="${MAFFT:-$(command -v mafft)}"
else
  MAFFT="${MAFFT:-/software/team301/mafft-7.525-with-extensions/core/mafft}"
  export MAFFT_BINARIES="${MAFFT_BINARIES:-/software/team301/mafft-7.525-with-extensions/core}"
fi
AMAS="${AMAS:-$(command -v AMAS.py || command -v amas || echo /software/team301/AMAS/amas/AMAS.py)}"
IQTREE2="${IQTREE2:-$(command -v iqtree2 || echo /software/team301/iqtree-2.4.0-Linux-intel/bin/iqtree2)}"
GFATOOLS="${GFATOOLS:-$(command -v gfatools || echo /software/team301/gfatools/gfatools)}"
PYTHON3="${PYTHON3:-python3}"

# Gene calling, RNA editing, trans-splicing and GFF3 export (previously
# denovo_annotation/editing/trans_splicing/unitig_coords/gff_export here)
# are handled by the `poa` package now - see README.md's "Analysis suite"
# section. `poa` resolves its own tools (nhmmscan, nhmmer, tRNAscan-SE,
# orfedit, transsplice) independently; nothing here needs to know about
# them.

# Bandage - needs QT_QPA_PLATFORM=offscreen for headless rendering (no Xvfb
# needed). On this cluster, installed as a self-contained conda env rather
# than a bioconda `bandage` package (no bioconda access from compute nodes
# when it was set up) - a plain `mamba install bandage` (as in
# environment.yml) works fine for external installs and needs none of the
# BANDAGE_ENV machinery below.
if command -v Bandage >/dev/null 2>&1; then
  BANDAGE="${BANDAGE:-$(command -v Bandage)}"
else
  BANDAGE_ENV="${BANDAGE_ENV:-/nfs/users/nfs_m/mb39/miniconda3/envs/bandage}"
  BANDAGE="${BANDAGE:-${BANDAGE_ENV}/bin/Bandage}"
fi
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"

# GraphAligner (only used by linearize/'s gfatk-resolve fallback tier). On
# this cluster, resolved to a stable shpc-wrapper absolute path (equivalent
# to `module load graphaligner/...`) so non-interactive bsub jobs don't
# depend on the modules system being sourced in their shell - a plain
# `mamba install graphaligner` (as in environment.yml) puts it on PATH
# directly for external installs.
GRAPHALIGNER="${GRAPHALIGNER:-$(command -v GraphAligner || echo /software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/graphaligner/1.0.19--hdcf5f25_1/bin/GraphAligner)}"

for tool_var in GFATK MINIMAP2 SAMTOOLS SEQKIT BLASTN MAKEBLASTDB CDHIT_EST MAFFT IQTREE2 GFATOOLS BANDAGE GRAPHALIGNER; do
  tool_path="${!tool_var}"
  if [[ ! -x "${tool_path}" ]]; then
    echo "[warn] tool_paths.sh: ${tool_var}=${tool_path} is not executable/found" >&2
  fi
done

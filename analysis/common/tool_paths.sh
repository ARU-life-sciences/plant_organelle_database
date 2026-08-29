#!/usr/bin/env bash
# Single source of truth for absolute tool paths used across analysis/.
# Source this from any bash script: `source "$(dirname "${BASH_SOURCE[0]}")/../../common/tool_paths.sh"`

GFATK="${GFATK:-$(command -v gfatk || echo "${HOME}/.cargo/bin/gfatk")}"
MINIMAP2="${MINIMAP2:-/software/team301/minimap2/minimap2}"
SAMTOOLS="${SAMTOOLS:-/software/team301/samtools/samtools}"
SEQKIT="${SEQKIT:-/software/team301/seqkit}"
BLASTN="${BLASTN:-/software/team301/ncbi-blast-2.16.0+/bin/blastn}"
MAKEBLASTDB="${MAKEBLASTDB:-/software/team301/ncbi-blast-2.16.0+/bin/makeblastdb}"
CDHIT_EST="${CDHIT_EST:-/software/team301/cdhit/cd-hit-est}"
# The top-level /software/team301/mafft wrapper is broken (hardcoded to look for
# helper binaries under /usr/local/libexec/mafft, which doesn't exist here).
# The working install is mafft-7.525-with-extensions/core/mafft, and it needs
# MAFFT_BINARIES pointed at that same core/ dir (that's where countlen,
# addsingle etc. actually live) - verified directly, not from install docs.
MAFFT="${MAFFT:-/software/team301/mafft-7.525-with-extensions/core/mafft}"
export MAFFT_BINARIES="${MAFFT_BINARIES:-/software/team301/mafft-7.525-with-extensions/core}"
AMAS="${AMAS:-/software/team301/AMAS/amas/AMAS.py}"
IQTREE2="${IQTREE2:-/software/team301/iqtree-2.4.0-Linux-intel/bin/iqtree2}"
GFATOOLS="${GFATOOLS:-/software/team301/gfatools/gfatools}"
PYTHON3="${PYTHON3:-python3}"

# Bandage (v0.9.0) + a bundled BLAST+ - installed as a self-contained conda env
# (no bioconda access from this cluster's compute nodes, so this was copied from
# a working install rather than `conda create`d fresh). Needs QT_QPA_PLATFORM=
# offscreen for headless rendering (no Xvfb needed) - verified the OTHER Bandage
# install on this cluster (/software/tola/images/bandage-0.8.1.sif, a singularity
# image) has broken headless SVG text rendering, so use this env specifically,
# not that image.
BANDAGE_ENV="${BANDAGE_ENV:-/nfs/users/nfs_m/mb39/miniconda3/envs/bandage}"
BANDAGE="${BANDAGE:-${BANDAGE_ENV}/bin/Bandage}"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-offscreen}"

for tool_var in GFATK MINIMAP2 SAMTOOLS SEQKIT BLASTN MAKEBLASTDB CDHIT_EST MAFFT IQTREE2 GFATOOLS BANDAGE; do
  tool_path="${!tool_var}"
  if [[ ! -x "${tool_path}" ]]; then
    echo "[warn] tool_paths.sh: ${tool_var}=${tool_path} is not executable/found" >&2
  fi
done

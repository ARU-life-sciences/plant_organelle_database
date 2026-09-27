#!/usr/bin/env bash
# Single source of truth for absolute tool paths used across analysis/.
# Source this from any bash script: `source "$(dirname "${BASH_SOURCE[0]}")/../../common/tool_paths.sh"`

GFATK="${GFATK:-$(command -v gfatk || echo "${HOME}/.cargo/bin/gfatk")}"
# orfedit: RNA-editing-tolerant gene boundary corrector (analysis/editing/) -
# same external-Rust-tool pattern as GFATK above.
ORFEDIT="${ORFEDIT:-$(command -v orfedit || echo "${HOME}/.cargo/bin/orfedit")}"
# transsplice: trans-spliced gene reconstruction (analysis/trans_splicing/) -
# depends on orfedit as a Rust library, same external-tool pattern otherwise.
TRANSSPLICE="${TRANSSPLICE:-$(command -v transsplice || echo "${HOME}/.cargo/bin/transsplice")}"
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

# denovo_annotation/ - oatk-independent core gene calling (nhmmscan against
# oatkDB's own gene-family databases, decoupled from oatk's assembler) plus
# tRNA/rRNA calling.
NHMMSCAN="${NHMMSCAN:-/software/team301/hmmer-3.4/src/nhmmscan}"
HMM_TO_GFF="${HMM_TO_GFF:-$(command -v hmm_to_gff || echo "${HOME}/.cargo/bin/hmm_to_gff")}"
FILTER_TBLOUT="${FILTER_TBLOUT:-$(command -v filter_tblout || echo "${HOME}/.cargo/bin/filter_tblout")}"
TRNASCAN="${TRNASCAN:-/software/team301/tRNAscan-SE/tRNAscan-SE}"
BARRNAP="${BARRNAP:-/software/team301/barrnap/bin/barrnap}"
OATKDB_MITO_FAM="${OATKDB_MITO_FAM:-/software/team301/OatkDB/viridiplantae_mito_v20250217.fam}"
# Not built yet as of 2026-09-27 - oatkdb's own NCBI/edirect fetch fails with
# persistent SSL errors from normal LSF compute-node egress (confirmed: works
# fine from farm22-head2 and from the `transfer` queue, fails from `long`/
# `normal` compute nodes) - build via the `transfer` queue, not a plain bsub,
# when this is next attempted.
OATKDB_PLTD_FAM="${OATKDB_PLTD_FAM:-/software/team301/OatkDB/viridiplantae_pltd_v20260927.fam}"

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

# GraphAligner (only used by linearize/'s gfatk-resolve fallback tier) - a shpc
# container wrapper, equivalent to `module load graphaligner/1.0.19--hdcf5f25_1`
# but resolved to a stable absolute path so non-interactive bsub jobs don't
# depend on the modules system being sourced in their shell.
GRAPHALIGNER="${GRAPHALIGNER:-/software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/graphaligner/1.0.19--hdcf5f25_1/bin/GraphAligner}"

for tool_var in GFATK ORFEDIT TRANSSPLICE MINIMAP2 SAMTOOLS SEQKIT BLASTN MAKEBLASTDB CDHIT_EST MAFFT IQTREE2 GFATOOLS BANDAGE GRAPHALIGNER NHMMSCAN HMM_TO_GFF FILTER_TBLOUT TRNASCAN BARRNAP; do
  tool_path="${!tool_var}"
  if [[ ! -x "${tool_path}" ]]; then
    echo "[warn] tool_paths.sh: ${tool_var}=${tool_path} is not executable/found" >&2
  fi
done

# A plant organellar genome database

This repository contains exploratory assemblies of plant mitochondrial and plastid genomes generated with `oatk` from Darwin Tree of Life long-read datasets (PacBio).

The goal is not just to “get a mitogenome”, but to build a comparative resource for studying:

- Organelle genome structure (graphs, repeats, multipartite forms)
- Gene content variation
- DNA exchange between mitochondria, plastids, and the nucleus
- Assembly graph topology (dead-ends, bubbles, repeats)

These are real research-grade assemblies, and many are structurally complex — especially the mitochondria.

## What to do with this data?

This is an ideal resource for short bioinformatics projects in R, Python, or Bash, and especially useful if you want to build skills in:

- minimap2
- GFA parsing
- Repeat analysis
- Comparative genomics
- Data visualisation

## Analysis suite (`analysis/`)

The exercises below are still good starting points, but a lot of this is
now automated, incremental, and QC-gated across the whole dataset in
[`analysis/`](analysis/README.md):

- **qc_basic_stats** — is a given assembly actually trustworthy? (Don't
  skip this — see below.)
- **annotation** — oatk's own gene calls, reshaped into one tidy table.
- **repeats** — repeat families and putative recombination-mediating
  repeats (the mechanism behind multipartite mitochondrial structure).
- **synteny** — the MTPT exercise just below, automated for every species,
  plus dotplots and a pairwise mode for comparing species in a genus.
- **phylogeny** — a partitioned ML tree per organelle from single-copy
  marker genes.
- **orf_scan** — non-core ORFs screened against Pfam for transposable
  element / mitovirus content (opt-in — heavier than the rest, see its
  README).

Every module works on the full ~1250-species dataset by default, but also
takes a `--species-list` for a mini-project on a handful of species (a
genus, a family, whatever question you're asking). Run the lightweight
ones together with `analysis/run_all.sh` — see
[`analysis/README.md`](analysis/README.md) for usage.

**Before trusting any of it**: with ~1250 semi-automated assemblies, some
are empty, fragmented, or otherwise suspect. `analysis/qc_basic_stats`
computes a `pass`/`flag`/`fail` per species and every other module filters
against it by default — see
[`analysis/qc_basic_stats/README.md`](analysis/qc_basic_stats/README.md)
for what's actually being checked and why.

## DNA transfer between mito and plastid

Detect DNA transfer between plastid and mitochondria (MTPTs)

A classic phenomenon in plants is plastid DNA inserted into mitochondria (often called MTPTs).

This is now automated for the whole dataset in `analysis/synteny/` (PAFs,
dotplots, and a summary table per species) — see
[`analysis/synteny/README.md`](analysis/synteny/README.md). The manual
command below is still a good first exercise for understanding what's
actually happening under the hood.

Simple workflow:

```bash
minimap2 -x asm5 mito.fasta plastid.fasta > mito_vs_plastid.paf
```

From this you can:

- Quantify total aligned length
- Identify insertion breakpoints
- Compare across species
- Test whether MTPT burden correlates with genome size

In R you could:

- Summarise alignment coverage
- Plot transfer length vs genome size
- Compare clades (monocots vs dicots)

## Notes

### A note on duplicate runs

A handful of species have both a bare (`Species.mito.gfa`) and a
parameter-infixed (`Species.k1001.s31.c80.mito.gfa`) run sitting side by
side in `data/`, because `src/03_move_assemblies.sh` doesn't guard against
a re-run leaving old and new files together. If you're working with these
files directly, don't assume the newer-looking (infixed) file is the
better one — for at least one species in this dataset it's actually a
failed, empty rerun, while the original bare file is the valid assembly.
`analysis/common/species_discovery.py` handles this correctly (picks
whichever run has the most complete, non-empty file set); if you're
writing your own script against `data/` directly, check file sizes.

### A note on default parameters
Default parameters are:

k = 1001
s = 31
c = 100

This of course needs to be tuned for some assemblies, but we are doing some exploratory stuff here.

### Revising the outputs

For memory failures, we increase coverage and decrease max GB.

```bash
src/07_revision_submit_oatk_from_list.sh \
  --list rev/lists/to_resubmit_memfix.txt \
  --coverage 120 \
  --max-gb 15 \
  --suffix memfix_c120
```

For dead-ends we lower coverage, which introduces rarer sequence.

```bash
src/07_revision_submit_oatk_from_list.sh \
  --list rev/lists/to_resubmit_deadend.txt \
  --coverage 40 \
  --max-gb 15 \
  --suffix deadend_c40  
```

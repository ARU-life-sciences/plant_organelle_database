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

## DNA transfer between mito and plastid

Detect DNA transfer between plastid and mitochondria (MTPTs)

A classic phenomenon in plants is plastid DNA inserted into mitochondria (often called MTPTs).

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

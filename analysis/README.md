# Analysis suite

A lightweight, incremental set of analyses on top of the oatk mito/plastid
assemblies in `data/mito/` and `data/plastid/`. Built for short
undergraduate bioinformatics projects — each module is cheap enough to run
on the whole dataset (~1250 species) on a login node or a single modest
LSF job, and cheap enough to re-run after new species or assembly fixes
land.

No nuclear genome data is used anywhere in this suite — that keeps
everything here fast. If you want mito/nuclear or plastid/nuclear
questions (NUMTs, transposon families, etc.), see the (data-frozen,
heavier) `mito_structural_variation` repo instead.

## Modules, in dependency order

```
qc_basic_stats  →  annotation  →  ┬─ repeats
                                   ├─ synteny
                                   └─ phylogeny
qc_basic_stats  →  orf_scan          (opt-in, not in run_all.sh - see below)
qc_basic_stats  →  topology_plots    (opt-in, not in run_all.sh - see below)
```

1. **[qc_basic_stats](qc_basic_stats/README.md)** — is this assembly real
   and complete enough to trust? Per-species graph/contig stats, gene
   completeness, and a QC gate (`qc_summary.tsv`) that every module below
   filters against.
2. **[annotation](annotation/README.md)** — reshapes oatk's own gene calls
   into one tidy table and extracts per-gene FASTAs for the marker-gene
   set used by phylogeny. No new annotation is performed.
3. **[repeats](repeats/README.md)** — self-alignment + clustering to find
   repeat families and flag putative recombination-mediating repeats near
   contig ends (the mechanism behind multipartite mitochondrial genomes).
4. **[synteny](synteny/README.md)** — whole-contig `minimap2` alignments:
   automated mito-vs-plastid MTPT detection for every species, plus an
   opt-in pairwise mode for comparing species within a genus/family.
5. **[phylogeny](phylogeny/README.md)** — a partitioned ML tree per
   organelle from single-copy marker genes (mafft + IQ-TREE).
6. **[orf_scan](orf_scan/README.md)** — `ORFfinder` + full Pfam-A `hmmscan`
   on ORFs *not* already explained by oatk's core-gene annotation, to find
   transposable-element and mitovirus-derived content. This is
   meaningfully heavier than everything else (a multi-hour full-dataset
   run, not minutes) so it's **opt-in, not part of `run_all.sh`** — see
   its README for how to run it, including on the full dataset.
7. **[topology_plots](topology_plots/README.md)** — gene-labeled `Bandage`
   assembly-graph plots for visually triaging QC `fail`/`flag` species (no
   manual editing, unlike the reference figures in `mito_structural_variation`
   this reproduces). Cheap (~1s/render) — also **opt-in, not part of
   `run_all.sh`**, since its default target (QC-flagged species) depends on
   `qc_basic_stats` having already run.

## Running it

For a mini-project on a handful of species (put one species name per line
in a file):

```bash
analysis/run_all.sh --species-list my_species.txt
```

For the full dataset, wrap it in one LSF job (no step needs more than
modest single-node parallelism — this is deliberately not a
~1250-job array):

```bash
bsub -n 8 -q normal -o analysis_run.out -e analysis_run.err \
  analysis/run_all.sh --organelle both
```

Useful flags (all scripts under `src/` also accept these individually):

- `--species-list FILE` — restrict to these species (default: all)
- `--organelle mito|pltd|both` — default `both`
- `--qc-status pass[,flag[,fail]]` — which QC tiers to include downstream
  of the QC gate (default `pass`). Try `pass,flag,fail` as an exercise in
  "what does bad QC actually look like."
- `--force` — recompute everything, ignoring existing results
- `--skip-heavy` — stop after `annotation` (skip repeats/synteny/phylogeny)

Every step is incremental: species already present in a result table are
skipped unless their underlying assembly file is newer, or `--force` is
given. Re-running after a handful of species get fixed/added only
recomputes those species.

## `common/`

Shared infrastructure used by every module:

- `tool_paths.sh` — the one place every absolute tool path is defined.
  If a tool moves, this is the only file that needs to change.
- `species_discovery.py` / `.sh` — resolves each species' files
  robustly: filenames may or may not carry a `.k1001.s31.c<N>.` parameter
  infix, and a few species have more than one candidate run (see its
  docstring — naively preferring the newer/infixed file is actively wrong
  for at least one real case in this dataset).
- `io_utils.py` — atomic, incremental TSV read/write helpers.
- `minimap_utils.py` — shared PAF running/parsing, used by both
  `qc_basic_stats`' contamination pre-screen and `synteny`'s full MTPT
  module, so both agree on what "aligned length" means.
- `args.sh` — shared bash flag parsing (`--species-list`, `--organelle`,
  `--force`, `--qc-status`).

## A note on the existing `src/01-07_*.sh` runners

Those scripts (assembly generation and revision) are unchanged — this
suite only reads their output. Two things worth knowing:

- A handful of species have both a bare (`Species.mito.gfa`) and a
  parameter-infixed (`Species.k1001.s31.c80.mito.gfa`) run sitting side by
  side, because `03_move_assemblies.sh` doesn't guard against a re-run
  leaving old and new files together. `species_discovery.py` handles this
  by picking whichever run has the most complete, non-empty file set —
  not by naming or recency.
- If `analysis/run_all.sh` runs while an LSF revision batch
  (`07_revision_submit_oatk_from_list.sh`) is still in flight, it may see
  a partially-written species directory. This self-heals on the next
  re-run (the mtime-based skip logic will pick it up once the revision
  finishes) — a `fail`/`flag` for that one species mid-revision isn't a
  sign anything is broken.

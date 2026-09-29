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
qc_basic_stats  →  linearize         (opt-in, not in run_all.sh - see below)

denovo_annotation  →  ┬─ editing         ┬→ gff_export
                       └─ trans_splicing ┘   ↑
                       unitig_coords ────────┘
```
(the second block is the "deep annotation" suite — new, opt-in, not in
`run_all.sh`, and not yet run at full-dataset scale for every module; see
each README for validated smoke-test numbers before a full run)

1. **[qc_basic_stats](qc_basic_stats/README.md)** — is this assembly real
   and complete enough to trust? Per-species graph/contig stats, gene
   completeness, and a QC gate (`qc_summary.tsv`) that every module below
   filters against.
2. **[annotation](annotation/README.md)** — reshapes oatk's own gene calls
   into one tidy table and extracts per-gene FASTAs for the marker-gene
   set used by phylogeny. No new annotation is performed — see
   `denovo_annotation` below for a from-scratch alternative that's been
   shown to measurably out-perform oatk's own calls for several genes.
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
8. **[linearize](linearize/README.md)** — an alternative, independent
   linearization of every assembly graph via `gfatk linear` (cheap,
   graph-only, no new dependencies), plus an opt-in `gfatk resolve` fallback
   using real PacBio HiFi read-path evidence (via `GraphAligner`) for
   species where oatk's own Pathfinder produced nothing at all
   (`no_resolved_ctg_fasta`). Validated circuits get promoted back into
   `data/` (`05_promote_resolve.py`), guarded by a coverage gate and an
   automatic check that a promotion never silently discards a genuinely
   complete oatk Pathfinder genome (`07_restore_pathfinder.py` — built
   after exactly that happened once for 44 plastids; see its README).
   **Opt-in, not part of `run_all.sh`** — the `gfatk linear` pass is cheap
   enough to run whenever, but the read-evidence fallback needs LSF
   submission and touches raw read data outside this repo, so it's never
   run implicitly.

### Deep annotation suite (new)

These go beyond oatk's own targeted gene search — reconstructing gene
structure oatk's calling can't see on its own (trans-splicing, RNA
editing), and cross-checking oatk's calls against an independent,
from-scratch pipeline. All are opt-in, not in `run_all.sh`, and validated
so far on smoke-test subsets (~14-15 species) plus, for some, a
full-dataset run — see each README for exactly what's been validated at
what scale before trusting a full run yourself.

9. **[denovo_annotation](denovo_annotation/README.md)** — calls genes
   directly against the assembled contigs (`nhmmscan` + `tRNAscan-SE` +
   `barrnap`), independent of oatk's own internal calling. **Validated,
   not hypothesized**, to find real genes oatk's bundled database misses
   entirely (oatk's default mito family database is trained on
   Acrogymnospermae — gymnosperms — and this dataset is overwhelmingly
   not gymnosperms) — e.g. `nad5`: oatk finds 1 fragment/species, direct
   `nhmmscan` finds 4, consistently. Its `gene_calls.tsv` is a drop-in
   alternative anchor for `editing`/`trans_splicing` below
   (`--gene-calls`).
10. **[editing](editing/README.md)** — corrects gene boundaries for
    post-transcriptional C-to-U RNA editing using a purpose-built Rust
    tool (`orfedit`), including edits that create/destroy start/stop
    codons a DNA-level caller can't see. Validated: mitochondrial editing
    density several-fold higher than plastid, matching the literature's
    expected direction and rough magnitude.
11. **[trans_splicing](trans_splicing/README.md)** — reconstructs genes
    whose exons are transcribed separately and spliced at the RNA level
    (`nad1`/`nad2`/`nad5`/`rps3`, plus several cis-spliced multi-exon
    genes), using a purpose-built Rust tool (`transsplice`). Validated
    against 5 real GenBank reference mitogenomes and, independently,
    against `Arabidopsis_thaliana`'s own oatk-assembled genome in this
    dataset compared to its curated RefSeq protein — with honestly-reported
    partial results (some genes reconstruct well, some don't yet, and why).
12. **[unitig_coords](unitig_coords/README.md)** — maps every annotated
    feature back onto the raw assembly-graph unitigs, not just the final
    resolved contig. Needed for `gff_export`'s unitig-level GFF and for
    QC diagnostics (`qc_basic_stats/src/06_low_depth_paths.py`) that need
    to reason about the graph directly.
13. **[gff_export](gff_export/README.md)** — merges oatk's own calls,
    `denovo_annotation`, `editing`, and `trans_splicing` into one
    standard, hierarchical GFF3 per species (real `gene`/`mRNA`/`CDS`
    structure, not flat single-line features), with an explicit
    `annotation_tier` (`reconstructed`/`edited`/`raw`) per gene so you
    always know how much correction a call went through. Also emits a
    unitig-coordinate GFF alongside the contig-coordinate one.

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

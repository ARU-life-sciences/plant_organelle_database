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

## Setup

```bash
mamba env create -f environment.yml
conda activate plant_organellar_database
./install.sh   # this project's own Rust tools (gfatk, orfedit, transsplice, hmm_to_gff, filter_tblout)
```

`environment.yml` covers every third-party bioinformatics tool the
pipeline needs (versions pinned to what was actually validated, checked
via `<tool> --version` on the reference install, not copied from a doc).
`install.sh` builds this project's own Rust CLI tools via `cargo`, then
verifies `analysis/common/tool_paths.sh` resolves cleanly. Reference
databases (oatkDB's mito/plastid gene-family HMMs, `editing`/
`trans_splicing`'s profile sets) are separate from both — see each
module's README for what's built and what isn't yet (e.g.
`denovo_annotation/README.md`'s "Known gap: no plastid database yet").

## Analysis suite (`analysis/`)

The exercises below are still good starting points, but a lot of this is
now automated, incremental, and QC-gated across the whole dataset in
[`analysis/`](analysis/README.md). Two layers:

**Core suite** (run on `data/*.ctg.fasta` / oatk's own gene calls):

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
  element / mitovirus content.
- **topology_plots** — gene-labeled Bandage plots of the assembly graph,
  for visually triaging QC `fail`/`flag` species.
- **linearize** — an independent second opinion on path resolution
  (`gfatk linear`/`gfatk resolve`) for species oatk's own Pathfinder
  couldn't resolve at all, promoted back into `data/` once validated.

**Deep annotation suite** (new — reconstructs what oatk's targeted gene
search can't see on its own):

- **denovo_annotation** — calls genes directly (`nhmmscan` against a
  Viridiplantae-wide family database, `tRNAscan-SE`, `barrnap`),
  independent of oatk's assembler. Measurably finds more/better gene
  fragments than oatk's own bundled calls (oatk's default mito database is
  gymnosperm-trained; most of this dataset isn't gymnosperms) — see its
  README's validated comparison.
- **editing** — corrects gene boundaries for post-transcriptional C-to-U
  RNA editing (routine in plant mitochondria, sparser in plastids),
  including edits that create/destroy start/stop codons that a DNA-level
  caller can't see.
- **trans_splicing** — reconstructs genes whose exons are transcribed
  separately and spliced together at the RNA level (`nad1`/`nad2`/`nad5`/
  `rps3` trans-spliced; several more cis-spliced multi-exon genes),
  validated against real GenBank reference genomes.
- **unitig_coords** — maps every annotated feature back onto the raw
  assembly-graph unitigs, not just the resolved contig — needed wherever a
  gene straddles a repeat copy or a graph junction.
- **gff_export** — merges all of the above into one standard, hierarchical
  GFF3 per species, with an explicit `annotation_tier` (`reconstructed`/
  `edited`/`raw`) so you always know how much correction a given gene call
  went through.

Every module works on the full ~1250-species dataset by default, but also
takes a `--species-list` for a mini-project on a handful of species (a
genus, a family, whatever question you're asking). Run the lightweight
core-suite ones together with `analysis/run_all.sh` — see
[`analysis/README.md`](analysis/README.md) for usage. The deep annotation
suite is new, opt-in, and not yet run at full-dataset scale for every
module — see each one's README for validated smoke-test numbers before
trusting a full run.

**Before trusting any of it**: with ~1250 semi-automated assemblies, some
are empty, fragmented, or otherwise suspect. `analysis/qc_basic_stats`
computes a `pass`/`flag`/`fail` per species and every other module filters
against it by default — see
[`analysis/qc_basic_stats/README.md`](analysis/qc_basic_stats/README.md)
for what's actually being checked and why.

## Where to start: best data sources for example projects

Picking a species (or genus) to build an exercise around matters — some
are much better teaching examples than others, for specific reasons:

| project idea | start here | why |
|---|---|---|
| **A clean, fully-worked example** — trace one genome end to end (assembly graph → resolved genome → gene calls → GFF3) | `Acer_campestre` (mito) | The worked example in [`gff_export/README.md`](analysis/gff_export/README.md) — real `ccmB` RNA-editing calls and a real `nad2` trans-spliced gene model, already laid out line by line. |
| **Ground-truth validation** — "how do we know any of this is right?" | `Arabidopsis_thaliana` (mito) | The one species in this dataset with an independently curated RefSeq reference genome (NC_037304.1) to check against — used throughout `trans_splicing`/`editing`/`gff_export` as the actual accuracy benchmark, not a guess. 54/55 real gene features recovered; single-exon edited proteins at 99.5-100% identity. |
| **Multichromosomal/multipartite mitogenome structure** | `Dactylorhiza_fuchsii` (mito) | 14 real, independently-abundance-verified circular chromosomes — 9 of them match the published DTOL assembly (OZ487954-62) base-for-base. The wider orchid clade (`Ophrys_apifera`, `Orchis_mascula`, `Gymnadenia_conopsea`) shows the same 14-17-circle pattern, good for a comparative angle. |
| **Why assembly graphs are hard** — repeat-rich, fragmented, then resolved | `Azolla_filiculoides` (mito) | Started as 175 disconnected single-node fragments (oatk's Pathfinder found nothing joinable); `analysis/linearize`'s read-evidence fallback resolved a real 238-segment, 530kb genome. A dramatic, well-documented before/after (see `qc_basic_stats/README.md`'s "unjoined contigs" section and `linearize/README.md`). |
| **Trans-splicing / discontinuous genes** | any of `nad1`/`nad2`/`nad5`/`rps3` across several angiosperm species | `analysis/trans_splicing` reconstructs these directly; `reconstructed_junctions.tsv` tells you, per species, which junctions are genuinely cis vs. trans-spliced — not inherited from a reference topology. |
| **RNA editing density/biology** | any mito vs. plastid gene pair, e.g. `nad9` vs. `atpB` | `analysis/editing` gives real predicted edit counts/density per gene; mitochondrial editing runs several-fold denser than plastid — a good comparative exercise. |
| **Comparative genomics within a genus** | pick your own — every genus with 3+ sequenced species gets its own genus-relative QC baseline (`qc_basic_stats`'s size-outlier check) | Use `--species-list` on any module to restrict to one genus/family. |

For anything not covered here, `analysis/qc_basic_stats/results/qc_summary.tsv`
is the fastest way to browse: filter to `status=="pass"` for a safe
default, or deliberately pick `flag`/`fail` species (with their
`status_reasons`) if the *problem itself* is the teaching point.

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

#### Revision plan v2 (driven by `analysis/qc_basic_stats`)

The two recipes above predate `analysis/` and were based on
`meta/dead_ends.tsv` (mito-only, dead-ends and memory-kills only).
`src/08_revision_plan_v2.py` supersedes that planning step for anything
`analysis/qc_basic_stats`'s richer QC gate can see (both organelles, 6
distinct failure reasons) — it reads `qc_summary.tsv` and writes
`meta/revision_plan_v2.tsv`, one row per species with a specific
`--coverage`/`--max-gb` recipe and a **written justification**, not just a
number. Re-run it any time `qc_summary.tsv` changes.

Scope: every `fail`-tier species, plus `no_resolved_ctg_fasta` species
(flag-tier only under the current thresholds, but judged to need the same
treatment — see below). The other flag-only reasons (`non_circular`,
`size_outlier`, `elevated_cross_organelle_alignment` 15-35%) are
deliberately **not** included — those are mostly real biological
variation, not assembly defects, and reassembling on that basis alone
would be guessing rather than fixing anything.

| reason | species | recipe | why |
|---|---|---|---|
| `dead_end_nodes` | 39 | `-c 40 --max-gb 15` | established convention (above): lower coverage admits rarer k-mers that can bridge a dead-end path |
| `fragmented` (n_subgraphs≥10) | 21 | `-c 40 --max-gb 15` | same graph-connectivity mechanism as dead-ends. **Least validated mapping** — fragmentation could in principle also come from too much low-frequency *noise*, which would call for *raising* coverage instead. Check whether the pilot's fragmented-species actually improve before scaling this one to the rest. |
| `no_resolved_ctg_fasta` | 55 | `-c 40 --max-gb 15` | the graph exists but Pathfinder couldn't resolve a clean path through it; re-running Pathfinder alone on the *same* graph is deterministic and would just fail the same way again, so this needs the same graph-quality fix as dead-ends/fragmentation, not a cheaper shortcut |
| `missing_gfa` | 38 | `-c 80 --max-gb 16` (defaults) | assembly never completed; first attempt is just a normal run, no speculative adjustment. (If a memory-kill is found in that species' logs instead, use the memfix recipe above, not this one.) |
| `core_gene_pct` <50% | 4 | `-c 80 --max-gb 25` | low completeness looks like too little input data, not a graph-resolution problem, so raise the data budget rather than the coverage threshold |
| `high_cross_organelle_alignment` >35% | 101 | `-c 80 --max-gb 16` (defaults) | not obviously a connectivity/data problem — a baseline re-check to rule out "was the original assembly just noisy", not a targeted fix. Some of these are likely real biology (e.g. `Juncus_inflexus`, visually confirmed as a clean plastid graph despite this flag via `analysis/topology_plots`) — species still flagged after one attempt should be documented as probably genuine, not repeatedly reassembled on no particular basis. |

A species failing for more than one reason (46/258) gets exactly one
recipe, chosen by priority (most clearly mechanical first): `dead_end_nodes`
› `fragmented` › `no_resolved_ctg_fasta` › `missing_gfa` › `core_gene_pct`
› `high_cross_organelle_alignment`.

#### Promoting a rerun: never on file size/recency alone

`src/09_promote_revision.py` compares a completed `rev/out/` rerun against
whatever is currently in `data/` using the same metrics `qc_basic_stats`
uses (dead-end nodes, subgraph count, resolved-contig-fasta presence), and
only promotes it if it's **objectively better** — never because it's
newer or the file happens to be a different size.

This matters because a first attempt at automating this (species_discovery's
own "prefer the larger, non-empty candidate" tie-break, designed to reject
a truly-empty failed rerun) silently picked the **worse** of two genuinely
valid runs in many cases — fixing dead-ends/fragmentation by lowering
coverage often produces a cleaner *and smaller* assembly, which that
heuristic penalised. Concretely: of 123 already-completed reruns found
sitting unpromoted in `rev/out/` from before this pipeline existed, only
59/109 comparable mito results and 7/120 plastid results were actually
better - the other 38 mito / 6 plastid would have made things **worse**
if promoted on file size alone.

Superseded old files are moved to a `superseded/` subdirectory within the
species folder — never deleted, never left in place to compete (that
subdirectory is invisible to `species_discovery`'s non-recursive glob).

#### Pilot results (54 species, ~10 per category, first batch of new reassemblies)

The first pass at comparing `rev/out/` reruns against `data/` (both here
and in the historical-backlog promotion above) required BOTH dead-end
count and subgraph count to not worsen before calling a rerun "better".
Verified on real data that this is too strict: a rerun targeting
dead-ends can genuinely fix them while only slightly nudging subgraph
count, and the blanket rule mislabelled real fixes as "worse". Also, the
first version of this table was generated by comparing against a QC
snapshot taken *before* the historical-backlog promotion had run, so
several "test" species had already been fixed by that promotion and
weren't actually testing anything. `src/09_promote_revision.py` was
fixed to (a) take an explicit `--primary-metric` per category (improving
the metric the recipe actually targets is enough for "better", with a
new `mixed` outcome - flagged for manual review, never auto-promoted -
when the *other* metric regresses by more than a couple of units, which
does happen for real) and (b) always compare against a freshly-regenerated
`revision_plan_v2.tsv`. Corrected numbers:

| reason | mito | plastid |
|---|---|---|
| `no_resolved_ctg_fasta` | 0/9 better | **9/9 better** |
| `dead_end_nodes` | **1/9 better**, 1 `mixed`, 1 worse, 6 already-fixed (`same`) | 1/10 better |
| `fragmented` | **0/10 better**, 1 `mixed`, 2 worse, 7 unchanged | 0/10 better |
| `missing_gfa` | 1/4 comparable (6/10 still produced no assembly at all) | 0/7 better |
| `core_gene_pct` | 1/3 better | 0/3 better |
| `high_cross_organelle_alignment` | 0/10 better | 0/10 better |

Reading this, now that it's measuring the right thing:
- `no_resolved_ctg_fasta` still works essentially as predicted (plastid
  side, where the gap overwhelmingly is) - **worth scaling to the rest**.
- `dead_end_nodes` genuinely does work when there's still a real problem to
  fix (`Agrostis_gigantea`: 6→0 dead-ends) - the earlier "0/9" reading was
  an artifact of testing against already-fixed species, not a failure of
  the recipe. **Worth scaling**, but expect a `mixed` outcome sometimes
  (`Amsinckia_menziesii`: dead-ends 7→4, but subgraphs exploded 1→25 - a
  real trade-off, not noise) that needs a human to look at rather than
  auto-promote.
- `fragmented` is the one recipe that doesn't hold up as a general rule:
  7/10 mito species showed **zero change** in subgraph count at all, 2 got
  worse, and only one (`Clematis_viticella`, the original example that
  motivated this check) improved dramatically (10→1 subgraphs) - but even
  that picked up 3 new dead-ends, landing it in `mixed`, not a clean win.
  Lower coverage occasionally fixes fragmentation dramatically but usually
  does nothing - **not worth scaling as-is** to the remaining ~21 species;
  a different oatk parameter (e.g. `-a`/arc-coverage or `--max-bubble`,
  not `-c`/coverage) is probably needed for most of these, or they may
  need to stay documented as genuinely difficult assemblies.
- `high_cross_organelle_alignment` showing **zero** change in 10/10
  attempts, unaffected by any of these fixes, remains real evidence for
  "some of this is genuine biology, not an assembly defect", not just a
  guess.
- `missing_gfa` species mostly still fail to assemble even on retry with
  default parameters, suggesting a deeper data issue than a simple
  parameter tweak for those - worth checking raw read volume/quality
  before just retrying again.

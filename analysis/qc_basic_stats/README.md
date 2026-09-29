# qc_basic_stats

**The question:** is this assembly real and complete enough to trust? With
~1250 species assembled semi-automatically, some fraction will be empty,
fragmented, or otherwise suspect — before building an exercise (or a
phylogeny!) on top of one, it's worth knowing that up front rather than
discovering it later as an unexplained outlier.

This module computes per-species metrics and combines them into a single
QC gate that every other module filters against by default.

## What it computes

| script | output | what |
|---|---|---|
| `01_gfa_stats.py` | `gfa_stats.tsv` | `gfatk stats` per subgraph: node/edge counts, dead-ends, length, GC%, coverage. A genome can legitimately have >1 subgraph (multipartite mitochondria are real biology). |
| `02_contig_stats.py` | `contig_stats.tsv` | Per-contig length/circularity/GC% from `.ctg.fasta` headers + `seqkit`. Falls back to `gfatools gfa2fa` (unitig sequences, not a resolved path) for the ~5% of plastid assemblies missing a `.ctg.fasta`, labelled `contig_source=gfa_fallback`. |
| `03_gene_matrix.py` | `gene_matrix_{mito,pltd}.tsv` | Long-format `(species, gene) → hits/score/aligned length/n_contigs`, straight from oatk's own `.ctg.bed` gene calls. |
| `04_core_gene_list.py` | `core_genes_{mito,pltd}.tsv` | Per-gene prevalence across non-empty assemblies — **computed from the data, never a hardcoded gene list.** |
| `05_qc_summary.py` | `qc_summary.tsv` | Combines everything into one `pass`/`flag`/`fail` per species×organelle, with a `status_reasons` column explaining why. |
| `06_low_depth_paths.py` | `low_depth_paths.tsv` | Reports-only (never changes `status`): flags mito assembly graphs carrying low-depth side paths embedded in the graph, which can make a `gfatk` circuit walk the real genome several times over. See "Low-depth embedded paths" below. |

## Reading `qc_summary.tsv`

Each check is either a **flag** (worth a look, not disqualifying) or a
**fail** (don't trust this assembly by default):

- **missing/empty assembly** → fail
- **dead-end nodes** in the graph → flag if any, fail if >2
- **fragmentation** (`n_fragment_subgraphs`, disconnected components in
  the assembly graph **that aren't real chromosomes** — see
  "Multichromosomal mitogenomes" below) → flag if ≥5, fail if ≥10. The
  dataset shows a clean natural break: most mito species have 1-4
  fragment-subgraphs (consistent with real single-chromosome or modest
  multipartite biology), then a long tail with 5-20+. Confirmed against
  `topology_plots` output: e.g. `Clematis_viticella` (10 fragment
  subgraphs) visibly shows "lots of small linear segments" — a failed
  assembly, not previously caught since nothing checked this before.
- **non-circular** → flag only (many real plant mitochondria are
  legitimately multipartite/linear — this is informational, not a defect).
  Note: this reads the *resolved contig's* `circular=` flag from
  `contig_stats.tsv`, not `gfatk`'s raw-graph-level flag — verified these
  can disagree (a plastid's raw unitig graph can be "non-circular" while
  its resolved path still closes a genuine circle via a repeat, which is
  exactly the normal single-copy-repeat plastid structure).
- **genome size outlier** → flag if more than 4 MADs from the
  genus-level median (falling back to the dataset-wide median when a
  species has fewer than 3 sequenced congeners). Genus-relative because
  organelle genome size varies enormously across plant lineages — a fixed
  absolute cutoff would be meaningless. Size is the *resolved* contigs'
  total length where a resolved assembly exists, not the raw graph's — a
  handful of plastid graphs carry a second organism's plastid at 8-40x
  lower coverage (`Carduus_nutans`, `Cratoneuron_filicinum`,
  `Isopterygiella_pulchella`, `Malus_domestica` — real contamination, not
  a bug), which doubled their *graph* size but not their resolved genome.
- **core-gene completeness** below 80%/50% → flag/fail. The core-gene set
  is whichever genes are present in ≥90% of non-empty assemblies for that
  organelle (see `core_genes_*.tsv`) — not a textbook list, so it
  naturally adapts if oatk's gene family database changes.
- **no resolved contig fasta** (only the `gfa_fallback`, or an
  `unjoined`-tagged variant — see below) → flag.
- **cross-organelle alignment fraction** above 15%/35% → flag/fail, with
  one exemption — see below.

### Multichromosomal mitogenomes: a closed circle is a chromosome, not a fragment

Real plant mitogenomes aren't always one molecule — some lineages
routinely carry many small, independent circular chromosomes (*Silene
noctiflora*: 59-63 of them; Wu et al. 2015, PNAS 112:10185). The original
fragmentation check didn't distinguish these from a genuinely broken,
disconnected graph, and wrongly excluded 15 complete mitogenomes from
QC-`pass` (and therefore from every downstream annotation step). Fixed
2026-09-28: a subgraph now counts toward fragmentation *unless* it's
closed circular **and** its depth sits within `CHROMOSOME_DEPTH_BAND`
(3-fold) of the species' own length-weighted median subgraph depth — real
chromosome abundance varies roughly 2-fold, so a circle at plausible depth
is a chromosome; a circle at wildly different depth (e.g. a contaminant
plastid sitting in a mito graph, at plastid-level depth) still counts as
a fragment. `n_fragment_subgraphs` is the column the thresholds actually
apply to (`n_subgraphs` is kept as the raw count).

Verified directly, not just against theory: in this dataset, the orchids
`Ophrys_apifera`/`Orchis_mascula`/`Gymnadenia_conopsea`/`Dactylorhiza_fuchsii`
each have 14-17 closed single-unitig circles at 0.6-1.95x of each other's
depth, with genes spread 1-3 per circle and the same gene linkages
recurring across all four species — not noise. `Dactylorhiza_fuchsii`'s
14 circles match the published DTOL assembly (OZ487954-62) base-for-base
for 9 of them; the other 5 carry `ccmC`/`ccmFc`/`nad4`/`nad7`, genes the
published 9-chromosome assembly doesn't have at all.

### A complete plastid found inside a mito assembly isn't the plastid's fault

The cross-organelle contamination check (below) used to penalize *both*
organelles symmetrically when their assemblies overlapped. But a
genuinely complete, circular, core-gene-complete plastid genome that also
turns up (partially) inside a species' *mito* assembly is either real
plastid-derived sequence or the mito graph having picked up plastid
contigs — either way, it isn't evidence the plastid assembly itself is
wrong. Fixed: a plastid with `any_circular=true` and `core_gene_pct>=0.90`
is now exempt from this check entirely; only the mito side's own
plastid-like share is judged. Verified this matters at real scale: judging
mito on the larger of the two shares (the old logic) would have wrongly
demoted 193 otherwise-fine mitogenomes whose own plastid-like content was
only 0-15% (median 5.6%) — ordinary MTPT-scale transfer, not
contamination.

### Low-depth embedded paths (`06_low_depth_paths.py`, report-only)

Some mito assembly graphs carry small, low-depth side paths (well under
the ~2-fold abundance range of real chromosomes) embedded at both ends
into the main graph — found on `Silene_vulgaris`: ~206kb across 117 small
pieces at ~6% of genome depth, which caused a `gfatk` circuit to walk the
real 351kb genome 3-8 times through them (1.25Mb reported instead of the
real ~351kb). This script flags that pattern (embedded low-depth
sequence, by depth and topology, never by gene content alone —
gene-*free* sequence is not itself suspicious, since real multichromosomal
genomes carry gene-free chromosomes too) without changing `status` itself.
Species it flags get a `# linearisation uncertain` note in `gff_export`'s
output, recommending the unitig-level GFF (which doesn't depend on a path
choice) over the contig-level one.

### `no_resolved_ctg_fasta(unjoined)`: a non-empty `.ctg.fasta` isn't proof anything resolved

Found by checking a specific user question ("ferns definitely wouldn't
linearize properly, right?") against `Azolla_filiculoides`: its
`.mito.ctg.fasta` has 175 contigs, but **every one of them is `nv=1`** —
one contig per raw graph segment, meaning Pathfinder never joined a single
pair of nodes into a path. The file being non-empty made `has_ctg_fasta`
true, so this was invisible to the QC gate entirely.

Scanning the whole dataset for the same signature (all contigs `nv=1`,
≥5 contigs — the ≥5 floor matches the fragmentation check's own precedent
that real multipartite plant mitogenomes are 2-4 subgenomic circles, so
5+ never-joined pieces isn't plausible biology) found **78 species**, heavily
concentrated in ferns/pteridophytes (`Vandenboschia_speciosa`,
`Asplenium_marinum`, `Polystichum_lonchitis`, `Dryopteris_filix_mas`,
`Pteridium_aquilinum`, and most other ferns in the dataset) — consistent
with their reputation for large, repeat-rich mitogenomes that trip up
graph-topology heuristics. `has_ctg_fasta` now treats this case as *not*
resolved, same as a genuinely empty file, with `n_contigs`/
`pct_unjoined_contigs` columns added to `qc_summary.tsv` so the two
subtypes stay distinguishable (`no_resolved_ctg_fasta` vs.
`no_resolved_ctg_fasta(unjoined)`).

This had a second, non-obvious effect: the cross-organelle contamination
check was previously running its mito-vs-plastid `minimap2` alignment
against these 78 species' fragmented mito *junk* too (since it only
gated on the old, wrong `has_ctg_fasta`), producing spurious high
alignment fractions against noise — 20 species' **plastid** status
corrected from `flag`/`fail` back to `pass` once their mito side was
correctly excluded from that comparison. Their plastid genomes were never
actually contaminated; they were being compared against garbage.

See [`../linearize/README.md`](../linearize/README.md) for how these 78
species (plus the pre-existing `no_resolved_ctg_fasta` ones) get a second
attempt at linearization.

### Why the contamination check is a minimap2 alignment, not a gene-name check

An earlier version of this check looked for plastid-only gene names
turning up in a species' mito annotation (or vice versa) as a
contamination signal. That was verified empirically to never fire: oatk
only ever searches a bin against its own family HMM database
(`acrogymnospermae_mito.fam` / `angiosperm_pltd.fam`), so a plastid gene
name can never appear in a `.mito.ctg.bed` regardless of what's actually
in the sequence.

Instead, for every species with both organelles assembled, this runs one
`minimap2 -x asm5` alignment between its mito and plastid contigs and
checks what fraction of one organelle's assembly aligns to the other.
Genuine MTPT (mitochondrial-plastid DNA transfer) is a small, isolated
insert — a few percent at most. A large aligned fraction instead suggests
a mis-binned contig. This is a coarse pre-screen; see
[`../synteny/README.md`](../synteny/README.md) for the full per-species
MTPT picture (PAFs, dotplots, a proper summary table).

## Try this

Run with `--qc-status pass,flag,fail` instead of the default `pass` and
compare the annotation/repeats/synteny/phylogeny output — what changes,
and does it match what the QC reasons predicted? This is a genuinely good
exercise in "what does bad input data actually do to downstream analysis,"
not just an abstract warning.

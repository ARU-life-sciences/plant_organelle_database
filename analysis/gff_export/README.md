# gff_export

**The question:** `denovo_annotation` produces one combined GFF per
species, but it's flat - isolated single-line features, no `gene`/`mRNA`/
`CDS` hierarchy, and it doesn't reflect anything `editing` or
`trans_splicing` found (their results live in separate TSVs). There was
no single file you could open and see "here is `nad2`, here are its 5
exons, here is where it's trans-spliced, here is where it's RNA-edited".
This module builds that file: one standardized, GFF3-spec-correct,
hierarchical annotation per species, merging all three modules by
precedence.

## Precedence: which module's call wins, per gene

1. **`trans_splicing`** - if the gene is one of `nad1`/`nad2`/`nad5`/
   `rps3` and a reconstruction with >=1 filled exon exists (partial is
   still more informative than a raw fragment).
2. **`editing`** - if this gene has an RNA-editing-corrected call.
3. **`denovo_annotation`** - the raw best-scoring hit, otherwise (also
   the *only* tier for tRNA/rRNA genes - neither of the other two modules
   touches those gene classes).

Recorded per gene as `annotation_tier=reconstructed|edited|raw` - an
explicit confidence signal, not implicit in which script you happened to
run last.

## How GFF3 actually encodes this (confirmed, not assumed)

**Multi-exon / trans-spliced genes**: the spec's real mechanism is
multiple `CDS` rows sharing one `ID`, each `Parent=`-linked to one `mRNA`.
Nothing in the spec requires those rows to share a strand or even a
seqid - which is exactly what a trans-spliced gene needs, since a real
exon can be on the opposite strand from its neighbours (confirmed
directly in the 5 GenBank reference genomes `trans_splicing` fetched:
their `CDS` lines use nested `join(complement(a..b),c..d,...)` location
strings for precisely this reason - GFF3's equivalent is just multiple
same-`ID` rows instead of one nested string). One real consequence worth
knowing, not a bug: the parent `gene`/`mRNA` feature's own start/end must
span every exon's extremes, so a trans-spliced gene's `gene` row can
cover a huge region (e.g. `nad2` below spans 233kb) that's mostly
unrelated intervening sequence - inherent to representing trans-splicing
in GFF3 at all, and its strand is written as `.` (mixed/not applicable)
rather than a misleading single value.

**RNA editing sites**: no single certain, widely-adopted SO term was
available with full confidence, so this uses `sequence_alteration` (a
real, general term - "this position differs from a reference") as a
pragmatic, documented choice, not a claim that it's the one definitive
standard.

**Phase column** (GFF3 field 8, required for `CDS`): written as `0` for
every `CDS` here. For `editing`/`trans_splicing`-derived exons this is a
verified *guarantee*, not a computed-and-hoped-for value - both tools'
DP only ever emits codon-aligned boundaries by construction (see
`orfedit`/`transsplice`). For tier-3 raw `denovo_annotation`-only calls
it's an unverified *assumption* - those coordinates are HMMER envelope
bounds with no codon-boundary enforcement, so treat phase there as
best-effort, not confirmed.

## A real coordinate bug this caught

Building this surfaced an off-by-one between the two `gene_calls.tsv`-
shaped tables: `annotation/`'s (oatk's own `.ctg.bed`) is 0-based
half-open (BED-style), but `denovo_annotation/`'s was copying coordinates
straight from `nhmmscan`/`tRNAscan-SE`/`barrnap`'s native GFF3 output
(1-based inclusive) with no conversion - a silent 1bp inconsistency
between two tables documented as interchangeable. Fixed at the source
(`denovo_annotation/src/04_build_gene_calls.py`, `start - 1` at the one
place raw GFF coordinates enter that table) and regenerated. Given the
~60bp flank used everywhere downstream, this almost certainly didn't
change any of `denovo_annotation`'s already-reported comparison numbers,
but it needed fixing before coordinates here could be trusted at bp
resolution.

## Worked example: `Acer_campestre`, real output

```
gene   382417 383028 -  annotation_tier=edited        # ccmB
mRNA   382417 383028 -  n_exons=1;n_edits=11
CDS    382417 383028 -  ID=cds-ccmB
sequence_alteration  382457 ... (x11, one per real predicted edit site)

gene   442932 676371 .  annotation_tier=reconstructed # nad2
mRNA   442932 676371 .  n_exons=5;n_edits=8
CDS  442932-443090 +  exon_number=3
CDS  670221-670610 -  exon_number=2
CDS  671813-671965 -  exon_number=1
CDS  674100-674672 +  exon_number=4
CDS  676123-676371 +  exon_number=5

gene  rrn26  216102-219285 -  annotation_tier=raw
rRNA  216102-219285 -
```

`ccmB`'s 11 `sequence_alteration` children match exactly the 11 real
edits already documented for this species in `editing/README.md` - a
useful cross-check that the merge is reading the right rows, not just
that it runs without error. `nad2`'s exon strands (1,2 on `-`; 3,4,5 on
`+`) match `trans_splicing`'s own `reconstructed_junctions.tsv` for this
species exactly (cis 1->2, trans 2->3, trans 3->4, cis 4->5).

Verified structurally sound, not just eyeballed: every `Parent=` in every
one of the 14 smoke-test files resolves to a real `ID=` in the same file
(0 orphan references).

## Assembly provenance and unitig coordinates

Every ctg gets a `##sequence-region` pragma and a `region` feature
recording which linearisation produced it - `resolver=oatk_pathfinder`
or `resolver=gfatk_resolve` (the `linearize` fallback promoted into
`data/` by `05_promote_resolve.py`; read straight from the ctg.fasta
header) - plus GFF3's own `Is_circular=true` where applicable.

Using `unitig_coords`' placement map, every feature also records where it
sits on the raw assembly-graph unitigs:

```
ctg000001c ... CDS 442932 443090 + ... exon_number=3;unitig_loc=u72:1107-1265:-
ctg000001c ... gene 105511 107090 + ... unitig_span=u309,u1593
```

`unitig_loc` (1-based inclusive, strand relative to the unitig) when the
feature lies inside one unitig placement; `unitig_span` listing the
unitigs it crosses when it spans a junction between placements - flagged,
not guessed (`unitig_span=unplaced` if it sits on ctg sequence no
verified placement covers). A position inside a link overlap belongs to
two placements; the one where it's core is preferred.

`results/unitig/<species>.<organelle>.unitig.gff` is the same annotation
on the unitig sequences themselves (`unitig_coords/results/unitig_fasta/`).
A junction-crossing feature becomes several rows sharing one `ID` (GFF3's
discontinuous-feature mechanism, as for multi-exon CDS), each tagged
`unitig_split=<ctg>:<start>-<end>` with the ctg stretch it came from, and
CDS phase recomputed per piece. This view collapses repeat copies onto
one sequence - which is why the ctg-level GFF stays primary - so each
unitig's `region` row carries `n_ctg_copies` to keep that visible.

For `gfatk_resolve` contigs there is no recorded unitig walk, so
placements come from exact sequence matching; a match found only after
trimming the unitig's ends (typically one of several near-identical
repeat variants) claims only its verified middle.

Verified on a 7-species test set spanning both resolvers, linear and
circular contigs, cross-contig trans-splicing and an origin-wrapping
single-unitig genome: every `unitig_loc` and every split piece is
byte-identical between ctg and unitig sequence (bar SNP/1bp-indel
differences confined to link overlaps, where the ctg carries one
neighbour's copy), every split feature's pieces sum to its full length,
and every `Parent=` resolves in both files.

## Two merge bugs fixed (2026-09-28)

- **Cross-organelle leakage.** Editing calls, edits and
  `trans_splicing` reconstructions were matched on species + gene only.
  Genes present in both genomes (`atpA`, `rpl16`, `rps19`, `rps4`,
  `rps14`, ...) could pick up the *other* organelle's call - 3911 edited
  gene models in 673 species had the other genome's coordinates and edit
  sites, and plastid `rps3` could take the mito trans-splicing
  reconstruction. Every table is now split by species AND organelle
  before lookup (with organelle matched, all 54278 editing calls
  correspond exactly to their raw best hit).
- **Cross-contig trans-spliced genes.** 258 of 3558 reconstructions have
  exons on more than one contig; the `gene`/`mRNA` row used the first
  exon's contig but spanned every exon's coordinates, running off the
  end of the contig. These now get one `gene`/`mRNA` row per contig
  (shared `ID`), and each edit site sits on its own exon's contig.

## Running

Stage 4 of `denovo_annotation/work/run_full_suite.sh`; not in
`run_all.sh`. Depends on `denovo_annotation`'s `gene_calls.tsv` and
benefits from `editing`/`trans_splicing` having been run first, but
degrades gracefully to tier-3-only output if they haven't (every
upstream table is loaded with its expected columns even when the file
doesn't exist yet). The unitig layer needs `unitig_coords` run first and
is skipped per species when its map is absent.

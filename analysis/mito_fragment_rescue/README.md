# mito_fragment_rescue

Standalone Rust tool (not wired into `run_all.sh` / the promotion pipeline).
Builds offline with `cargo build --release` (no external crates).

## Problem it addresses

oatk classifies each unitig as mito OR plastid by comparing the *total* HMM
gene score for each organelle on that unitig (winner-take-all; see the
revision investigation this session). For a `missing_gfa` species, a unitig
can carry real, genuine mitochondrial gene hits alongside a much stronger
plastid score, so the whole unitig - including the real mito content - gets
classified as plastid and excluded from `mito.gfa` entirely. The mito
content is still sitting in oatk's own raw nhmmscan output
(`*.annot_mito.txt`), just never extracted.

Several gene families are genuinely homologous between plant mito and
plastid genomes (shared bacterial ancestry: ribosomal proteins, ATP synthase
subunits, `nad`~`ndh`, `cob`~`petB`), so a naive "any mito HMM hit counts"
approach produces false positives from HMM cross-reactivity against the
plastid's own genes, not real mito signal.

## Method

For each candidate species, for every historical oatk rerun that has:
- a non-empty `<prefix>.annot_mito.txt`
- a same-prefix `<prefix>.annot_pltd.txt`
- a same-prefix `<prefix>.utg.final.gfa` (preferred) or `<prefix>.utg.gfa`

...it parses mito hits restricted to a curated **mito-exclusive gene list**
(genes with no plastid-encoded paralog at all: `nad1-9`, `cox1-3`, `cob`,
`mttB`, `matR`, `ccmB/C/Fc/Fn`, `sdh3/4`), e-value <= 1e-6, then discards any
hit whose coordinates overlap (50bp slop) a same-run plastid HMM hit on the
same unitig. The cross-check is always against the *same oatk run's own*
plastid annotation, not a separately-run "official" plastid call - unitig
IDs are just per-run sequential labels, so comparing across different
assembly graphs is not meaningful even if an ID happens to match.

It then picks, per species, the single candidate run with the most distinct
clean genes, and validates that every kept hit's coordinates actually fit
inside its unitig's sequence length in that run's graph file. **This
validation matters**: at least two species (`Sagina_nodosa`,
`Myosotis_sylvatica`) have a `*.utg.gfa` whose unitig numbering/sequences no
longer correspond to the graph their own `annot_mito.txt` was computed
against (likely overwritten by a later rerun under the same filename). Any
candidate that fails this check is discarded outright rather than risk
silently extracting the wrong sequence - this is a real pipeline
data-integrity hazard (stale graph file reused under the same prefix after
a rerun), not just a theoretical concern.

Surviving hits are merged into fragments per (gene, unitig) with the same
50bp slop (keeping distant HSPs of the same gene on the same unitig as
separate fragments - this does not stitch multi-exon genes together), then
extracted directly from the GFA's `S` line sequence (reverse-complemented
if the alignment is on the minus strand).

## Output

`results/mito_fragments.{fasta,tsv}` - every fragment is explicitly labelled
`fragmentary_partial_recovery;not_a_resolved_genome` in both the FASTA
header and the TSV `note` column. **This is not a `mito.gfa` and nothing
here should be promoted as a resolved mitochondrial genome** - it is
gene-level evidence recovered from an otherwise plastid-classified unitig,
for manual review.

## Current result (all 24 species currently `missing_gfa` for mito)

Run against the full, current population - not a hand-picked subset. The
candidate list is every species with `organelle=mito, status=fail,
status_reasons~=missing_gfa` in `analysis/qc_basic_stats/results/qc_summary.tsv`
(regenerate it with the command below; two species that were historically on
the older static list, `Carex_atrofusca` and `Ranunculus_circinatus`, have
since been independently promoted via the revision pipeline and are excluded
since they're no longer `missing_gfa`).

```
awk -F'\t' '$2=="mito" && $14=="fail" && $15 ~ /missing_gfa/{print $1}' \
    analysis/qc_basic_stats/results/qc_summary.tsv | sort
```

Only 4 of the 24 produce any output; the rest come back empty (0 genes, 0
fragments) once validated against same-run plastid calls and graph-sequence
consistency:

| species | genes | fragments | bp |
|---|---|---|---|
| Sagina_nodosa | 18 | 31 | 15386 |
| Silene_flos_cuculi | 19 | 32 | 14955 |
| Hydrocharis_morsus_ranae | 4 | 4 | 1775 |
| Myosotis_sylvatica | 1 | 1 | 72 |

Sagina_nodosa and Silene_flos_cuculi have strong, reproducible signal (same
gene clusters recur across independent coverage reruns). Three further
species (`Gagea_bohemica`, `Gagea_lutea`, `Lilium_martagon`,
`Tulipa_sylvestris`) have no eligible rerun at all (no saved annotation +
matching graph file pair under `rev/out/`) and are reported as skipped.

## Usage

```
cargo build --release
./target/release/mito_fragment_rescue --root <repo_root> [SPECIES...] \
    [--species-list FILE] [--evalue 1e-6] [--slop 50] [--out-dir DIR]
```

With no species given, defaults to the 7 species identified during initial
investigation (a subset of the 24 above). For a full, current run, generate
a species list with the `awk` command above and pass it via `--species-list`.

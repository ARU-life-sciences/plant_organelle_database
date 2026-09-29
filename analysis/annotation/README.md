# annotation

**The question:** what genes does this genome have, and where? oatk
already answers this during assembly (that's what `.ctg.bed` and
`.annot_{mito,pltd}.txt` are) — this module doesn't re-annotate anything,
it just reshapes oatk's own output into a form the other modules (and
students) can use directly.

If you're coming from the `mito_structural_variation` repo: that one runs
a full from-scratch pipeline (ORFfinder → hmmscan → tRNAscan-SE → barrnap
→ Infernal) to annotate genes. This module doesn't reimplement that — it
just reshapes oatk's own calls.

**Update, no longer true that oatk's own calling is "good enough"**: see
[`../denovo_annotation/README.md`](../denovo_annotation/README.md), which
runs that same from-scratch approach directly against this dataset and
**measurably out-performs** oatk's bundled calls for several genes (oatk's
default mito family database is trained on Acrogymnospermae — gymnosperms
— and this dataset is overwhelmingly not gymnosperms). Concretely: `nad5`
oatk finds 1 fragment/species vs. 4 via direct `nhmmscan`; some QC-`pass`
species have zero oatk calls for `nad1`/`nad2`/`nad5`/`rps3` entirely
where direct calling finds 4-6 fragments each. This module (`annotation/`)
is still the right place for a quick, zero-extra-compute reshape of
what's already in `data/`; `denovo_annotation/` is the better source when
gene completeness itself matters.

## What it computes

- `01_gene_table.py` → `gene_calls.tsv`: every `.ctg.bed` gene call across
  QC-eligible species, concatenated with species/organelle columns added.
  `.ctg.bed` is already GFF-shaped (`seq_name align_from align_to
  gene_name score strand`), so this is a direct reshape.
- `02_gene_fasta_extract.py` → `genes/{mito,pltd}/<gene>.fasta`: one FASTA
  per phylogeny marker gene (see [`qc_basic_stats`](../qc_basic_stats/README.md)'s
  `core_genes_*.tsv` — present in ≥90% of species and single-copy in ≥80%
  of those; single-copy threshold is 80%, not 90%, because at full-dataset
  scale (~1200 species) the best mito genes top out around 89% single-copy —
  real biology, plant mitochondrial genomes duplicate genes across
  multipartite structure far more than plastids do, not a bug), one
  sequence per species, extracted via `samtools faidx` and strand-corrected.

## Multi-copy / multi-exon genes

Some real genes (e.g. plant mitochondrial `nad1/2/4/5/7`, `cox2`,
`ccmFc/Fn`, `rps3` — trans-spliced or multi-exon in plants) show up as
several separate `.ctg.bed` hits per species. Extracting "the gene" as a
single contiguous FASTA record for these would silently concatenate
unrelated exon fragments, which is wrong. Rather than guess at exon
boundaries, this module simply doesn't extract per-gene FASTAs for genes
below the single-copy threshold — see
[`../phylogeny/README.md`](../phylogeny/README.md) for exactly which genes
that excludes and why it matters for tree-building.

**Update**: [`../trans_splicing/README.md`](../trans_splicing/README.md)
now does the thing this section says not to guess at — reconstructing
these genes properly (per-species cis/trans junction classification, not
inherited from a reference topology), validated against real GenBank
references. It hasn't been fed back into this module or `phylogeny`'s
marker-gene set yet (see that README's "Try this" — a real next step, not
done in this phase), so the exclusion here still stands for now.

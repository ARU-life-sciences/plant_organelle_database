# unitig_coords

**The question:** a resolved `.ctg.fasta` genome doesn't record which raw
GFA unitigs (in what order, orientation, and copy number) built it —
that information only exists implicitly in the graph. Anything that wants
to reason about the graph directly (which repeat copy is a feature really
on, has a low-depth side path been walked through more than once, does a
gene sit inside a link overlap) needs that placement map. This module
builds it: an unambiguous mapping from every unitig onto the coordinates
of the `.ctg.fasta` it contributed to.

Downstream of [`linearize`](../linearize/README.md) (produces the
`.ctg.fasta`/`.gfa` this module reads) and gated by
[`qc_basic_stats`](../qc_basic_stats/README.md) (`--qc-status`, default
`pass`). Translation runs **ctg → unitig only, deliberately**: existing
ctg-coordinate results (gene calls, etc.) get relabelled into unitig
coordinates using this map, rather than independently re-calling genes at
unitig level, since a feature spanning a join between two placed unitigs
would otherwise fragment, and a multi-copy unitig would make "which copy"
ambiguous. See [`gff_export`](../gff_export/README.md) for the consumer —
its unitig-level GFF is built this way.

## What it computes

| script | output | what |
|---|---|---|
| `00_build_unitig_map.py` | `results/unitig_map/<species>.<organelle>.tsv`, `results/unitig_fasta/<species>.<organelle>.unitig.fasta` | per-unitig placement onto ctg coordinates, plus the raw (non-linearised) unitig sequences themselves |
| `01_linearization_qc.py` | `results/linearization_qc.tsv` | two independent completeness signals (graph-side, ctg-side) for how fully the linearisation accounts for the graph it came from |

## Method

### `00_build_unitig_map.py`: two tiers, by how the ctg was built

**Tier 1 — oatk's own Pathfinder** (the majority of species). The
`.ctg.fasta` header already records the answer: `path=u66+,u68+,u69+,...`.
The script walks that path and subtracts each junction's overlap (read
straight from the GFA's `L` lines, e.g. `L u66 + u68 + 1304M` = 1304bp) to
reconstruct exact ctg coordinates for every placement. Verified by hand on
`Arabidopsis_thaliana` mito: its 8 segments (371,177bp total) minus 7
internal + 1 circular-closing overlap (2,349bp total) equals exactly
368,828bp — the header's own `length=`, bit for bit. Correctly handles
reverse-complement placement and genuinely multi-copy unitigs (a unitig
appearing twice in the path — real recombination-mediated repeat
structure, not noise; an earlier substring-search-only version of this
script missed one of two real placements for a repeated unitig, which the
`path=` field catches).

**Tier 2 — `gfatk resolve` fallback** (species oatk's Pathfinder couldn't
resolve, promoted by `linearize/05_promote_resolve.py`). These contigs'
headers don't carry an ordered/oriented walk (`resolve_summary.tsv` only
records a segment *count*), so there's no authoritative path to parse.
Falls back to exact sequence matching: each unitig's raw sequence (or its
reverse complement) either appears in the ctg byte-exact, or it doesn't,
with circularity handled by searching the ctg doubled on itself and a
small symmetric trim tolerance at each end for the overlap a circular
closure removes. Confirmed on `Achillea_maritima` mito: unitig `u6`'s raw
202,131bp is exactly 36bp longer than the trimmed circular ctg's
202,095bp. A Tier-1 species whose path turns out to be unwalkable (a
missing segment/link) drops to this tier rather than emit a partial map —
logged per-species, and counted separately in the run summary so it's
never silent.

**On-the-origin wraparound** (both tiers): on a circular ctg, a placement
can wrap the origin, so `ctg_end` can exceed the ctg's own length — this
is *expected*, not a bug (Tier 1: the first unitig of every circular path,
since oatk trims the closing overlap off its start, shifting the ctg
origin partway into it; Tier 2: any seqmatch hit across the join).
Consumers should take positions modulo ctg length.

**Overlap regions carry one neighbour's copy, and it can differ by 1bp**:
adjacent Tier-1 placements overlap by their link's bp, and within that
overlap the ctg only stores one neighbour's version of the sequence —
verified at full dataset scale to differ from the other copy by a SNP or
1bp indel in 17/1714 path-mapped species, always confined to the overlap
itself (every placement's non-overlap *core* was byte-exact). A ctg
position inside a link overlap belongs to two placements and is only
approximately placed in either — prefer whichever placement it's core to.

### `01_linearization_qc.py`: two signals, kept separate on purpose

- **Graph-side**: fraction of total *unique* unitig content (each unitig
  counted once regardless of copy number) placed anywhere in the ctg.
  Zero-placement unitigs are graph content the linearisation dropped
  entirely.
- **Ctg-side**: fraction of the ctg's own length accounted for by a
  placed unitig — should be near 100% when placements tile without
  overlap; low values flag either untrimmed overlap (not modelled by the
  exact-match approach) or real ctg sequence that didn't come from any
  whole unitig.

`n_multicopy_unitigs` is reported as its own column, deliberately not
folded into either score: a unitig placed more than once is real
recombination/dosage signal, not a defect — exactly why the ctg-level GFF
stays `gff_export`'s primary output (a unitig-level GFF collapses both
copies onto one sequence, losing that dosage information; see its
README's "Assembly provenance and unitig coordinates" section). The
`qc_basic_stats`-derived columns (`qc_status`/`any_circular`/
`core_gene_pct`/`n_contigs`) are included as corroborating context only —
not blended into one score, since each measures something genuinely
different and conflating them would hide which one actually explains a
low value.

## Reading the output

`unitig_map/<species>.<organelle>.tsv` — one row per placement (a
multi-copy unitig gets >1 row; an unplaced unitig gets none — the
`unitig_fasta/` file is the complete list of what exists, which is why
that's written for every species regardless of placement success):

```
unitig  unitig_length  ctg_seqid  ctg_start  ctg_end  strand  source  match_trim
```

`ctg_start`/`ctg_end` are 0-based half-open, matching `.ctg.bed`'s
convention elsewhere in this pipeline. `source` is `path` (Tier 1) or
`seqmatch` (Tier 2). `match_trim` (seqmatch only, 0 for `path`) is how
many bp at each end weren't independently verified — only
`[ctg_start+trim, ctg_end-trim)` is trustworthy for a trimmed hit.

`linearization_qc.tsv` — one row per species×organelle:

```
species organelle map_source ctg_length_bp n_unitigs total_unitig_bp
pct_unitig_content_placed n_orphan_unitigs n_multicopy_unitigs
pct_ctg_explained_by_unitigs qc_status any_circular core_gene_pct n_contigs
```

`map_source` is the sorted, comma-joined set of sources seen for that
species (`path`, `seqmatch`, or both — mixed is possible when some
unitigs needed the seqmatch fallback), or `none` if nothing placed.

## Running

```bash
python3 analysis/unitig_coords/src/00_build_unitig_map.py --organelle both --jobs 8
python3 analysis/unitig_coords/src/01_linearization_qc.py --organelle both
```

Both accept `--species-list`/`--qc-status` (default `pass`, filtered
against `qc_basic_stats/results/qc_summary.tsv` — if that file is
missing, `00_build_unitig_map.py` warns and proceeds unfiltered rather
than silently processing nothing). `00_build_unitig_map.py` also takes
`--jobs` (thread pool size) and `--force`; a map is otherwise rebuilt
automatically — without `--force` — whenever it's older than its GFA or
`.ctg.fasta` inputs, since a stale map (e.g. after
`linearize/07_restore_pathfinder.py` swaps a genome in place) is worse
than no map. `01_linearization_qc.py` has neither flag — it's a cheap
single-pass summary, fine to always rerun in full.

Both scripts are pure Python (stdlib only) and don't call any external
tool via `tool_paths.sh` — the only dependency is
`analysis/common/species_discovery.py` for per-species file resolution.

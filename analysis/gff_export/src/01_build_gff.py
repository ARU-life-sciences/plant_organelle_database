#!/usr/bin/env python3
"""Build one standardized, hierarchical GFF3 per species, merging the
three annotation modules by precedence:

1. `trans_splicing` reconstruction (`nad1`/`nad2`/`nad5`/`rps3` only, if a
   reconstruction with >=1 filled exon exists - partial is still more
   informative than a single raw fragment).
2. `editing`'s RNA-editing-corrected single-exon call, if this gene has one.
3. `denovo_annotation`'s raw best-scoring hit, otherwise (also the only
   tier for tRNA/rRNA genes - neither of the other two modules touches
   those gene classes).

Which tier won is recorded per gene as `annotation_tier=reconstructed|
edited|raw` - an explicit confidence signal, not just implicit in which
script you happened to run.

GFF3 mechanics used here (see `../README.md` for the full explanation):
- Protein-coding genes are `gene` -> `mRNA` -> one or more `CDS` rows
  sharing one `ID` (that's the spec's actual mechanism for a discontinuous/
  multi-exon feature - it does not require those rows to share a strand or
  seqid, which is exactly what a trans-spliced gene needs).
- tRNA/rRNA are `gene` -> `tRNA`/`rRNA` directly - no `CDS`/`mRNA` layer,
  they're not protein-coding.
- Predicted RNA-editing sites are child `sequence_alteration` features
  under the relevant `CDS` `ID` (a real, general Sequence Ontology term
  for "this position differs from a reference" - a pragmatic choice, not
  a claim that this is THE canonical way to encode RNA editing in GFF3).
- `phase` (GFF3 column 8) is written as `0` for every `CDS` - a verified
  guarantee for `editing`/`trans_splicing`-derived exons (their DP only
  ever emits codon-aligned boundaries by construction) but an unverified
  *assumption* for tier-3 raw `denovo_annotation`-only calls (HMMER
  envelope coordinates aren't codon-boundary-enforced) - see README.

Assembly provenance and unitig coordinates (see README):
- Each ctg gets a `##sequence-region` pragma and a `region` feature
  carrying `resolver=oatk_pathfinder|gfatk_resolve` (which linearisation
  produced it - read from the ctg.fasta header) and GFF3's `Is_circular`.
- Every feature gets `unitig_loc=<unitig>:<start>-<end>:<strand>` (1-based
  inclusive, on the raw unitig from `unitig_coords`) when it lies inside
  one unitig placement, or `unitig_span=<u>,<u>,...` when it crosses a
  junction between placements - flagged explicitly, not guessed. Needs
  `unitig_coords/results/unitig_map/`; skipped per species if absent.
- A second, unitig-level GFF places the same features on the unitig
  sequences themselves. A junction-crossing feature is split into pieces
  sharing one `ID` (GFF3's discontinuous-feature mechanism), CDS phase
  recomputed per piece. This view collapses repeat copies onto one
  sequence - the ctg-level GFF stays primary; each unitig's `region` row
  records `n_ctg_copies` so that dosage isn't silently lost.

Output:
  analysis/gff_export/results/<species>.<organelle>.gff (ctg-level, primary)
  analysis/gff_export/results/unitig/<species>.<organelle>.unitig.gff
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = ANALYSIS_DIR.parent
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import species_discovery as sd  # noqa: E402

TRANS_SPLICING_GENES = {"nad1", "nad2", "nad5", "rps3"}
UNITIG_MAP_DIR = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_map"
RESOLVER_RE = re.compile(r"\bresolver=(\S+)")
CIRCULAR_RE = re.compile(r"\bcircular=(true|false)")
FLIP_STRAND = {"+": "-", "-": "+", ".": "."}


@dataclass
class Segment:
    """One linear stretch of a unitig placement on a ctg (a placement that
    wraps a circular ctg's origin is two segments)."""
    unitig: str
    unitig_length: int
    strand: str
    ctg_lo: int
    ctg_hi: int
    unitig_off: int  # forward-walk offset into the (oriented) unitig at ctg_lo


@dataclass
class SeqContext:
    lengths: dict[str, int] = field(default_factory=dict)
    resolver: dict[str, str] = field(default_factory=dict)
    circular: dict[str, bool] = field(default_factory=dict)
    segments: dict[str, list[Segment]] = field(default_factory=dict)  # by ctg, sorted
    tiles: dict[str, list[Segment]] = field(default_factory=dict)     # disjoint partition
    unitig_lengths: dict[str, int] = field(default_factory=dict)
    unitig_copies: dict[str, int] = field(default_factory=dict)


def load_seq_context(ctg_fasta: Path, map_path: Path) -> SeqContext:
    ctx = SeqContext()
    name = None
    with open(ctg_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                name = header.split()[0]
                ctx.lengths[name] = 0
                m = RESOLVER_RE.search(header)
                ctx.resolver[name] = m.group(1) if m else ("oatk_pathfinder" if "path=" in header else "unknown")
                c = CIRCULAR_RE.search(header)
                ctx.circular[name] = c is not None and c.group(1) == "true"
            elif name is not None:
                ctx.lengths[name] += len(line.strip())
    if not map_path.exists():
        return ctx

    for row in pd.read_csv(map_path, sep="\t").itertuples():
        clen = ctx.lengths.get(row.ctg_seqid)
        if not clen:
            continue
        # A trimmed seqmatch hit (near-identical repeat variant) only
        # verifies its middle, so only that is claimed - an unverified end
        # would place features on sequence the unitig doesn't carry.
        trim = int(getattr(row, "match_trim", 0) or 0)
        a, b = int(row.ctg_start) + trim, int(row.ctg_end) - trim
        segs = ctx.segments.setdefault(row.ctg_seqid, [])
        if a >= clen:  # the verified core starts past the origin
            a, b, off = a - clen, b - clen, trim
        else:
            off = trim
        segs.append(Segment(row.unitig, int(row.unitig_length), row.strand, a, min(b, clen), off))
        if b > clen:  # wraps the origin (see unitig_coords/src/00_build_unitig_map.py)
            segs.append(Segment(row.unitig, int(row.unitig_length), row.strand, 0, b - clen, off + clen - a))
        ctx.unitig_lengths[row.unitig] = int(row.unitig_length)
        ctx.unitig_copies[row.unitig] = ctx.unitig_copies.get(row.unitig, 0) + 1

    # Adjacent placements overlap by their GFA link; tiles hand each
    # overlap to the earlier placement so a split feature's pieces never
    # double-count bases.
    for ctg, segs in ctx.segments.items():
        segs.sort(key=lambda s: (s.ctg_lo, s.ctg_hi))
        tiles, covered_to = [], 0
        for s in segs:
            lo = max(s.ctg_lo, covered_to)
            if lo < s.ctg_hi:
                tiles.append(Segment(s.unitig, s.unitig_length, s.strand, lo, s.ctg_hi,
                                     s.unitig_off + (lo - s.ctg_lo)))
                covered_to = s.ctg_hi
        ctx.tiles[ctg] = tiles
    return ctx


def to_unitig(seg: Segment, s: int, e: int, strand: str) -> tuple[int, int, str]:
    """ctg [s,e) inside `seg` -> (unitig start, unitig end, strand), 0-based
    half-open on the unitig's own forward sequence."""
    o_s = seg.unitig_off + (s - seg.ctg_lo)
    o_e = o_s + (e - s)
    if seg.strand == "+":
        return o_s, o_e, strand
    return seg.unitig_length - o_e, seg.unitig_length - o_s, FLIP_STRAND[strand]


def locate(ctx: SeqContext, ctg: str, s: int, e: int) -> tuple[Segment | None, list[Segment]]:
    """(containing segment or None, tiles the feature crosses in ctg order).
    A containing segment clear of every other placement (its core) is
    preferred over one where the feature sits inside a link overlap."""
    segs = ctx.segments.get(ctg, [])
    containing = [g for g in segs if g.ctg_lo <= s and e <= g.ctg_hi]
    if containing:
        def in_core(g):
            return not any(o is not g and o.ctg_lo < e and s < o.ctg_hi for o in segs)
        return next((g for g in containing if in_core(g)), containing[0]), []
    return None, [t for t in ctx.tiles.get(ctg, []) if t.ctg_lo < e and s < t.ctg_hi]


def add_unitig_layer(lines: list[str], ctx: SeqContext) -> tuple[list[str], list[str]]:
    """Returns (ctg-level lines with unitig attributes added, unitig-level lines)."""
    ctg_out, unitig_out = [], []
    for line in lines:
        f = line.split("\t")
        ctg, s, e, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        if ctg not in ctx.segments:
            ctg_out.append(line)
            continue
        seg, crossed = locate(ctx, ctg, s, e)
        if seg is not None:
            us, ue, ustrand = to_unitig(seg, s, e, strand)
            f[8] += f";unitig_loc={seg.unitig}:{us + 1}-{ue}:{ustrand}"
            unitig_out.append("\t".join([seg.unitig, f[1], f[2], str(us + 1), str(ue), f[5], ustrand, f[7], f[8]]))
        else:
            f[8] += f";unitig_span={','.join(t.unitig for t in crossed) or 'unplaced'}"
            for t in crossed:
                ps, pe = max(s, t.ctg_lo), min(e, t.ctg_hi)
                us, ue, ustrand = to_unitig(t, ps, pe, strand)
                phase = f[7]
                if f[2] == "CDS" and phase != ".":
                    # bases of this CDS upstream of the piece, in transcript direction
                    k = (ps - s if strand != "-" else e - pe) - int(phase)
                    phase = str((3 - k % 3) % 3)
                unitig_out.append("\t".join([t.unitig, f[1], f[2], str(us + 1), str(ue), f[5], ustrand, phase,
                                             f[8] + f";unitig_split={ctg}:{ps + 1}-{pe}"]))
        ctg_out.append("\t".join(f))
    return ctg_out, unitig_out


def ctg_header(ctx: SeqContext, seqids: list[str]) -> tuple[list[str], list[str]]:
    pragmas = [f"##sequence-region {c} 1 {ctx.lengths[c]}" for c in seqids]
    regions = ["\t".join([c, "gff_export", "region", "1", str(ctx.lengths[c]), ".", "+", ".",
                          gff_attrs({"ID": f"region-{c}", "resolver": ctx.resolver[c],
                                     "Is_circular": "true" if ctx.circular[c] else None})]) for c in seqids]
    return pragmas, regions


def unitig_header(ctx: SeqContext, unitigs: list[str]) -> tuple[list[str], list[str]]:
    pragmas = [f"##sequence-region {u} 1 {ctx.unitig_lengths[u]}" for u in unitigs]
    regions = ["\t".join([u, "gff_export", "region", "1", str(ctx.unitig_lengths[u]), ".", "+", ".",
                          gff_attrs({"ID": f"region-{u}", "n_ctg_copies": ctx.unitig_copies[u]})])
               for u in unitigs]
    return pragmas, regions


def assemble(pragmas: list[str], regions: list[str], features: list[str]) -> str:
    """GFF3 text with each seqid's region row leading its features."""
    by_seqid: dict[str, list[str]] = defaultdict(list)
    for line in features:
        by_seqid[line.split("\t", 1)[0]].append(line)
    body = []
    for region in regions:
        seqid = region.split("\t", 1)[0]
        body.append(region)
        body.extend(sorted(by_seqid.pop(seqid, []), key=lambda l: int(l.split("\t")[3])))
    for seqid in sorted(by_seqid):  # features on a seqid with no region row (shouldn't happen)
        body.extend(by_seqid[seqid])
    return "\n".join(["##gff-version 3", *pragmas, *body]) + "\n"


def best_hit_per_gene(gene_calls: pd.DataFrame) -> pd.DataFrame:
    if gene_calls.empty:
        return gene_calls
    idx = gene_calls.groupby(["species", "organelle", "gene"])["score"].idxmax()
    return gene_calls.loc[idx]


def gff_attrs(d: dict) -> str:
    return ";".join(f"{k}={v}" for k, v in d.items() if v is not None and v != "")


def emit_rna_gene(lines: list[str], contig: str, start: int, end: int, strand: str,
                   score: float, gene: str, feature: str) -> None:
    gid = f"gene-{gene}"
    lines.append("\t".join([contig, "gff_export", "gene", str(start + 1), str(end), f"{score:.3f}", strand, ".",
                             gff_attrs({"ID": gid, "Name": gene, "annotation_tier": "raw"})]))
    lines.append("\t".join([contig, "gff_export", feature, str(start + 1), str(end), f"{score:.3f}", strand, ".",
                             gff_attrs({"ID": f"{feature.lower()}-{gene}", "Parent": gid, "Name": gene})]))


def emit_protein_coding_gene(lines: list[str], gene: str, tier: str, parts: list[dict],
                              edits: list[dict], gene_score: float, n_edits_summary: int) -> None:
    """parts: list of {contig, start, end, strand, score, exon_number} in exon order.
    edits: list of {genomic_pos, resulting_aa, [slot]} to attach as sequence_alteration
    children - `slot` (trans_splicing) names the exon, and so the contig, an edit sits on.

    A trans-spliced gene's exons can sit on different contigs; its gene/mRNA
    then get one row per contig (same ID - GFF3's discontinuous-feature
    mechanism), each spanning only that contig's exons."""
    gene_strand = parts[0]["strand"] if len(parts) == 1 else "."  # mixed-strand for real trans-spliced genes
    gid, mid, cid = f"gene-{gene}", f"mRNA-{gene}", f"cds-{gene}"

    by_contig: dict[str, list[dict]] = defaultdict(list)
    for p in parts:
        by_contig[p["contig"]].append(p)
    for contig, cparts in by_contig.items():
        span = [str(min(p["start"] for p in cparts) + 1), str(max(p["end"] for p in cparts))]
        lines.append("\t".join([contig, "gff_export", "gene", *span, f"{gene_score:.3f}", gene_strand, ".",
                                 gff_attrs({"ID": gid, "Name": gene, "annotation_tier": tier})]))
        lines.append("\t".join([contig, "gff_export", "mRNA", *span, f"{gene_score:.3f}", gene_strand, ".",
                                 gff_attrs({"ID": mid, "Parent": gid, "Name": gene,
                                            "n_exons": len(parts), "n_edits": n_edits_summary})]))
    for p in parts:
        lines.append("\t".join([p["contig"], "gff_export", "CDS", str(p["start"] + 1), str(p["end"]),
                                 f"{p['score']:.3f}", p["strand"], "0",
                                 gff_attrs({"ID": cid, "Parent": mid, "exon_number": p["exon_number"]})]))
    by_exon = {p["exon_number"]: p for p in parts}
    for i, e in enumerate(edits, start=1):
        pos = int(e["genomic_pos"])
        exon = by_exon.get(int(e["slot"]), parts[0]) if pd.notna(e.get("slot")) else parts[0]
        lines.append("\t".join([exon["contig"], "gff_export", "sequence_alteration", str(pos + 1), str(pos + 1),
                                 ".", exon["strand"], ".",
                                 gff_attrs({"ID": f"sa-{gene}-{i}", "Parent": cid,
                                            "edited_from": "C", "edited_to": "T",
                                            "resulting_aa": e.get("resulting_aa", "")})]))


def build_species_gff(sub: pd.DataFrame, editing_calls: pd.DataFrame, editing_edits: pd.DataFrame,
                       ts_genes: pd.DataFrame, ts_exons: pd.DataFrame, ts_edits: pd.DataFrame) -> list[str]:
    """Every table must already be restricted to one species AND one
    organelle - genes like atpA/rpl16/rps3 exist in both genomes, and
    matching on species+gene alone let one organelle's editing call or
    reconstruction leak into the other's GFF."""
    lines: list[str] = []

    for row in sub.itertuples():
        gene = row.gene
        if gene.startswith("trn"):
            emit_rna_gene(lines, row.contig_id, row.start, row.end, row.strand, row.score, gene, "tRNA")
            continue
        if gene.startswith("rrn"):
            emit_rna_gene(lines, row.contig_id, row.start, row.end, row.strand, row.score, gene, "rRNA")
            continue

        if gene in TRANS_SPLICING_GENES:
            g = ts_genes[ts_genes.gene == gene]
            if not g.empty and int(g.iloc[0].n_filled) >= 1:
                grow = g.iloc[0]
                exons = ts_exons[ts_exons.gene == gene].sort_values("slot")
                parts = [{"contig": e.contig, "start": int(e.start), "end": int(e.end), "strand": e.strand,
                          "score": e.score, "exon_number": int(e.slot)} for e in exons.itertuples()]
                edits = ts_edits[ts_edits.gene == gene].to_dict("records")
                emit_protein_coding_gene(lines, gene, "reconstructed", parts, edits,
                                          grow.whole_gene_score, int(grow.whole_gene_n_edits))
                continue

        e = editing_calls[editing_calls.gene == gene]
        if not e.empty:
            erow = e.iloc[0]
            part = {"contig": row.contig_id, "start": int(erow.corrected_start), "end": int(erow.corrected_end),
                    "strand": row.strand, "score": erow.score, "exon_number": 1}
            edits = editing_edits[editing_edits.gene == gene].to_dict("records")
            emit_protein_coding_gene(lines, gene, "edited", [part], edits, erow.score, int(erow.n_edits))
            continue

        part = {"contig": row.contig_id, "start": int(row.start), "end": int(row.end),
                "strand": row.strand, "score": row.score, "exon_number": 1}
        emit_protein_coding_gene(lines, gene, "raw", [part], [], row.score, 0)

    def sort_key(line: str):
        f = line.split("\t")
        return (f[0], int(f[3]))
    lines.sort(key=sort_key)
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--gene-calls", default=None, help="raw gene_calls.tsv to use for tier-3 fallback + "
                     "tRNA/rRNA + contig/strand lookups (default: denovo_annotation/results/gene_calls.tsv)")
    args = ap.parse_args()

    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    raw_path = Path(args.gene_calls) if args.gene_calls else ANALYSIS_DIR / "denovo_annotation" / "results" / "gene_calls.tsv"
    if not raw_path.exists():
        print(f"[err] missing {raw_path} - run denovo_annotation/src/04_build_gene_calls.py first "
              f"(or pass --gene-calls)", file=sys.stderr)
        sys.exit(1)
    raw = pd.read_csv(raw_path, sep="\t")
    if species_filter is not None:
        raw = raw[raw.species.isin(species_filter)]
    raw_best = best_hit_per_gene(raw)

    def load(path: Path, columns: list[str]) -> pd.DataFrame:
        # Always returns a DataFrame with `columns` present (empty if the
        # upstream file doesn't exist yet, e.g. plastid before editing/
        # trans_splicing ever cover it) so every `.species`/`.gene` filter
        # downstream is safe regardless of which modules have been run.
        if not path.exists():
            return pd.DataFrame(columns=columns)
        df = pd.read_csv(path, sep="\t")
        return df[df.species.isin(species_filter)] if species_filter is not None and not df.empty else df

    editing_calls = load(ANALYSIS_DIR / "editing" / "results" / "edited_gene_calls.tsv",
                          ["species", "organelle", "gene", "corrected_start", "corrected_end", "score", "n_edits"])
    editing_edits = load(ANALYSIS_DIR / "editing" / "results" / "predicted_edits.tsv",
                          ["species", "organelle", "gene", "genomic_pos", "resulting_aa"])
    ts_genes = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_genes.tsv",
                     ["species", "organelle", "gene", "n_filled", "whole_gene_score", "whole_gene_n_edits"])
    ts_exons = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_exons.tsv",
                     ["species", "organelle", "gene", "slot", "contig", "start", "end", "strand", "score"])
    ts_edits = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_edits.tsv",
                     ["species", "organelle", "gene", "slot", "genomic_pos", "resulting_aa"])

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    out_dir = ANALYSIS_DIR / "gff_export" / "results"
    unitig_out_dir = out_dir / "unitig"
    unitig_out_dir.mkdir(parents=True, exist_ok=True)

    tables = [editing_calls, editing_edits, ts_genes, ts_exons, ts_edits]
    empty = [t.iloc[0:0] for t in tables]
    split = [{k: g for k, g in t.groupby(["species", "organelle"])} if not t.empty else {} for t in tables]
    raw_split = {k: g for k, g in raw_best.groupby(["species", "organelle"])}

    n_written = n_unitig = n_no_fasta = 0
    for organelle in organelles:
        species_here = sorted(raw_best[raw_best.organelle == organelle].species.unique())
        ctg_fastas = {r.species: r.ctg_fasta for r in
                      sd.discover_all(sd.repo_data_root(ROOT_DIR, organelle), organelle,
                                      species_filter=set(species_here))
                      if r.status == "ok"}
        for species in species_here:
            key = (species, organelle)
            lines = build_species_gff(raw_split[key], *[s.get(key, e) for s, e in zip(split, empty)])
            if not lines:
                continue
            out_path = out_dir / f"{species}.{organelle}.gff"
            if species not in ctg_fastas:
                n_no_fasta += 1
                out_path.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
                n_written += 1
                continue

            ctx = load_seq_context(Path(ctg_fastas[species]), UNITIG_MAP_DIR / f"{species}.{organelle}.tsv")
            ctg_lines, unitig_lines = add_unitig_layer(lines, ctx)
            out_path.write_text(assemble(*ctg_header(ctx, list(ctx.lengths)), ctg_lines))
            n_written += 1
            if unitig_lines:
                unitigs = sorted({l.split("\t", 1)[0] for l in unitig_lines}, key=lambda u: (len(u), u))
                (unitig_out_dir / f"{species}.{organelle}.unitig.gff").write_text(
                    assemble(*unitig_header(ctx, unitigs), unitig_lines))
                n_unitig += 1
    print(f"[info] build_gff: {n_written} species -> {out_dir} ({n_unitig} with unitig-level GFF -> "
          f"{unitig_out_dir}; {n_no_fasta} without a resolvable ctg.fasta, no region/unitig layer)",
          file=sys.stderr)


if __name__ == "__main__":
    main()

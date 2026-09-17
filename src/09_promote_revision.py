#!/usr/bin/env python3
"""Compare a completed rev/out/ revision rerun against the currently-promoted
data/mito or data/plastid files for the same species, and promote it ONLY
if it's objectively better (fewer dead-end nodes, fewer disconnected
subgraphs) - never on file size/recency alone.

Why this exists (a real bug this script fixes, found the hard way): a
species_discovery.py tie-break (prefer the larger ctg.fasta when both
candidates are otherwise complete) was designed to reject a truly-empty
failed rerun, but silently picks the WORSE of two genuinely valid runs
when the better one happens to be smaller - which is common, since fixing
dead-ends/fragmentation by lowering coverage often produces a cleaner,
smaller assembly. Promoting blindly on file presence, or leaving both
candidates to compete under that heuristic, is not safe. This script
instead: (1) explicitly compares the two using the same metrics
qc_basic_stats uses, (2) only copies in a rerun that is verifiably better,
and (3) moves the superseded old files into a `superseded/` subdirectory
within the species folder - NOT deleted (nothing here is ever destroyed),
and NOT left in the same directory to compete (species_discovery's glob is
non-recursive, so files here are correctly invisible to it).

Verified against the first real batch: of 123 already-completed-but-never-
promoted historical reruns sitting in rev/out/, only 59/109 comparable
mito results and 7/120 plastid results were actually better - the rest
(38 mito, 6 plastid) would have made things WORSE if promoted blindly.

Usage:
  src/09_promote_revision.py --list rev/lists/pilot_v2_deadend_c40.txt --organelle mito
  src/09_promote_revision.py --list rev/lists/pilot_v2_deadend_c40.txt --organelle both
  src/09_promote_revision.py --all-revout --organelle both   # scan every species under rev/out/

Output: prints a promoted/worse/same/incomparable summary; writes the full
per-species comparison to meta/promotion_log_<timestamp>.tsv (append-only
record of what was compared and decided, when).
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parent.parent
GFATK = shutil.which("gfatk") or str(Path.home() / ".cargo" / "bin" / "gfatk")

SUFFIXES = {
    "mito": [".mito.gfa", ".mito.ctg.bed", ".mito.ctg.fasta", ".mito.bed", ".annot_mito.txt"],
    "pltd": [".pltd.gfa", ".pltd.ctg.bed", ".pltd.ctg.fasta", ".pltd.bed", ".annot_pltd.txt"],
}
DATA_DIRNAME = {"mito": "mito", "pltd": "plastid"}


def gfatk_metrics(gfa_path: Path) -> dict | None:
    proc = subprocess.run([GFATK, "stats", str(gfa_path)], capture_output=True, text=True, timeout=60)
    if proc.returncode != 0:
        return None
    text = proc.stdout
    m = re.search(r"Total number of subgraphs:\s*(\d+)", text)
    n_sub = int(m.group(1)) if m else None
    if n_sub == 0:
        return {"n_subgraphs": 0, "dead_end_nodes": None}
    dead_ends = [int(x) for x in re.findall(r"Dead-end nodes:\s*(\d+)", text)]
    return {"n_subgraphs": n_sub, "dead_end_nodes": max(dead_ends) if dead_ends else None}


SECONDARY_REGRESSION_THRESHOLD = 2  # absolute increase beyond which a secondary
                                     # metric counts as a *meaningful* regression,
                                     # not just tie-break noise (+1 subgraph is
                                     # common/trivial; +24, as verified on real
                                     # data for one species, clearly is not).


def verdict(old: dict, new: dict, primary_metric: str | None = None) -> str:
    """primary_metric: None for the original balanced comparison (both
    dead_end_nodes and n_subgraphs must not worsen) - used for categories with
    no single clear target metric (core_gene_pct, high_cross_organelle_alignment).
    Otherwise "dead_end_nodes" or "n_subgraphs": verified on real data that
    requiring BOTH metrics to improve is too strict - a rerun aimed at fixing
    dead-ends can legitimately improve dead-ends while only slightly nudging
    subgraph count, and blanket-labelling that "worse" mislabels a real fix
    (e.g. one pilot species: dead-ends 6->0, subgraphs 3->4, wrongly called
    "worse" under the old logic). Returns a fourth outcome, "mixed", for the
    case where the PRIMARY metric improved but the secondary one regressed by
    more than SECONDARY_REGRESSION_THRESHOLD - flagged for manual review, but
    NOT promoted automatically (verified on real data this genuinely happens:
    one pilot species improved dead-ends 7->4 while its subgraph count
    exploded 1->25, which is a real problem, not noise).
    """
    # Resolved-contig-fasta presence is checked FIRST and independently of
    # everything else: verified on real data that a rerun can fix exactly
    # this (Pathfinder now emits a .ctg.fasta where it previously couldn't)
    # while dead_end_nodes/n_subgraphs on the raw graph look identical - a
    # graph-metrics-only comparison would otherwise silently miss this whole
    # class of real improvement (9/10 in the first no_resolved_ctg_fasta pilot).
    old_has_fasta = old.get("has_ctg_fasta", True)  # default True: most callers don't check this dimension
    new_has_fasta = new.get("has_ctg_fasta", True)
    if new_has_fasta and not old_has_fasta:
        return "better"
    if old_has_fasta and not new_has_fasta:
        return "worse"

    od = 999 if old.get("dead_end_nodes") is None else old["dead_end_nodes"]
    os_ = 999 if old.get("n_subgraphs") is None else old["n_subgraphs"]
    nd = 999 if new.get("dead_end_nodes") is None else new["dead_end_nodes"]
    ns = 999 if new.get("n_subgraphs") is None else new["n_subgraphs"]

    if primary_metric in ("dead_end_nodes", "n_subgraphs"):
        primary_old, primary_new = (od, nd) if primary_metric == "dead_end_nodes" else (os_, ns)
        secondary_old, secondary_new = (os_, ns) if primary_metric == "dead_end_nodes" else (od, nd)
        if primary_new < primary_old:
            if secondary_new > secondary_old + SECONDARY_REGRESSION_THRESHOLD:
                return "mixed"
            return "better"
        if primary_new > primary_old:
            return "worse"
        # primary unchanged - fall through to the secondary metric alone
        if secondary_new < secondary_old:
            return "better"
        if secondary_new > secondary_old:
            return "worse"
        return "same"

    # No primary metric for this category (core_gene_pct, high_cross_organelle_alignment):
    # original balanced comparison, both metrics must not worsen.
    if nd <= od and ns <= os_ and (nd < od or ns < os_):
        return "better"
    if nd > od or ns > os_:
        return "worse"
    return "same"


def current_metrics(data_root: Path, species: str, organelle: str) -> tuple[dict, str | None]:
    """Best-effort: ask gfatk directly about whatever gfa species_discovery
    would currently pick (mirrors its own candidate selection, minimally)."""
    sys.path.insert(0, str(ROOT_DIR / "analysis" / "common"))
    import species_discovery as sd  # noqa: E402
    r = sd.resolve_species(data_root / species, organelle)
    if not r.gfa:
        return {"n_subgraphs": None, "dead_end_nodes": None, "has_ctg_fasta": False}, None
    m = gfatk_metrics(Path(r.gfa)) or {"n_subgraphs": None, "dead_end_nodes": None}
    m["has_ctg_fasta"] = r.status == "ok"
    return m, r.gfa


def promote(species: str, organelle: str, new_gfa: Path) -> None:
    data_root = ROOT_DIR / "data" / DATA_DIRNAME[organelle]
    species_dir = data_root / species
    species_dir.mkdir(parents=True, exist_ok=True)
    gfa_suffix = SUFFIXES[organelle][0]
    new_prefix = new_gfa.name[: -len(gfa_suffix)]
    src_dir = new_gfa.parent

    # archive whatever's currently there (exact suffix match only - never a loose
    # prefix glob, which can wrongly sweep up an unrelated longer-prefixed file)
    import species_discovery as sd  # already on sys.path from current_metrics()
    old_prefixes = [p for p in sd._candidate_prefixes(species_dir, gfa_suffix) if p != new_prefix]
    if old_prefixes:
        archive_dir = species_dir / "superseded"
        archive_dir.mkdir(exist_ok=True)
        for old_prefix in old_prefixes:
            for suf in SUFFIXES[organelle]:
                f = species_dir / f"{old_prefix}{suf}"
                if f.exists():
                    shutil.move(str(f), str(archive_dir / f.name))

    for suf in SUFFIXES[organelle]:
        f = src_dir / f"{new_prefix}{suf}"
        if f.exists() and f.stat().st_size > 0:
            shutil.copy2(f, species_dir / f.name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", default=None, help="species list file (one per line)")
    ap.add_argument("--all-revout", action="store_true", help="scan every species directory under rev/out/")
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--primary-metric", default=None, choices=[None, "dead_end_nodes", "n_subgraphs"],
                     help="the metric this batch's recipe specifically targets (e.g. dead_end_nodes for a "
                          "deadend_c40 rerun, n_subgraphs for a fragmented_c40 one) - improving this is enough "
                          "for 'better' even if the other metric ticks up slightly; omit for categories with no "
                          "single clear target (core_gene_pct, high_cross_organelle_alignment), which keeps the "
                          "stricter balanced comparison (both metrics must not worsen)")
    ap.add_argument("--dry-run", action="store_true", help="compare and report only, don't copy/archive anything")
    args = ap.parse_args()

    if args.all_revout:
        species_list = sorted(p.name for p in (ROOT_DIR / "rev" / "out").iterdir() if p.is_dir())
    elif args.list:
        species_list = [l.strip() for l in open(args.list) if l.strip() and not l.startswith("#")]
    else:
        print("[err] pass --list FILE or --all-revout", file=sys.stderr)
        sys.exit(1)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    rows = []
    for organelle in organelles:
        data_root = ROOT_DIR / "data" / DATA_DIRNAME[organelle]
        gfa_suffix = SUFFIXES[organelle][0]
        for sp in species_list:
            revout_dir = ROOT_DIR / "rev" / "out" / sp
            cands = [g for g in revout_dir.glob(f"*{gfa_suffix}") if g.stat().st_size > 0] if revout_dir.is_dir() else []
            if not cands:
                continue
            new_gfa = max(cands, key=lambda g: g.stat().st_size)
            new_m = gfatk_metrics(new_gfa)
            if new_m is None:
                continue
            fasta_suffix = SUFFIXES[organelle][2]  # ".mito.ctg.fasta" / ".pltd.ctg.fasta"
            new_prefix = new_gfa.name[: -len(gfa_suffix)]
            new_fasta = new_gfa.parent / f"{new_prefix}{fasta_suffix}"
            new_m["has_ctg_fasta"] = new_fasta.exists() and new_fasta.stat().st_size > 0
            old_m, old_gfa = current_metrics(data_root, sp, organelle)
            v = verdict(old_m, new_m, primary_metric=args.primary_metric)
            rows.append({
                "species": sp, "organelle": organelle, "verdict": v,
                "old_gfa": old_gfa, "old_dead_end_nodes": old_m.get("dead_end_nodes"),
                "old_n_subgraphs": old_m.get("n_subgraphs"),
                "new_gfa": str(new_gfa), "new_dead_end_nodes": new_m.get("dead_end_nodes"),
                "new_n_subgraphs": new_m.get("n_subgraphs"),
            })
            if v == "better" and not args.dry_run:
                promote(sp, organelle, new_gfa)

    df = pd.DataFrame(rows)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = ROOT_DIR / "meta" / f"promotion_log_{ts}.tsv"
    df.to_csv(log_path, sep="\t", index=False)
    print(df.groupby(["organelle", "verdict"]).size() if not df.empty else "no comparable species found",
          file=sys.stderr)
    print(f"[info] {'(dry run) ' if args.dry_run else ''}log -> {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

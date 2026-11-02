#!/usr/bin/env python3
"""Reshape oatk's own .ctg.bed gene calls into one tidy long table.

No new annotation is performed here - oatk already annotates every contig
(that's what .ctg.bed / .annot_*.txt already are); this just adds
species/organelle provenance columns and applies the QC filter so
downstream modules read one trusted table instead of re-deriving this
themselves.

Output: analysis/annotation/results/gene_calls.tsv
Columns: species organelle contig_id start end gene score strand
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = ANALYSIS_DIR.parent
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import species_discovery as sd  # noqa: E402
import io_utils  # noqa: E402

BED_COLS = ["contig_id", "start", "end", "gene", "score", "strand"]
COLUMNS = ["species", "organelle"] + BED_COLS


def qc_eligible_species(qc_status: list[str]) -> set[tuple[str, str]]:
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    if not qc_path.exists():
        print(f"[err] missing {qc_path} - run qc_basic_stats/src/05_qc_summary.py first", file=sys.stderr)
        sys.exit(1)
    qc = pd.read_csv(qc_path, sep="\t")
    sub = qc[qc["status"].isin(qc_status)]
    return set(zip(sub["species"], sub["organelle"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail")
    ap.add_argument("--force", action="store_true", help="accepted for CLI consistency with run_all.sh; "
                     "this script always does a full recompute (cheap), so it's a no-op")
    args = ap.parse_args()

    qc_status = args.qc_status.split(",")
    eligible = qc_eligible_species(qc_status)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    all_rows = []
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        n_used = 0
        for r in resolved:
            if (r.species, organelle) not in eligible or not r.ctg_bed:
                continue
            df = pd.read_csv(r.ctg_bed, sep="\t", comment="#", names=BED_COLS)
            if df.empty:
                continue
            df["species"] = r.species
            df["organelle"] = organelle
            all_rows.append(df[COLUMNS])
            n_used += 1
        print(f"[info] gene_table: {organelle}: {n_used} QC-eligible species with gene calls", file=sys.stderr)

    out_path = ANALYSIS_DIR / "annotation" / "results" / "gene_calls.tsv"
    # Always write what this run actually computed to its own path first;
    # merge_species_scoped_tsv then decides how that relates to what's
    # already in out_path (replace entirely if unscoped, species-scoped
    # merge otherwise - see its docstring for why that distinction matters:
    # a --species-list-scoped run must never overwrite other species' rows).
    out_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(columns=COLUMNS)
    new_path = ANALYSIS_DIR / "annotation" / "work" / "gene_calls.new.tsv"
    new_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(new_path, sep="\t", index=False)
    io_utils.merge_species_scoped_tsv(out_path, new_path, species_filter)
    print(f"[info] gene_table: {len(out_df)} rows "
          f"({'species-scoped merge' if species_filter else 'full rewrite'}) -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

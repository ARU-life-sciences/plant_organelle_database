#!/usr/bin/env python3
"""Revision plan v2: reads analysis/qc_basic_stats/results/qc_summary.tsv
(the rich, both-organelle QC gate) rather than the old mito-only
meta/dead_ends.tsv that 05_revision_plan.sh used, and assigns each failing
species a specific oatk re-run recipe (coverage / max-gb) with a written
justification - not just a bare parameter change.

Scope: species with a status of 'fail' on at least one organelle, PLUS
species whose only issue is no_resolved_ctg_fasta (a flag-tier-only
reason under the current QC thresholds - it can never trigger 'fail' - but
judged to need the same reassembly treatment: the graph exists but oatk's
Pathfinder couldn't resolve a path through it, which won't self-resolve by
re-running Pathfinder alone on the same graph). Only reasons judged
fixable by re-running oatk with different parameters are in scope (see
PARAMETER STRATEGY below) - the other flag-only reasons (non_circular,
size_outlier, elevated_cross_organelle_alignment 15-35%) are NOT included:
those are mostly real biological variation, not assembly defects, and
reassembling on that basis alone would be guessing.

A species can fail for more than one reason (46/335 do). Rather than try
to combine parameter changes, each species gets ONE recipe, chosen by
PRIORITY (most clearly a graph-connectivity problem first, most ambiguous
last) - see PRIORITY below.

============================== PARAMETER STRATEGY ==============================

dead_end_nodes (64 species as sole reason)
  -> coverage=40, max_gb=15 (suffix "deadend_c40")
  Already-established project convention (see top-level README's
  "Revising the outputs" section, written before this script existed):
  lowering the minimum k-mer coverage threshold admits rarer k-mers into
  the graph, which can complete/bridge a path that a stricter threshold
  left as a dead end. Costs some graph noise/complexity in exchange.

fragmented (n_subgraphs>=10; 34 species as sole reason)
  -> coverage=40, max_gb=15 (suffix "fragmented_c40")
  Same underlying mechanism as dead-ends: the graph failed to connect
  regions that should be joined. Treated with the same fix on the same
  rationale. This is the least-validated mapping of the six - fragmentation
  could in principle also be caused by too much low-frequency NOISE
  (in which case *raising* coverage would help instead) rather than too
  little signal. Flagged here explicitly so the pilot's outcome for this
  category gets checked before scaling to the full 34+ species - if lower
  coverage doesn't measurably reduce n_subgraphs in the pilot, the
  opposite adjustment should be tried for the remainder, not assumed.

no_resolved_ctg_fasta (52 species as sole reason)
  -> coverage=40, max_gb=15 (suffix "unresolved_c40")
  The assembly graph exists but oatk's Pathfinder step couldn't resolve a
  clean/circular path through it well enough to emit a .ctg.fasta.
  Re-running Pathfinder alone (oatk -G <existing.gfa>) on the SAME graph
  with the SAME thresholds would almost certainly reproduce the same
  failure - it's deterministic. So despite this being the cheapest-looking
  fix on paper, it needs the same graph-quality intervention as dead-ends/
  fragmentation (re-assemble with lower coverage to get a cleaner graph
  Pathfinder can actually resolve), not just a re-run of the last step.

missing_gfa (34 species as sole reason)
  -> coverage=80, max_gb=16 (suffix "missing_c80", i.e. the ORIGINAL defaults)
  The assembly never completed at all. With no prior attempt to diagnose,
  the first re-run should just be a normal default-parameter attempt, not
  a speculative adjustment. (If a memory-kill is later found in this
  species' logs, use the project's existing established mem-fix recipe
  instead: coverage=120, max_gb=15 - see top-level README - not this one.)

core_gene_pct (<50%; 12 species as sole reason)
  -> coverage=80 (default), max_gb=25 (suffix "incomplete_g25")
  Low core-gene recovery suggests the assembly is missing large stretches
  of real sequence - most directly addressed by giving the assembler MORE
  raw input data (raising the max-gb budget from the original run's cap),
  not by changing the coverage threshold, which governs graph resolution/
  cleanliness rather than how much data it has to work with in the first
  place.

high_cross_organelle_alignment (>35%; 93 species as sole reason)
  -> coverage=80, max_gb=16 (suffix "crossorg_c80", i.e. defaults - a
     baseline re-check, not a targeted fix)
  Unlike the other five, this isn't obviously a graph-connectivity/
  completeness problem, so no coverage/data-volume adjustment is
  obviously "the fix". A single default-parameter re-assembly serves to
  rule out "was the original assembly just unusually noisy/mis-binned"
  (visually confirmed for one example, Juncus_inflexus, to be a
  completely clean plastid graph despite this flag - i.e. sometimes this
  reflects a genuinely large MTPT event, not an assembly defect at all).
  Species still flagged after this one attempt should be documented as
  likely real rather than repeatedly reassembled with different
  parameters on no particular basis.

PRIORITY (for the 46 species with more than one reason - pick ONE recipe,
most-clearly-mechanical first, most-ambiguous last, rather than combining
parameter changes):
  dead_end_nodes > fragmented > no_resolved_ctg_fasta > missing_gfa
  > core_gene_pct > high_cross_organelle_alignment

Output: meta/revision_plan_v2.tsv
Columns: species organelle_reasons primary_reason all_reasons action coverage max_gb suffix rationale
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parent.parent

# priority order: earliest = most clearly a graph-mechanical problem
PRIORITY = ["dead_end_nodes", "fragmented", "no_resolved_ctg_fasta",
            "missing_gfa", "core_gene_pct", "high_cross_organelle_alignment"]

RECIPES = {
    "dead_end_nodes": dict(coverage=40, max_gb=15, suffix="deadend_c40",
        rationale="lower coverage admits rarer k-mers, can bridge a dead-end path (established project convention)"),
    "fragmented": dict(coverage=40, max_gb=15, suffix="fragmented_c40",
        rationale="same graph-connectivity mechanism as dead-ends; LEAST validated mapping, check pilot outcome before scaling"),
    "no_resolved_ctg_fasta": dict(coverage=40, max_gb=15, suffix="unresolved_c40",
        rationale="Pathfinder failure on existing graph is deterministic - needs a cleaner re-assembled graph, not just a Pathfinder re-run"),
    "missing_gfa": dict(coverage=80, max_gb=16, suffix="missing_c80",
        rationale="assembly never completed; first attempt uses original defaults, no speculative adjustment"),
    "core_gene_pct": dict(coverage=80, max_gb=25, suffix="incomplete_g25",
        rationale="low completeness suggests too little input data, not a graph-resolution problem; raise max-gb"),
    "high_cross_organelle_alignment": dict(coverage=80, max_gb=16, suffix="crossorg_c80",
        rationale="not obviously a connectivity/data-volume problem; baseline re-check only, may well be real biology (e.g. Juncus_inflexus)"),
}


def main():
    qc_path = ROOT_DIR / "analysis" / "qc_basic_stats" / "results" / "qc_summary.tsv"
    df = pd.read_csv(qc_path, sep="\t")

    def matched(reasons):
        if pd.isna(reasons):
            return []
        return [r for r in PRIORITY if r in reasons]

    # In scope: every 'fail' row, PLUS 'flag' rows whose only relevant issue is
    # no_resolved_ctg_fasta (a flag-tier-only reason - it can never trigger
    # 'fail' under the current QC thresholds - but per-species review judged
    # these need the same reassembly treatment as the fail-tier reasons, not
    # the "leave alone" treatment given to the other flag-only reasons like
    # non_circular/size_outlier/elevated_cross_organelle_alignment(15-35%)).
    in_scope = df[(df["status"] == "fail") |
                  ((df["status"] == "flag") & df["status_reasons"].str.contains("no_resolved_ctg_fasta", na=False))]

    per_species: dict[str, dict] = {}
    for _, row in in_scope.iterrows():
        hits = matched(row["status_reasons"])
        if not hits:
            continue  # a fail for a reason not in scope here (shouldn't happen given current checks)
        entry = per_species.setdefault(row["species"], {"organelle_reasons": [], "all_reasons": set()})
        entry["organelle_reasons"].append(f"{row['organelle']}:{row['status_reasons']}")
        entry["all_reasons"].update(hits)

    rows = []
    for species, entry in sorted(per_species.items()):
        all_reasons = entry["all_reasons"]
        primary = next(r for r in PRIORITY if r in all_reasons)
        recipe = RECIPES[primary]
        rows.append({
            "species": species,
            "organelle_reasons": ";".join(entry["organelle_reasons"]),
            "primary_reason": primary,
            "all_reasons": ",".join(sorted(all_reasons)),
            "action": f"resubmit_{recipe['suffix']}",
            "coverage": recipe["coverage"],
            "max_gb": recipe["max_gb"],
            "suffix": recipe["suffix"],
            "rationale": recipe["rationale"],
        })

    out_path = ROOT_DIR / "meta" / "revision_plan_v2.tsv"
    out_df = pd.DataFrame(rows, columns=["species", "organelle_reasons", "primary_reason", "all_reasons",
                                          "action", "coverage", "max_gb", "suffix", "rationale"])
    out_df.to_csv(out_path, sep="\t", index=False)

    print(f"[info] revision_plan_v2: {len(out_df)} species -> {out_path}", file=sys.stderr)
    print(out_df["primary_reason"].value_counts(), file=sys.stderr)


if __name__ == "__main__":
    main()

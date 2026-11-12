"""How much of an assembly graph a resolved contig set accounts for.

Used to gate gfatk-resolve promotions (`linearize/src/05_promote_resolve.py`,
`07_restore_pathfinder.py`): a circuit that walks only part of the graph
(Cyperaceae plastids came out 34-55 kb of a ~150-230 kb genome) must not
replace anything. Placement uses the same sequence-matching logic
`unitig_coords/src/00_build_unitig_map.py` used before that module was
retired in favor of the standalone `poa` package (793c2a17) - inlined here
(gfa/fasta parsing, revcomp, tolerant doubled-haystack substring search)
since this gate still needs it and poa doesn't expose an equivalent. A
number here means the same thing `pct_unitig_content_placed` used to mean
in `unitig_coords/results/linearization_qc.tsv`.
"""
from __future__ import annotations

from pathlib import Path

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_MIN_SEQMATCH_BP = 50
MAX_TRIM_BP = 2000


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def read_gfa(gfa_path: Path) -> dict[str, str]:
    """Returns segment_id -> sequence."""
    segs: dict[str, str] = {}
    with open(gfa_path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                segs[parts[1]] = parts[2]
    return segs


def read_ctg_fasta(fasta_path: Path) -> dict[str, str]:
    """Returns seqid -> sequence."""
    seqs: dict[str, list[str]] = {}
    name = None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def _hits(haystack_doubled: str, hay_len: int, sub: str) -> list[int]:
    positions = []
    start = 0
    while True:
        p = haystack_doubled.find(sub, start)
        if p == -1 or p >= hay_len:
            return positions
        positions.append(p)
        start = p + 1


def find_all_with_tolerance(haystack_doubled: str, hay_len: int, needle: str,
                             max_trim: int) -> tuple[list[int], int]:
    """Search a self-doubled haystack for `needle`, allowing the smallest
    symmetric trim (0..max_trim) off both ends that yields a hit (handles a
    circular join's overlap trim). Hits are monotone in trim - a trimmed
    needle is a substring of every less-trimmed one - so binary search
    finds the smallest trim in ~log2(max_trim) scans. Returns (start
    positions of the *untrimmed* needle modulo hay_len, trim used) - only
    the needle's middle [trim, len - trim) is verified at those positions."""
    max_trim = min(max_trim, (len(needle) - _MIN_SEQMATCH_BP) // 2)
    if max_trim < 0:
        return [], 0
    exact = _hits(haystack_doubled, hay_len, needle)
    if exact:
        return exact, 0
    if max_trim == 0 or not _hits(haystack_doubled, hay_len, needle[max_trim:len(needle) - max_trim]):
        return [], 0
    lo, hi = 0, max_trim  # lo misses, hi hits
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if _hits(haystack_doubled, hay_len, needle[mid:len(needle) - mid]):
            hi = mid
        else:
            lo = mid
    return [(p - hi) % hay_len for p in _hits(haystack_doubled, hay_len, needle[hi:len(needle) - hi])], hi


def unitig_content_placed(gfa: Path, ctg_fasta: Path) -> float:
    """Fraction (0-1) of the graph's unique unitig bp found in `ctg_fasta`
    (forward or reverse complement; circular joins and trimmed repeat
    variants handled by the tolerant search above)."""
    segs = read_gfa(gfa)
    ctgs = read_ctg_fasta(ctg_fasta)
    total = sum(len(s) for s in segs.values())
    if not total or not ctgs:
        return 0.0
    doubled = {cid: c + c for cid, c in ctgs.items()}
    placed = 0
    for useq in segs.values():
        max_trim = min(MAX_TRIM_BP, len(useq) // 4)
        rc = revcomp(useq)
        if any(find_all_with_tolerance(doubled[cid], len(c), needle, max_trim)[0]
               for cid, c in ctgs.items() for needle in (useq, rc)):
            placed += len(useq)
    return placed / total

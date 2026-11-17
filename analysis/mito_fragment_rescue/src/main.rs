// mito_fragment_rescue
//
// For species where oatk's own HMM annotation shows real mitochondrial gene
// hits sitting on a unitig that its winner-take-all classifier assigned
// entirely to the plastid (so mito.gfa comes out empty), pull out just the
// gene-level fragments that are NOT attributable to plastid/mito HMM
// cross-reactivity, and write them out as a labelled, fragmentary rescue set.
//
// This does not produce a mito.gfa or a genome assembly. It produces
// sequence fragments for specific genes, cross-checked per-run against that
// same run's own plastid annotation (same unitig graph, so coordinates are
// directly comparable - this is why the check is done per-candidate-run
// rather than against a separately-run "official" plastid annotation).
//
// Usage:
//   mito_fragment_rescue --root <repo_root> [--species-list FILE | SPECIES...]
//                         [--evalue E] [--slop N] [--out-dir DIR]

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

const DEFAULT_EVALUE_THRESH: f64 = 1e-6;
const DEFAULT_SLOP: i64 = 50;

fn mito_exclusive_genes() -> HashSet<&'static str> {
    [
        "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6", "nad7", "nad9", "cox1", "cox2",
        "cox3", "cob", "mttB", "matR", "ccmB", "ccmC", "ccmFc", "ccmFn", "sdh3", "sdh4",
    ]
    .into_iter()
    .collect()
}

#[derive(Clone, Debug)]
struct Hit {
    gene: String,
    query: String,
    lo: i64,
    hi: i64,
    evalue: f64,
    score: f64,
}

fn parse_annot(path: &Path) -> Vec<Hit> {
    let mut out = Vec::new();
    let text = match fs::read_to_string(path) {
        Ok(t) => t,
        Err(_) => return out,
    };
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 15 {
            continue;
        }
        let gene = parts[0].to_string();
        let query = parts[2].to_string();
        let (lo, hi) = match (parts[6].parse::<i64>(), parts[7].parse::<i64>()) {
            (Ok(a), Ok(b)) => (a, b),
            _ => continue,
        };
        let (evalue, score) = match (parts[12].parse::<f64>(), parts[13].parse::<f64>()) {
            (Ok(a), Ok(b)) => (a, b),
            _ => continue,
        };
        out.push(Hit {
            gene,
            query,
            lo,
            hi,
            evalue,
            score,
        });
    }
    out
}

fn overlaps(a_lo: i64, a_hi: i64, b_lo: i64, b_hi: i64, slop: i64) -> bool {
    let (lo1, hi1) = (a_lo.min(a_hi), a_lo.max(a_hi));
    let (lo2, hi2) = (b_lo.min(b_hi), b_lo.max(b_hi));
    lo1 - slop <= hi2 && lo2 - slop <= hi1
}

fn revcomp(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            other => other,
        } as char)
        .collect()
}

fn parse_gfa_sequences(path: &Path) -> HashMap<String, String> {
    let mut seqs = HashMap::new();
    let text = match fs::read_to_string(path) {
        Ok(t) => t,
        Err(_) => return seqs,
    };
    for line in text.lines() {
        if !line.starts_with("S\t") {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        seqs.insert(parts[1].to_string(), parts[2].to_uppercase());
    }
    seqs
}

struct Candidate {
    prefix: String,
    graph_path: PathBuf,
    clean_hits: Vec<Hit>,
    seqs: HashMap<String, String>,
}

fn find_candidates(species_dir: &Path, evalue_thresh: f64, slop: i64) -> Vec<Candidate> {
    let mito_set = mito_exclusive_genes();
    let mut candidates = Vec::new();

    let entries = match fs::read_dir(species_dir) {
        Ok(e) => e,
        Err(_) => return candidates,
    };

    let mut annot_files: Vec<PathBuf> = Vec::new();
    for entry in entries.flatten() {
        let path = entry.path();
        let name = match path.file_name().and_then(|n| n.to_str()) {
            Some(n) => n,
            None => continue,
        };
        if name.ends_with(".annot_mito.txt") {
            if let Ok(meta) = entry.metadata() {
                if meta.len() > 0 {
                    annot_files.push(path);
                }
            }
        }
    }

    for mito_path in annot_files {
        let name = mito_path.file_name().unwrap().to_str().unwrap();
        let prefix = name.trim_end_matches(".annot_mito.txt").to_string();

        let pltd_path = species_dir.join(format!("{}.annot_pltd.txt", prefix));
        if !pltd_path.exists() {
            continue; // can't rule out clashes without the same-run plastid calls
        }

        let final_gfa = species_dir.join(format!("{}.utg.final.gfa", prefix));
        let raw_gfa = species_dir.join(format!("{}.utg.gfa", prefix));
        let graph_path = if final_gfa.exists() {
            final_gfa
        } else if raw_gfa.exists() {
            raw_gfa
        } else {
            continue; // no sequence source for this run
        };

        let mito_hits = parse_annot(&mito_path);
        let pltd_hits = parse_annot(&pltd_path);

        let mut pltd_by_query: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for h in pltd_hits.iter().filter(|h| h.evalue <= evalue_thresh) {
            pltd_by_query
                .entry(h.query.clone())
                .or_default()
                .push((h.lo, h.hi));
        }

        let mut clean_hits = Vec::new();
        for h in mito_hits
            .into_iter()
            .filter(|h| mito_set.contains(h.gene.as_str()) && h.evalue <= evalue_thresh)
        {
            let clash = pltd_by_query
                .get(&h.query)
                .map(|spans| spans.iter().any(|&(lo, hi)| overlaps(h.lo, h.hi, lo, hi, slop)))
                .unwrap_or(false);
            if !clash {
                clean_hits.push(h);
            }
        }

        // Guard against a graph file that no longer corresponds to the same
        // unitig numbering/sequences the annotation was computed against
        // (observed in practice: a later rerun overwrote utg.gfa under the
        // same filename after annot_mito.txt had already been written).
        // If any clean hit's coordinates don't fit inside its unitig's
        // sequence in this graph file, the pairing is untrustworthy and the
        // whole candidate is discarded rather than risk emitting wrong
        // sequence.
        let seqs = parse_gfa_sequences(&graph_path);
        let consistent = clean_hits.iter().all(|h| {
            let (lo, hi) = (h.lo.min(h.hi), h.lo.max(h.hi));
            lo >= 1 && seqs.get(&h.query).map(|s| hi as usize <= s.len()).unwrap_or(false)
        });
        if !consistent {
            eprintln!(
                "WARN: discarding candidate {} for {}: annotation coordinates do not fit the sequence lengths in {} (stale/mismatched graph file)",
                prefix, species_dir.display(), graph_path.display()
            );
            continue;
        }

        candidates.push(Candidate {
            prefix,
            graph_path,
            clean_hits,
            seqs,
        });
    }

    candidates
}

fn choose_best(candidates: &[Candidate]) -> Option<usize> {
    let mut best_idx: Option<usize> = None;
    let mut best_key: (usize, bool, usize) = (0, false, 0); // (n_distinct_genes, is_final_gfa, n_hits)
    for (i, c) in candidates.iter().enumerate() {
        let distinct: HashSet<&str> = c.clean_hits.iter().map(|h| h.gene.as_str()).collect();
        let is_final = c
            .graph_path
            .to_str()
            .map(|s| s.ends_with(".utg.final.gfa"))
            .unwrap_or(false);
        let key = (distinct.len(), is_final, c.clean_hits.len());
        if best_idx.is_none() || key > best_key {
            best_idx = Some(i);
            best_key = key;
        }
    }
    best_idx
}

struct Fragment {
    gene: String,
    unitig: String,
    start: i64,
    end: i64,
    strand: char,
    n_hsps: usize,
    best_evalue: f64,
    best_score: f64,
}

fn merge_fragments(hits: &[Hit], slop: i64) -> Vec<Fragment> {
    let mut by_gene_unitig: HashMap<(String, String), Vec<&Hit>> = HashMap::new();
    for h in hits {
        by_gene_unitig
            .entry((h.gene.clone(), h.query.clone()))
            .or_default()
            .push(h);
    }

    let mut fragments = Vec::new();
    for ((gene, unitig), mut group) in by_gene_unitig {
        group.sort_by_key(|h| h.lo.min(h.hi));
        let mut cur_lo = group[0].lo.min(group[0].hi);
        let mut cur_hi = group[0].lo.max(group[0].hi);
        let mut cur_members: Vec<&Hit> = vec![group[0]];

        let flush = |lo: i64, hi: i64, members: &[&Hit], frags: &mut Vec<Fragment>| {
            let best = members
                .iter()
                .min_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap())
                .unwrap();
            let strand = if best.lo <= best.hi { '+' } else { '-' };
            frags.push(Fragment {
                gene: gene.clone(),
                unitig: unitig.clone(),
                start: lo,
                end: hi,
                strand,
                n_hsps: members.len(),
                best_evalue: best.evalue,
                best_score: best.score,
            });
        };

        for h in group.into_iter().skip(1) {
            let (lo, hi) = (h.lo.min(h.hi), h.lo.max(h.hi));
            if lo - slop <= cur_hi {
                cur_hi = cur_hi.max(hi);
                cur_members.push(h);
            } else {
                flush(cur_lo, cur_hi, &cur_members, &mut fragments);
                cur_lo = lo;
                cur_hi = hi;
                cur_members = vec![h];
            }
        }
        flush(cur_lo, cur_hi, &cur_members, &mut fragments);
    }
    fragments
}

fn write_fasta_record(out: &mut String, header: &str, seq: &str) {
    out.push('>');
    out.push_str(header);
    out.push('\n');
    for chunk in seq.as_bytes().chunks(70) {
        out.push_str(std::str::from_utf8(chunk).unwrap());
        out.push('\n');
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut root = env::current_dir().unwrap();
    let mut species: Vec<String> = Vec::new();
    let mut evalue_thresh = DEFAULT_EVALUE_THRESH;
    let mut slop = DEFAULT_SLOP;
    let mut out_dir = PathBuf::from("analysis/mito_fragment_rescue/results");

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--root" => {
                i += 1;
                root = PathBuf::from(&args[i]);
            }
            "--species-list" => {
                i += 1;
                let text = fs::read_to_string(&args[i]).expect("failed to read species-list file");
                for line in text.lines() {
                    let l = line.trim();
                    if !l.is_empty() {
                        species.push(l.to_string());
                    }
                }
            }
            "--evalue" => {
                i += 1;
                evalue_thresh = args[i].parse().expect("bad --evalue");
            }
            "--slop" => {
                i += 1;
                slop = args[i].parse().expect("bad --slop");
            }
            "--out-dir" => {
                i += 1;
                out_dir = PathBuf::from(&args[i]);
            }
            other => {
                species.push(other.to_string());
            }
        }
        i += 1;
    }

    if species.is_empty() {
        species = [
            "Sagina_nodosa",
            "Silene_flos_cuculi",
            "Myosotis_sylvatica",
            "Hydrocharis_morsus_ranae",
            "Cerastium_arvense",
            "Convallaria_majalis",
            "Myosotis_alpestris",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
    }

    let out_dir = if out_dir.is_relative() {
        root.join(out_dir)
    } else {
        out_dir
    };
    fs::create_dir_all(&out_dir).expect("failed to create out-dir");

    let mut fasta_out = String::new();
    let mut tsv_out = String::from(
        "species\torganelle\tgene\tunitig\tstart\tend\tstrand\tn_hsps_merged\tbest_evalue\tbest_score\tsource_run\tnote\n",
    );

    println!(
        "{:<28}{:<45}{:>10}{:>10}{:>10}",
        "species", "chosen_run", "n_genes", "n_frags", "bp_total"
    );

    for sp in &species {
        let species_dir = root.join("rev").join("out").join(sp);
        let candidates = find_candidates(&species_dir, evalue_thresh, slop);
        if candidates.is_empty() {
            println!("{:<28}{:<45}{:>10}{:>10}{:>10}", sp, "(no eligible candidate run)", 0, 0, 0);
            continue;
        }
        let best_i = match choose_best(&candidates) {
            Some(i) => i,
            None => {
                println!("{:<28}{:<45}{:>10}{:>10}{:>10}", sp, "(choose_best failed)", 0, 0, 0);
                continue;
            }
        };
        let chosen = &candidates[best_i];
        let seqs = &chosen.seqs;
        let fragments = merge_fragments(&chosen.clean_hits, slop);

        let mut n_genes: HashSet<&str> = HashSet::new();
        let mut bp_total: i64 = 0;

        for frag in &fragments {
            let seq = match seqs.get(&frag.unitig) {
                Some(s) => s,
                None => continue, // unitig referenced in annotation but absent from graph file
            };
            let start = frag.start.max(1) as usize;
            let end = (frag.end as usize).min(seq.len());
            if start > end || start == 0 {
                continue;
            }
            let window = &seq[(start - 1)..end];
            let extracted = if frag.strand == '-' {
                revcomp(window)
            } else {
                window.to_string()
            };

            let header = format!(
                "{}__{}__{}:{}-{}({})__src={}__fragmentary_partial_recovery",
                sp, frag.gene, frag.unitig, frag.start, frag.end, frag.strand, chosen.prefix
            );
            write_fasta_record(&mut fasta_out, &header, &extracted);

            tsv_out.push_str(&format!(
                "{}\tmito\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3e}\t{:.1}\t{}\t{}\n",
                sp,
                frag.gene,
                frag.unitig,
                frag.start,
                frag.end,
                frag.strand,
                frag.n_hsps,
                frag.best_evalue,
                frag.best_score,
                chosen.prefix,
                "fragmentary_partial_recovery;not_a_resolved_genome;rescued_from_pltd-classified_unitig"
            ));

            n_genes.insert(frag.gene.as_str());
            bp_total += extracted.len() as i64;
        }

        println!(
            "{:<28}{:<45}{:>10}{:>10}{:>10}",
            sp,
            chosen.prefix,
            n_genes.len(),
            fragments.len(),
            bp_total
        );

    }

    fs::write(out_dir.join("mito_fragments.fasta"), fasta_out).expect("failed to write fasta");
    fs::write(out_dir.join("mito_fragments.tsv"), tsv_out).expect("failed to write tsv");
    println!("\nWrote {}/mito_fragments.{{fasta,tsv}}", out_dir.display());
}

#!/usr/bin/env bash
set -euo pipefail

ROOT="/data/tol/data/darwin"
OUT="meta/file.txt"

mkdir -p meta

# Find all *.trim.fasta.gz and write: species<TAB>full_path
# Species is the directory right after the clade (dicots/monocots/...)
fd -L -t f --full-path -g '**/genomic_data/*/pacbio/fasta/*.trim.fasta.gz' \
  "${ROOT}"/{monocots,dicots,vascular-plants,non-vascular-plants} \
  | sort \
  | awk -F/ 'BEGIN{OFS="\t"}
      {
        # path layout: .../darwin/<clade>/<Species>/genomic_data/<sample>/pacbio/fasta/<file>
        clade=$6;
        sp=$7;
        print sp, $0
      }' > "${OUT}"

echo "[OK] wrote ${OUT}"
echo "[OK] lines: $(wc -l < "${OUT}")"


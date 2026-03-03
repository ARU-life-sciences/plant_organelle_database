#!/usr/bin/env bash
set -euo pipefail

############################################
# Required
############################################
LIST=""
META="meta/file.txt"

############################################
# Defaults (override via flags)
############################################
K=1001
S=31
COVERAGE=120
MAX_FASTA_GB=15
MEM_MB=100000
THREADS=20
QUEUE="normal"

OUTROOT="rev/out"
LOGDIR="rev/logs"
SUFFIX="rev"

############################################
# OATK resources
############################################
OATK="/software/team301/oatk/oatk"
NHMMSCAN="/software/team301/hmmer-3.4/src/nhmmscan"
MITO_FAM="/software/team301/oatk/acrogymnospermae_mito.fam"
PLTD_FAM="/data/tol/resources/OatkDB/v20230210/angiosperm_pltd.fam"

############################################
usage() {
  cat <<EOF
Usage:
  $0 --list species.txt [options]

Options:
  --coverage <int>     (default: 120)
  --max-gb <int>       (default: 15)
  --mem-mb <int>       (default: 100000)
  --threads <int>      (default: 20)
  --queue <name>       (default: normal)
  --suffix <tag>       (default: rev)
  --meta <file>        (default: meta/file.txt)

Example:
  $0 --list rev/lists/to_resubmit_memfix.txt \
     --coverage 120 --max-gb 15 --suffix memfix_c120
EOF
}
############################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    --list) LIST="$2"; shift 2;;
    --coverage) COVERAGE="$2"; shift 2;;
    --max-gb) MAX_FASTA_GB="$2"; shift 2;;
    --mem-mb) MEM_MB="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --queue) QUEUE="$2"; shift 2;;
    --suffix) SUFFIX="$2"; shift 2;;
    --meta) META="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -n "$LIST" ]] || { echo "[ERROR] --list required"; exit 1; }
[[ -f "$LIST" ]] || { echo "[ERROR] list not found: $LIST"; exit 1; }
[[ -f "$META" ]] || { echo "[ERROR] meta not found: $META"; exit 1; }

mkdir -p "$OUTROOT" "$LOGDIR"

MAX_SIZE=$((MAX_FASTA_GB * 1024 * 1024 * 1024))

submitted=0

############################################
# Main loop
############################################
while IFS= read -r sp; do
  [[ -n "$sp" ]] || continue

  # Pull FASTA paths for this species
  mapfile -t files < <(awk -F'\t' -v s="$sp" '$1==s {print $2}' "$META")

  if (( ${#files[@]} == 0 )); then
    echo "[WARN] $sp: no FASTA entries in meta"
    continue
  fi

  ##########################################
  # Select FASTAs up to size cap
  ##########################################
  selected=()
  total=0
  tmp="$(mktemp)"
  trap 'rm -f "$tmp"' RETURN

  for f in "${files[@]}"; do
    [[ -f "$f" ]] || continue
    sz=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f")
    printf "%s\t%s\n" "$sz" "$f" >> "$tmp"
  done

  while IFS=$'\t' read -r sz f; do
    (( sz > MAX_SIZE )) && continue
    (( total + sz > MAX_SIZE )) && break
    selected+=("$f")
    total=$((total + sz))
  done < <(sort -n "$tmp")

  rm -f "$tmp"
  trap - RETURN

  if (( ${#selected[@]} == 0 )); then
    echo "[WARN] $sp: nothing selected under ${MAX_FASTA_GB}GB"
    continue
  fi

  ##########################################
  # Output paths
  ##########################################
  outdir="${OUTROOT}/${sp}"
  mkdir -p "$outdir"

  prefix="${sp}.k${K}.s${S}.c${COVERAGE}.${SUFFIX}"

  # Skip if already exists (avoid duplicates)
  if ls "${outdir}/${prefix}"*.mito.gfa* >/dev/null 2>&1; then
    echo "[SKIP] $sp: already has ${prefix}"
    continue
  fi

  ##########################################
  # Build oatk command
  ##########################################
  cmd=(
    "$OATK"
    -k "$K"
    -c "$COVERAGE"
    -t "$THREADS"
    --nhmmscan "$NHMMSCAN"
    -m "$MITO_FAM"
    -p "$PLTD_FAM"
    -o "${outdir}/${prefix}"
  )
  cmd+=( "${selected[@]}" )

  echo "[SUBMIT] $sp  files=${#selected[@]}  total_gb=$(python3 - <<PY
print(round(${total}/(1024**3),2))
PY
)  cov=$COVERAGE"

  bsub -n "$THREADS" -q "$QUEUE" \
    -R"span[hosts=1] select[mem>${MEM_MB}] rusage[mem=${MEM_MB}]" \
    -M"$MEM_MB" \
    -o "${LOGDIR}/${sp}.${SUFFIX}.out" \
    -e "${LOGDIR}/${sp}.${SUFFIX}.err" \
    "$(printf "%q " "${cmd[@]}")"

  ((submitted+=1))

done < "$LIST"

echo "[DONE] submitted=$submitted"

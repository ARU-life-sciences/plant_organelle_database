#!/usr/bin/env bash
set -euo pipefail

# to be run from the root dir

META="meta/file.txt"

# oatk parameters
K=1001
S=31
COVERAGE="${1:-80}"          # default coverage if not provided
MAX_FASTA_GB="${2:-16}"     # max total FASTA input per species (GB), default 100

# oatk resources (as per your command)
OATK="/software/team301/oatk/oatk"
NHMMSCAN="/software/team301/hmmer-3.4/src/nhmmscan"
MITO_FAM="/software/team301/oatk/acrogymnospermae_mito.fam"
PLTD_FAM="/data/tol/resources/OatkDB/v20230210/angiosperm_pltd.fam"

# outputs
OUTROOT="data/oatk_assemblies"   # you can change to data/mito later if you want
LOGDIR="logs"
mkdir -p "${OUTROOT}" "${LOGDIR}"

if [[ ! -f "${META}" ]]; then
  echo "[ERROR] Missing ${META}. Run src/01_make_meta.sh first." >&2
  exit 1
fi

MAX_SIZE=$((MAX_FASTA_GB * 1024 * 1024 * 1024))

# Build a temporary list per species via awk, then loop species
# This prints: species<TAB>file1<TAB>file2...
tmp="$(mktemp)"
trap 'rm -f "${tmp}"' EXIT

awk -F'\t' '
  { sp=$1; path=$2; files[sp]=files[sp] "\t" path }
  END { for (sp in files) print sp files[sp] }
' "${META}" | sort > "${tmp}"

submitted=0
skipped=0

while IFS=$'\t' read -r sp rest; do
  [[ -n "${sp}" ]] || continue

    # Skip ONLY if both mito AND plastid are already present
  mito_done=0
  pltd_done=0

  if ls "data/mito/${sp}"/*.mito.gfa >/dev/null 2>&1; then
    mito_done=1
  fi
  if ls "data/plastid/${sp}"/*.pltd.gfa >/dev/null 2>&1; then
    pltd_done=1
  fi

  if (( mito_done == 1 && pltd_done == 1 )); then
    ((skipped+=1))
    continue
  fi

  # Collect files (all remaining fields)
  # shellcheck disable=SC2206
  files=( $rest )

  # Filter by total size budget, but keep it deterministic:
  # sort by size (small->large), then take until MAX_SIZE
  # (more stable than shuf; easier to debug/teach)
  selected=()
  total=0

  sizes_tmp="$(mktemp)"
  for f in "${files[@]}"; do
    [[ -f "$f" ]] || continue
    sz="$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f")"
    printf "%s\t%s\n" "$sz" "$f" >> "${sizes_tmp}"
  done

  # IMPORTANT: use process substitution so the while runs in the current shell
  while IFS=$'\t' read -r sz f; do
    # skip single files bigger than limit
    if (( sz > MAX_SIZE )); then
      continue
    fi
    # stop when we'd exceed the limit
    if (( total + sz > MAX_SIZE )); then
      break
    fi
    selected+=("$f")
    total=$((total + sz))
  done < <(sort -n "${sizes_tmp}")

  rm -f "${sizes_tmp}"


  if (( ${#selected[@]} == 0 )); then
    echo "[WARN] No files selected for ${sp} (maybe all too big?)"
    continue
  fi

  # Prepare output dir
  outdir="${OUTROOT}/${sp}"
  mkdir -p "${outdir}"

  prefix="${sp}.k${K}.s${S}.c${COVERAGE}"

  # Size-tiered resource request. Calibrated 2026-11-02 against LSF's own
  # "Resource usage summary" for 15 already-completed species: actual max
  # memory tracks selected-input size at roughly 5-9x (e.g. 1.81GB->11.6GB,
  # 7.11GB->53.7GB, 14.51GB->98.8GB used), while a blanket 100GB/20-core ask
  # for every species (median actual use ~40GB) was starving jobs of the
  # increasingly scarce "whole free host with 20 cores + 100GB" combination
  # on a busy shared normal queue.
  #
  # Top tier bumped 100GB->130GB on 2026-11-05: of the first 44 species
  # submitted at >10GB input (then capped at 100GB), 15 (34%) were
  # TERM_MEMLIMIT-killed right around 100-103GB with no size-within-tier
  # pattern distinguishing them from the 29 that finished under 100GB -
  # it's real variance (repeat content etc.), not something a size
  # sub-split can predict, so the cap itself needed headroom instead.
  #
  # Compare in bytes, not GB-truncated-to-int, so e.g. a 2.9GB input isn't
  # silently rounded down into the <=2GB tier.
  if (( total <= 2*1024*1024*1024 )); then
    cores=6;  mbMem=24000
  elif (( total <= 5*1024*1024*1024 )); then
    cores=10; mbMem=48000
  elif (( total <= 10*1024*1024*1024 )); then
    cores=16; mbMem=80000
  else
    cores=20; mbMem=130000
  fi

  # Build oatk command
  cmd=(
    "${OATK}"
    -k "${K}"
    -c "${COVERAGE}"
    -t "${cores}"
    --nhmmscan "${NHMMSCAN}"
    -m "${MITO_FAM}"
    -p "${PLTD_FAM}"
    -o "${outdir}/${prefix}"
  )
  cmd+=( "${selected[@]}" )

  if (( mito_done == 1 && pltd_done == 0 )); then
    echo "[INFO] ${sp}: mito present, plastid missing -> will run"
  elif (( mito_done == 0 && pltd_done == 1 )); then
    echo "[INFO] ${sp}: plastid present, mito missing -> will run"
  else
    echo "[INFO] ${sp}: neither present -> will run"
  fi

  # Submit
  bsub -n "${cores}" -q normal \
    -R"span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" \
    -M"${mbMem}" \
    -o "${LOGDIR}/${sp}.out" \
    -e "${LOGDIR}/${sp}.err" \
    "$(printf "%q " "${cmd[@]}")"

  echo "[SUBMITTED] ${sp}  files=${#selected[@]}  total_bytes=${total}  tier=(cores=${cores},mem=${mbMem}MB)"
  ((submitted+=1))

done < "${tmp}"

echo "[OK] submitted=${submitted} skipped_done=${skipped}"


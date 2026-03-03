# make a presence/absence matrix of genes across species
awk -F'\t' '
BEGIN { OFS="\t" }
# Skip comments and unwanted features
/^#/ { next }

{
  # species = directory name; adjust if yours differs
  split(FILENAME, p, "/"); sp = p[1]

  # gene name from column 4; drop anything after a space (e.g. exon info)
  g = $4
  if (g ~ /^(trn|rrn|#)/) next
  sub(/ .*/, "", g)

  Species[sp] = 1
  Genes[g] = 1
  seen[sp SUBSEP g] = 1    # presence indicator
}

END {
  # sort genes and species for stable columns/rows
  ng = asorti(Genes, Gs)
  ns = asorti(Species, Ss)

  # header
  printf "species"
  for (i=1; i<=ng; i++) printf OFS "%s", Gs[i]
  printf "\n"

  # rows
  for (i=1; i<=ns; i++) {
    sp = Ss[i]
    printf "%s", sp
    for (j=1; j<=ng; j++) {
      g = Gs[j]
      printf OFS ((sp SUBSEP g) in seen ? 1 : 0)
    }
    printf "\n"
  }
}
' */*.pltd.bed > presence_absence.tsv

# Calculate percentage of species with each gene
awk -F'\t' '
NR==1 { for (i=2;i<=NF;i++) gene[i]=$i; next }
      { n++; for (i=2;i<=NF;i++) sum[i]+=($i>0) }
END   {
  print "Gene\tpct_complete"
  for (i=2;i<=NF;i++) printf "%s\t%.13f\n", gene[i], 100*sum[i]/n
}
' presence_absence.tsv | { read -r h; echo "$h"; sort -k2,2nr; } > gene_pct.tsv

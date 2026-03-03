# Get the plastid data

Where the plastid data occurs in the primary assemblies (this is most), it was just extracted:

```bash
/software/team301/seqkit grep -r -p 'Pltd' <fasta> > <{spp}_pltd.fasta>
```

Otherwise I've had to assemble using `oatk` and extract the linearised plastid genome. This was achieved using a modified version of the `../../src/memlim.bash`.

## Collect data on genes

All genes present. There are 80.

```bash
cut -f4 */*.pltd.bed | sort | uniq -c | sort -k1,1nr | rg -v "trn|rrn|#" | wc -l
```

Need to make a matrix of 250 x 80, presence/absence of the gene.

```bash
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
' */*.pltd.bed
```

## Sort these species out

./Azolla_filiculoides   18
./Adiantum_capillus_veneris     13
./Ginkgo_biloba 13
./Polystichum_setiferum 13
./Pteridium_aquilinum   13
./Humulus_lupulus       11
./Asplenium_adiantum_nigrum     8
./Asplenium_ceterach    8
./Asplenium_marinum     8
./Asplenium_scolopendrium       8

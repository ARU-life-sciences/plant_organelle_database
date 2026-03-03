#!/usr/bin/env bash
set -euo pipefail

src_dir="../oatk_assemblies"
dst_dir="../data/oatk_assemblies_pltd"

mkdir -p "$dst_dir"

# Find all plastid files
fd -e txt -e bed -e fasta -e gfa -g "*.pltd*" "$src_dir" | while read -r file; do
    species=$(basename "$(dirname "$file")")   # species name from parent dir
    dest="${dst_dir}/${species}"
    mkdir -p "$dest"
    echo "Moving $file → $dest/"
    mv "$file" "$dest/"
done


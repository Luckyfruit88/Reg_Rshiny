#!/bin/bash
# Make a gene-searchable copy of Reactome's Ensembl2Reactome_All_Levels.txt for
# the Variant interpretation tab: human rows keyed by Ensembl gene ID, in byte
# order of the first column so `look` finds one gene by binary search.
#
#   scripts/build_pathway_index.sh Ensembl2Reactome_All_Levels.txt DEST.tsv
#
# Columns are unchanged: gene ID, pathway ID, URL, pathway name, evidence code,
# species. The source is only read; DEST is written atomically.
set -euo pipefail
if [ "$#" -ne 2 ]; then echo "usage: $0 Ensembl2Reactome_All_Levels.txt DEST.tsv" >&2; exit 2; fi
src=$1; dest=$2
[ -r "$src" ] || { echo "cannot read $src" >&2; exit 1; }
if [ -s "$dest" ] && [ "$dest" -nt "$src" ]; then echo "up to date: $dest"; exit 0; fi
mkdir -p "$(dirname "$dest")"
work=$(mktemp -d "$(dirname "$dest")/.pathway_index.XXXXXX")
trap 'rm -rf "$work"' EXIT
LC_ALL=C awk -F'\t' '$1 ~ /^ENSG/ && $6 == "Homo sapiens"' "$src" | LC_ALL=C sort -u -S 1G -T "$work" > "$work/sorted.tsv"
[ -s "$work/sorted.tsv" ] || { echo "no human Ensembl gene row in $src" >&2; exit 1; }
chmod 0644 "$work/sorted.tsv"
mv "$work/sorted.tsv" "$dest"
echo "rows: $(wc -l < "$dest")  genes: $(cut -f1 "$dest" | uniq | wc -l)  $dest"

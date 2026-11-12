#!/bin/bash
# Make a position-searchable copy of a QTL summary-statistics CSV for the
# Variant interpretation tab. The copy is the same file with its data rows in
# byte order of the first column (SNP id such as chr1:12345:A:G:rs1), which lets
# `look` find one variant by binary search in well under a second.
#
#   scripts/build_qtl_index.sh SOURCE.csv DEST.csv
#
# The source is only read. DEST is written atomically; an existing DEST is kept
# unless it is older than the source. Needs about the source's size in free
# space twice over while it runs.
set -euo pipefail
if [ "$#" -ne 2 ]; then echo "usage: $0 SOURCE.csv DEST.csv" >&2; exit 2; fi
src=$1; dest=$2
[ -r "$src" ] || { echo "cannot read $src" >&2; exit 1; }
if [ -s "$dest" ] && [ "$dest" -nt "$src" ]; then echo "up to date: $dest"; exit 0; fi
mkdir -p "$(dirname "$dest")"
work=$(mktemp -d "$(dirname "$dest")/.qtl_index.XXXXXX")
trap 'rm -rf "$work"' EXIT
head -n 1 "$src" > "$work/sorted.csv"
case "$(head -n 1 "$src")" in SNP,*) ;; *) echo "first column must be SNP" >&2; exit 1;; esac
# A final line without a newline is an interrupted copy of a row: leave it out.
tail -n +2 "$src" | awk 'RT == "\n"' RS='\n' | LC_ALL=C sort -S 2G -T "$work" >> "$work/sorted.csv"
chmod 0644 "$work/sorted.csv"
mv "$work/sorted.csv" "$dest"
echo "rows: $(($(wc -l < "$dest") - 1))  $dest"

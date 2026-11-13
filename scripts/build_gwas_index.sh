#!/bin/bash
# Make a position-indexed copy of a GWAS Catalog associations file (the
# "all associations v1.0.2" TSV, GRCh38 positions) for the Variant
# interpretation tab.
#
#   scripts/build_gwas_index.sh SOURCE.tsv DEST.tsv.gz
#
# DEST holds one row per reported variant position with the columns named in
# its header line, sorted by chromosome and position, bgzip-compressed and
# tabix-indexed (DEST.tbi). The source is only read. Rows without a mapped
# position are left out; a row reporting several variants becomes one row each.
# Needs bgzip and tabix (htslib); set BGZIP / TABIX when they are not on PATH.
set -euo pipefail
if [ "$#" -ne 2 ]; then echo "usage: $0 SOURCE.tsv DEST.tsv.gz" >&2; exit 2; fi
src=$1; dest=$2
bgzip=${BGZIP:-bgzip}; tabix=${TABIX:-tabix}
[ -r "$src" ] || { echo "cannot read $src" >&2; exit 1; }
if [ -s "$dest" ] && [ -s "$dest.tbi" ] && [ "$dest" -nt "$src" ]; then echo "up to date: $dest"; exit 0; fi
mkdir -p "$(dirname "$dest")"
work=$(mktemp -d "$(dirname "$dest")/.gwas_index.XXXXXX")
trap 'rm -rf "$work"' EXIT
LC_ALL=C awk -F'\t' -v OFS='\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    n = split("CHR_ID|CHR_POS|SNPS|STRONGEST SNP-RISK ALLELE|DISEASE/TRAIT|MAPPED_TRAIT|P-VALUE|PVALUE_MLOG|OR or BETA|95% CI (TEXT)|RISK ALLELE FREQUENCY|MAPPED_GENE|PUBMEDID|FIRST AUTHOR|DATE|STUDY ACCESSION|INITIAL SAMPLE SIZE", name, "|")
    for (i = 1; i <= n; i++) if (!(name[i] in col)) { print "missing column: " name[i] > "/dev/stderr"; exit 1 }
    next
  }
  {
    nc = split($col["CHR_ID"], chrom, /[;x] ?/); np = split($col["CHR_POS"], pos, /[;x] ?/)
    ns = split($col["SNPS"], snp, /[;x,] ?/); na = split($col["STRONGEST SNP-RISK ALLELE"], allele, /[;x,] ?/)
    if (nc != np) next
    for (k = 1; k <= np; k++) {
      if (pos[k] !~ /^[0-9]+$/ || chrom[k] !~ /^([0-9]+|X|Y|MT)$/) continue
      c = (chrom[k] == "MT") ? "chrM" : "chr" chrom[k]
      line = c OFS pos[k] OFS (ns == np ? snp[k] : $col["SNPS"]) OFS (na == np ? allele[k] : $col["STRONGEST SNP-RISK ALLELE"])
      for (i = 5; i <= n; i++) line = line OFS $col[name[i]]
      print line
    }
  }' "$src" | LC_ALL=C sort -t "$(printf '\t')" -k1,1 -k2,2n -S 1G -T "$work" > "$work/body.tsv"
[ -s "$work/body.tsv" ] || { echo "no positioned association in $src" >&2; exit 1; }
{
  printf '#chrom\tpos\trsid\trisk_allele\ttrait\tmapped_trait\tp_value\tp_mlog\teffect\tci_text\trisk_af\tmapped_gene\tpubmed\tfirst_author\tdate\tstudy_accession\tinitial_sample\n'
  cat "$work/body.tsv"
} | "$bgzip" -c > "$work/out.tsv.gz"
"$tabix" -s 1 -b 2 -e 2 -c '#' "$work/out.tsv.gz"
chmod 0644 "$work/out.tsv.gz" "$work/out.tsv.gz.tbi"
mv "$work/out.tsv.gz.tbi" "$dest.tbi"
mv "$work/out.tsv.gz" "$dest"
touch "$dest.tbi"
echo "rows: $(wc -l < "$work/body.tsv")  $dest"

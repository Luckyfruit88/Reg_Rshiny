#!/bin/bash
# Make a position-indexed copy of cis-pQTL summary statistics (GRCh38) for the
# Variant interpretation tab.
#
#   scripts/build_pqtl_index.sh regenie SOURCE_DIR DEST.tsv.gz
#   scripts/build_pqtl_index.sh table   SOURCE.txt.gz DEST.tsv.gz
#
# regenie: SOURCE_DIR holds one <protein>.csv.gz per protein, as written by
#   REGENIE with GRCh38 positions added (columns CHROM, ALLELE0, ALLELE1,
#   A1FREQ, INFO, N, BETA, SE, LOG10P, rsid, POS38). ALLELE1 is the effect
#   allele. The protein name is the file name.
# table: SOURCE is one tab-separated file with the columns SNP, EA, OA, CHR,
#   BP (GRCh38), EAF, OlinkID, Gene, BETA, SE, P.
#
# DEST holds one row per variant and protein with the columns named in its
# header line, sorted by chromosome and position, bgzip-compressed and
# tabix-indexed (DEST.tbi). Sources are only read. Needs bgzip and tabix
# (htslib); set BGZIP / TABIX when they are not on PATH.
set -euo pipefail
if [ "$#" -ne 3 ]; then echo "usage: $0 regenie|table SOURCE DEST.tsv.gz" >&2; exit 2; fi
mode=$1; src=$2; dest=$3
bgzip=${BGZIP:-bgzip}; tabix=${TABIX:-tabix}
[ -r "$src" ] || { echo "cannot read $src" >&2; exit 1; }
mkdir -p "$(dirname "$dest")"
work=$(mktemp -d "${TMPDIR:-$(dirname "$dest")}/.pqtl_index.XXXXXX")
trap 'rm -rf "$work"' EXIT
tab=$(printf '\t')

case "$mode" in
  regenie)
    [ -d "$src" ] || { echo "not a directory: $src" >&2; exit 1; }
    for f in "$src"/*.csv.gz; do
      [ -r "$f" ] || continue
      protein=$(basename "$f" .csv.gz)
      gzip -dc "$f" | LC_ALL=C awk -F',' -v OFS='\t' -v protein="$protein" '
        NR == 1 {
          for (i = 1; i <= NF; i++) col[$i] = i
          n = split("CHROM|ALLELE0|ALLELE1|A1FREQ|INFO|N|BETA|SE|LOG10P|rsid|POS38", name, "|")
          for (i = 1; i <= n; i++) if (!(name[i] in col)) { print "missing column " name[i] " in " protein > "/dev/stderr"; exit 1 }
          next
        }
        $col["POS38"] ~ /^[0-9]+$/ && $col["BETA"] != "NA" {
          c = $col["CHROM"]; if (c == "23") c = "X"
          print "chr" c, $col["POS38"], $col["ALLELE0"], $col["ALLELE1"], $col["rsid"], protein, ".",
                $col["A1FREQ"], $col["INFO"], $col["N"], $col["BETA"], $col["SE"], $col["LOG10P"]
        }'
    done > "$work/unsorted.tsv"
    ;;
  table)
    gzip -dc "$src" | LC_ALL=C awk -F'\t' -v OFS='\t' '
      NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        n = split("SNP|EA|OA|CHR|BP|EAF|OlinkID|Gene|BETA|SE|P", name, "|")
        for (i = 1; i <= n; i++) if (!(name[i] in col)) { print "missing column " name[i] > "/dev/stderr"; exit 1 }
        next
      }
      $col["BP"] ~ /^[0-9]+$/ && $col["BETA"] != "NA" && $col["P"] != "NA" {
        c = $col["CHR"]; if (c == "23") c = "X"
        p = $col["P"] + 0
        mlog = (p > 0) ? -log(p) / log(10) : 400
        print "chr" c, $col["BP"], $col["OA"], $col["EA"], $col["SNP"], $col["Gene"], $col["OlinkID"],
              $col["EAF"], ".", ".", $col["BETA"], $col["SE"], sprintf("%.5g", mlog)
      }' > "$work/unsorted.tsv"
    ;;
  *) echo "unknown mode: $mode" >&2; exit 2 ;;
esac

[ -s "$work/unsorted.tsv" ] || { echo "no usable row in $src" >&2; exit 1; }
LC_ALL=C sort -t "$tab" -k1,1 -k2,2n -S "${SORT_MEMORY:-4G}" -T "$work" "$work/unsorted.tsv" > "$work/body.tsv"
{
  printf '#chrom\tpos\tother_allele\teffect_allele\trsid\tprotein\tassay_id\teffect_af\tinfo\tn\tbeta\tse\tp_mlog\n'
  cat "$work/body.tsv"
} | "$bgzip" -c > "$work/out.tsv.gz"
"$tabix" -s 1 -b 2 -e 2 -c '#' "$work/out.tsv.gz"
chmod 0644 "$work/out.tsv.gz" "$work/out.tsv.gz.tbi"
cp "$work/out.tsv.gz" "$dest.part" && mv "$dest.part" "$dest"
cp "$work/out.tsv.gz.tbi" "$dest.tbi.part" && mv "$dest.tbi.part" "$dest.tbi"
touch "$dest.tbi"
echo "rows: $(wc -l < "$work/body.tsv")  proteins: $(cut -f6 "$work/body.tsv" | LC_ALL=C sort -u | wc -l)  $dest"

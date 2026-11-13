#!/usr/bin/env python3
"""Make the Olink assay -> gene table used with the cis-pQTL copies.

    scripts/build_protein_map.py --out MAP.tsv [--ukbppp olink_protein_map_3k_v1.tsv] [--olink-ht olink-exht-annot.csv]

--ukbppp   the UKB-PPP protein map (tab-separated: Assay, olink_target_fullname, OlinkID,
           HGNC.symbol, ensembl_id, chr, gene_start, gene_end; GRCh38)
--olink-ht the Olink Explore HT assay annotation (comma-separated: OlinkID, Assay,
           Protein name, Gene name, Chr, Start, End, ENSG; GRCh38)

MAP gets one row per assay and gene with the columns
    assay  olink_id  protein_name  gene  gene_id  chrom  gene_start  gene_end  panel
The source files are only read.
"""
import argparse
import csv
import os
import sys

COLUMNS = ["assay", "olink_id", "protein_name", "gene", "gene_id", "chrom", "gene_start", "gene_end", "panel"]


def chrom_name(value):
    value = (value or "").strip()
    if not value:
        return ""
    return value if value.startswith("chr") else "chr" + value


def whole(value):
    value = (value or "").strip()
    return value if value.isdigit() else ""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ukbppp")
    parser.add_argument("--olink-ht")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    if not args.ukbppp and not args.olink_ht:
        sys.exit("give --ukbppp, --olink-ht or both")
    rows = []
    if args.ukbppp:
        with open(args.ukbppp, newline="", encoding="utf-8-sig") as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                rows.append([r.get("Assay", ""), r.get("OlinkID", ""), r.get("olink_target_fullname", ""), r.get("HGNC.symbol", ""),
                             r.get("ensembl_id", ""), chrom_name(r.get("chr")), whole(r.get("gene_start")), whole(r.get("gene_end")),
                             "Olink Explore 3072"])
    if args.olink_ht:
        with open(args.olink_ht, newline="", encoding="utf-8-sig") as handle:
            for r in csv.DictReader(handle):
                gene = (r.get("Gene name") or "").strip()
                if not gene:
                    continue
                rows.append([r.get("Assay", ""), r.get("OlinkID", ""), r.get("Protein name", ""), gene, (r.get("ENSG") or "").strip(),
                             chrom_name(r.get("Chr")), whole(r.get("Start")), whole(r.get("End")), "Olink Explore HT"])
    clean = [[field.replace("\t", " ").strip() for field in row] for row in rows if row[0] or row[3]]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out + ".part", "w", encoding="utf-8") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in clean:
            out.write("\t".join(row) + "\n")
    os.chmod(args.out + ".part", 0o644)
    os.replace(args.out + ".part", args.out)
    print(f"assay rows: {len(clean)}  {args.out}")


if __name__ == "__main__":
    main()

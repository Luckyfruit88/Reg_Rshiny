#!/usr/bin/env python3
"""Make a gene-sorted table of the conditions ClinVar ties to each gene, for the Variant interpretation tab.

    scripts/build_gene_disease_index.py --clinvar clinvar.vcf.gz --out GENE_DISEASE.tsv

For every gene (GENEINFO symbol) the table counts the ClinVar variants classified pathogenic or likely
pathogenic (conflicting classifications are left out) and lists the conditions they are reported for, most
variants first:

    gene  pathogenic_variants  conditions_total  conditions

`conditions` holds up to eight entries "condition (number of variants)" separated by "; ". Rows are sorted by
gene symbol in C order so that `look` finds a gene. The VCF is only read.
"""
import argparse
import collections
import gzip
import os
import re
import sys

SKIP = {"not_provided", "not_specified", "", "."}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clinvar", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-conditions", type=int, default=8)
    args = parser.parse_args()
    variants = collections.Counter()
    conditions = collections.defaultdict(collections.Counter)
    file_date = ""
    opener = gzip.open if args.clinvar.endswith(".gz") else open
    with opener(args.clinvar, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                if line.startswith("##fileDate="):
                    file_date = line.strip().split("=", 1)[1]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            significance = info.get("CLNSIG", "").lower()
            if "pathogenic" not in significance or "conflicting" in significance:
                continue
            genes = {entry.split(":")[0] for entry in info.get("GENEINFO", "").split("|") if entry}
            names = {name for name in re.split(r"[|]", info.get("CLNDN", "")) if name not in SKIP}
            for gene in genes:
                if not gene or "\t" in gene:
                    continue
                variants[gene] += 1
                for name in names:
                    conditions[gene][name.replace("_", " ")] += 1
    if not variants:
        sys.exit(f"no pathogenic or likely pathogenic record with a gene in {args.clinvar}")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out + ".part", "w", encoding="utf-8") as out:
        for gene in sorted(variants, key=lambda g: g.encode("utf-8")):
            top = sorted(conditions[gene].items(), key=lambda item: (-item[1], item[0]))
            listed = "; ".join(f"{name} ({count})" for name, count in top[: args.max_conditions])
            out.write(f"{gene}\t{variants[gene]}\t{len(top)}\t{listed}\n")
    os.chmod(args.out + ".part", 0o644)
    os.replace(args.out + ".part", args.out)
    print(f"genes: {len(variants)}  ClinVar file date: {file_date or 'unknown'}  {args.out}")


if __name__ == "__main__":
    main()

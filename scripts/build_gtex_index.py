#!/usr/bin/env python3
"""Make position-indexed copies of GTEx cis-eQTL results for the Variant interpretation tab.

    scripts/build_gtex_index.py --source DIR --pairs PAIRS.tsv.gz --genes GENES.tsv
                                [--bgzip bgzip] [--tabix tabix] [--work DIR] [--memory 8G]

DIR holds the GTEx Analysis release files, per tissue
    <Tissue>.<release>.eQTLs.signif_pairs.parquet   significant variant-gene pairs
    <Tissue>.<release>.eGenes.txt.gz                one row per tested gene

PAIRS gets one row per significant variant-gene pair and tissue, sorted by
chromosome and position, bgzip-compressed and tabix-indexed (PAIRS.tbi).
GENES gets one row per tested gene and tissue, sorted by gene id (without the
version) so that `look` finds a gene: whether the gene was tested in the tissue,
its q-value, its nominal p-value threshold and its lead variant.
The source files are only read. Needs pyarrow, sort, bgzip and tabix.
"""
import argparse
import csv
import glob
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

PAIR_COLUMNS = ["chrom", "pos", "ref", "alt", "gene_id", "tissue", "slope", "slope_se", "pval_nominal",
                "af", "pval_nominal_threshold", "start_distance"]
GENE_COLUMNS = ["gene", "gene_id", "gene_name", "biotype", "tissue", "qval", "pval_nominal_threshold",
                "lead_variant", "lead_slope", "lead_slope_se", "lead_pval_nominal", "lead_af", "afc"]


def tissue_of(path):
    return os.path.basename(path).split(".")[0]


def write_pairs(path, out):
    tissue = tissue_of(path)
    source = pq.ParquetFile(path)
    rows = 0
    for index in range(source.metadata.num_row_groups):
        table = source.read_row_group(index, columns=["phenotype_id", "variant_id", "start_distance", "af",
                                                      "pval_nominal", "slope", "slope_se", "pval_nominal_threshold"])
        parts = pc.split_pattern(table["variant_id"], "_")
        columns = [
            pc.list_element(parts, 0), pc.cast(pc.list_element(parts, 1), pa.int64()),
            pc.list_element(parts, 2), pc.list_element(parts, 3), table["phenotype_id"],
            pa.array([tissue] * table.num_rows, pa.string()), table["slope"], table["slope_se"],
            table["pval_nominal"], table["af"], table["pval_nominal_threshold"], table["start_distance"],
        ]
        pacsv.write_csv(pa.table(columns, names=PAIR_COLUMNS), out,
                        pacsv.WriteOptions(include_header=False, delimiter="\t", quoting_style="none"))
        rows += table.num_rows
    return rows


def write_genes(path, out):
    tissue = tissue_of(path)
    rows = 0
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_id = row["gene_id"]
            out.write("\t".join([
                gene_id.split(".")[0], gene_id, row.get("gene_name", ""), row.get("biotype", ""), tissue,
                row.get("qval", ""), row.get("pval_nominal_threshold", ""), row.get("variant_id", ""),
                row.get("slope", ""), row.get("slope_se", ""), row.get("pval_nominal", ""), row.get("af", ""),
                row.get("afc", ""),
            ]) + "\n")
            rows += 1
    return rows


def run(command, **kwargs):
    subprocess.run(command, check=True, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--bgzip", default="bgzip")
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--work", default=None)
    parser.add_argument("--memory", default="8G")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    pair_files = sorted(glob.glob(os.path.join(args.source, "*.eQTLs.signif_pairs.parquet")))
    gene_files = sorted(glob.glob(os.path.join(args.source, "*.eGenes.txt.gz")))
    if not pair_files or not gene_files:
        sys.exit(f"no GTEx eQTL files in {args.source}")
    release = sorted({re.sub(r"^[^.]+\.([^.]+)\..*$", r"\1", os.path.basename(p)) for p in pair_files})
    for target in (args.pairs, args.genes):
        os.makedirs(os.path.dirname(os.path.abspath(target)), exist_ok=True)
    work = tempfile.mkdtemp(prefix=".gtex_index.", dir=args.work or os.environ.get("TMPDIR") or os.path.dirname(os.path.abspath(args.pairs)))
    env = dict(os.environ, LC_ALL="C")
    try:
        unsorted = os.path.join(work, "pairs.unsorted.tsv")
        total = 0
        with open(unsorted, "wb") as out:
            for path in pair_files:
                rows = write_pairs(path, out)
                total += rows
                print(f"{tissue_of(path)}\t{rows}", flush=True)
        body = os.path.join(work, "pairs.tsv")
        run(["sort", "-t", "\t", "-k1,1", "-k2,2n", "-S", args.memory, "--parallel", str(args.threads), "-T", work,
             "-o", body, unsorted], env=env)
        os.remove(unsorted)
        packed = os.path.join(work, "pairs.tsv.gz")
        with open(packed, "wb") as out:
            writer = subprocess.Popen([args.bgzip, "-@", str(args.threads), "-c"], stdin=subprocess.PIPE, stdout=out)
            writer.stdin.write(("#" + "\t".join(PAIR_COLUMNS) + "\n").encode())
            with open(body, "rb") as handle:
                shutil.copyfileobj(handle, writer.stdin, 1 << 24)
            writer.stdin.close()
            if writer.wait() != 0:
                sys.exit("bgzip failed")
        os.remove(body)
        run([args.tabix, "-s", "1", "-b", "2", "-e", "2", "-c", "#", packed])

        genes_unsorted = os.path.join(work, "genes.unsorted.tsv")
        gene_rows = 0
        with open(genes_unsorted, "w") as out:
            for path in gene_files:
                gene_rows += write_genes(path, out)
        genes_body = os.path.join(work, "genes.tsv")
        run(["sort", "-t", "\t", "-k1,1", "-k5,5", "-S", "1G", "-T", work, "-o", genes_body, genes_unsorted], env=env)

        for built, target in ((packed, args.pairs), (packed + ".tbi", args.pairs + ".tbi"), (genes_body, args.genes)):
            shutil.copyfile(built, target + ".part")
            os.chmod(target + ".part", 0o644)
            os.replace(target + ".part", target)
        os.utime(args.pairs + ".tbi")
        print(f"release: {','.join(release)}  tissues: {len(pair_files)}  pairs: {total}  gene rows: {gene_rows}")
        print(f"gene columns: {' '.join(GENE_COLUMNS)}")
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()

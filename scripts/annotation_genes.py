#!/usr/bin/env python3
"""Gene records of a GENCODE index by stable Ensembl gene ID. Read-only.

Kept apart from build_annotation.py: prepared indexes are keyed by that
builder's checksum, so a lookup added here does not make anyone rebuild.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_annotation import open_readonly, source_state  # noqa: E402


def genes_by_id(db, ids):
    """Gene records for stable Ensembl gene IDs; the version suffix is ignored."""
    ids = sorted({part.split(".")[0] for part in ids.split(",") if part})
    if len(ids) > 2000 or any(not re.match(r"^[A-Za-z0-9_]+$", gene) for gene in ids):
        raise ValueError("Give at most 2,000 plain gene identifiers.")
    before = source_state(db)
    con, metadata = open_readonly(db)
    try:
        rows = []
        for gene in ids:
            # "." sorts just below "/": this range is every version of one stable ID.
            rows.extend(dict(row) for row in con.execute(
                "SELECT * FROM genes WHERE gene_id>=? AND gene_id<? ORDER BY gene_id", (gene + ".", gene + "/")))
    finally:
        con.close()
    if source_state(db) != before:
        raise ValueError("Annotation index changed while being read.")
    return {"genes": rows, "release": metadata["release"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True)
    parser.add_argument("--ids", required=True, help="Comma-separated Ensembl gene IDs, with or without version.")
    args = parser.parse_args()
    try:
        print(json.dumps(genes_by_id(args.db, args.ids), sort_keys=True, separators=(",", ":")))
    except Exception as error:
        print("Annotation error: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

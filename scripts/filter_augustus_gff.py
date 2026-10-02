#!/usr/bin/env python3
"""Filter AUGUSTUS GFF3 gene predictions by gene probability and length.

A gene is kept if its probability score is at least --min-score,
or if it is at least --min-length bp long (long genes are kept regardless of score).
"""

import argparse


def gene_passes(gene_line, min_score, min_length):
    """Return True if a gene line meets the score-or-length rule."""
    fields = gene_line.rstrip("\n").split("\t")
    start, end = int(fields[3]), int(fields[4])
    score = float(fields[5])
    length = end - start + 1
    return score >= min_score or length >= min_length


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input", help="AUGUSTUS GFF3 file")
parser.add_argument("output", help="Filtered GFF3 file")
parser.add_argument("--min-score", type=float, default=0.5,
                    help="Keep genes with at least this probability (default: 0.5)")
parser.add_argument("--min-length", type=int, default=300,
                    help="Keep genes at least this long in bp, whatever their score (default: 300)")
args = parser.parse_args()

kept = 0
removed = 0
block = None
gene_line = None

with open(args.input) as infile, open(args.output, "w") as outfile:
    for line in infile:
        if line.startswith("# start gene"):
            block = [line]
            gene_line = None
        elif block is not None:
            block.append(line)
            fields = line.split("\t")
            if len(fields) > 2 and fields[2] == "gene":
                gene_line = line
            if line.startswith("# end gene"):
                if gene_line is not None and gene_passes(gene_line, args.min_score, args.min_length):
                    outfile.writelines(block)
                    kept += 1
                else:
                    removed += 1
                block = None
        else:
            outfile.write(line)

print(f"Rule: keep genes with score >= {args.min_score} or length >= {args.min_length} bp")
print(f"Genes kept:    {kept:,}")
print(f"Genes removed: {removed:,}")

#!/usr/bin/env python3
"""Summarise antiSMASH region GenBank files into a table and print summary counts."""

import argparse
import glob
import os
from collections import Counter

from Bio import SeqIO


def read_region(path):
    """Return a dict describing the region feature in one antiSMASH region GenBank file."""
    record = SeqIO.read(path, "genbank")
    for feature in record.features:
        if feature.type == "region":
            products = feature.qualifiers.get("product", [])
            edge = feature.qualifiers.get("contig_edge", ["False"])[0] == "True"
            length = int(feature.location.end) - int(feature.location.start)
            name = os.path.basename(path).replace(".gbk", "")
            scaffold, region = name.split(".")
            scaffold_number = int(scaffold.split("_")[1])
            region_number = int(region.replace("region", ""))
            return {
                "label": f"{scaffold_number}.{region_number}",
                "scaffold_number": scaffold_number,
                "region_number": region_number,
                "length_bp": length,
                "products": "+".join(products),
                "contig_edge": edge,
            }
    return None


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("results_dir", help="antiSMASH results folder")
parser.add_argument("output", help="Output table (tab-separated)")
args = parser.parse_args()

paths = glob.glob(os.path.join(args.results_dir, "*.region*.gbk"))
regions = [read_region(path) for path in paths]
regions = [r for r in regions if r is not None]
regions.sort(key=lambda r: (r["scaffold_number"], r["region_number"]))

columns = ["label", "length_bp", "products", "contig_edge"]
with open(args.output, "w") as out:
    out.write("\t".join(columns) + "\n")
    for r in regions:
        out.write("\t".join(str(r[c]) for c in columns) + "\n")

total = len(regions)
on_edge = sum(r["contig_edge"] for r in regions)
print(f"Regions:             {total}")
print(f"On a scaffold edge:  {on_edge} ({on_edge / total:.0%})")
print(f"Not on an edge:      {total - on_edge}")
print()
print("Regions by type:")
for product, count in Counter(r["products"] for r in regions).most_common():
    print(f"  {count:3d}  {product}")

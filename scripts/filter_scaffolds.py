#!/usr/bin/env python3
"""Filter a FASTA file by sequence length and report summary statistics."""

import argparse


def read_fasta(path):
    """Yield (header, sequence) pairs from a FASTA file, one record at a time."""
    header = None
    seq_parts = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def write_fasta(handle, header, seq, width=60):
    """Write one FASTA record, wrapping the sequence at `width` characters per line."""
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def n50(lengths):
    """Return the N50: the length L such that sequences >= L contain half the total bases."""
    half = sum(lengths) / 2
    running_total = 0
    for length in sorted(lengths, reverse=True):
        running_total += length
        if running_total >= half:
            return length
    return 0


def summarise(label, lengths):
    """Print the number, total length and N50 of a set of sequences."""
    print(label)
    print(f"  Sequences:    {len(lengths):,}")
    print(f"  Total length: {sum(lengths):,} bp")
    print(f"  N50:          {n50(lengths):,} bp")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input", help="Input FASTA file")
parser.add_argument("output", help="Output FASTA file for sequences that pass the filter")
parser.add_argument("-m", "--min-length", type=int, default=500,
                    help="Minimum sequence length to keep (default: 500)")
args = parser.parse_args()

all_lengths = []
kept_lengths = []

with open(args.output, "w") as out:
    for header, seq in read_fasta(args.input):
        all_lengths.append(len(seq))
        if len(seq) >= args.min_length:
            kept_lengths.append(len(seq))
            write_fasta(out, header, seq)

summarise("Before filtering", all_lengths)
summarise(f"After filtering (>= {args.min_length} bp)", kept_lengths)

removed_count = len(all_lengths) - len(kept_lengths)
removed_bases = sum(all_lengths) - sum(kept_lengths)
print(f"Removed {removed_count:,} sequences ({removed_bases:,} bp)")

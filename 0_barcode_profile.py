#!/usr/bin/env python3

import argparse
import sys
from collections import Counter, defaultdict


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Analyze barcode composition and structure.")

    parser.add_argument("--input", required=True,
                        help="Input file containing barcode sequences.")
    parser.add_argument("--output", required=True,
                        help="Output file with barcode analysis results.")

    parser.add_argument("--no-header", action="store_true",
                        help="Indicate that the input file has no header row.")

    parser.add_argument("--barcode-column", type=int, default=1,
                        help="Column containing barcode sequences (1-based index). Default = 1.")

    return parser.parse_args()


def read_barcodes(input_file, barcode_column, no_header):
    """
    Read barcode sequences from the specified column of the input file.

    Returns:
        list of barcode sequences.
    """

    print("[INFO] Reading barcode sequences from input file...")

    barcodes = []
    column_index = barcode_column - 1

    with open(input_file, "r") as infile:

        if not no_header:
            next(infile, None)  # skip header row

        for line in infile:

            line = line.strip()

            if not line:
                continue

            fields = line.split("\t")

            # handle single-column files
            if len(fields) == 1:
                barcode = fields[0]
            else:
                if column_index >= len(fields):
                    print(f"[ERROR] Barcode column {barcode_column} does not exist in input.")
                    sys.exit(1)

                barcode = fields[column_index]

            barcodes.append(barcode.upper())

    print(f"[INFO] Loaded {len(barcodes)} barcodes.")

    return barcodes


def analyze_barcodes(barcodes):
    """
    Perform structural analysis of barcode sequences.
    """

    print("[INFO] Computing barcode statistics...")

    length_counts = Counter()
    base_counts = Counter()

    position_counts = defaultdict(lambda: Counter())

    for barcode in barcodes:

        length = len(barcode)
        length_counts[length] += 1

        for i, base in enumerate(barcode):

            base_counts[base] += 1
            position_counts[i][base] += 1

    return length_counts, base_counts, position_counts


def infer_structure(position_counts, total_barcodes):
    """
    Infer a naive structural motif based on positional variability.
    """

    motif = []

    for pos in sorted(position_counts):

        counts = position_counts[pos]

        bases = [b for b in "ACGT" if counts[b] > 0]

        if len(bases) == 1:
            motif.append(bases[0])
        elif len(bases) == 4:
            motif.append("N")
        else:
            motif.append("[" + "".join(sorted(bases)) + "]")

    return "".join(motif)


def write_output(output_file, barcodes, length_counts, base_counts, position_counts, motif):
    """
    Write analysis results to output file.
    """

    print("[INFO] Writing results to output file...")

    total_bases = sum(base_counts.values())

    with open(output_file, "w") as out:

        out.write("# Barcode analysis summary\n")

        out.write(f"Total_barcodes\t{len(barcodes)}\n")
        out.write(f"Total_bases\t{total_bases}\n")

        out.write("\n# Length_distribution\n")
        out.write("Length\tCount\n")

        for length in sorted(length_counts):
            out.write(f"{length}\t{length_counts[length]}\n")

        out.write("\n# Base_composition\n")
        out.write("Base\tCount\tFrequency\n")

        for base in "ACGT":
            count = base_counts[base]
            freq = count / total_bases if total_bases > 0 else 0
            out.write(f"{base}\t{count}\t{freq:.6f}\n")

        out.write("\n# Position_base_frequencies\n")
        out.write("Position\tA\tC\tG\tT\n")

        for pos in sorted(position_counts):

            counts = position_counts[pos]
            total = sum(counts[b] for b in "ACGT")

            freqs = []
            for base in "ACGT":
                if total > 0:
                    freqs.append(counts[base] / total)
                else:
                    freqs.append(0)

            out.write(
                f"{pos+1}\t"
                f"{freqs[0]:.6f}\t"
                f"{freqs[1]:.6f}\t"
                f"{freqs[2]:.6f}\t"
                f"{freqs[3]:.6f}\n"
            )

        out.write("\n# Inferred_structure\n")
        out.write(f"Motif\t{motif}\n")


def main():

    args = parse_arguments()

    print("[INFO] Starting barcode profile analysis...")

    barcodes = read_barcodes(args.input, args.barcode_column, args.no_header)

    if len(barcodes) == 0:
        print("[ERROR] No barcode sequences were found.")
        sys.exit(1)

    length_counts, base_counts, position_counts = analyze_barcodes(barcodes)

    motif = infer_structure(position_counts, len(barcodes))

    write_output(args.output, barcodes, length_counts, base_counts, position_counts, motif)

    print("[INFO] Analysis complete.")


if __name__ == "__main__":
    main()
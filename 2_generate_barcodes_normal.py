#!/usr/bin/env python3

import argparse
import json
import random


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Generate barcode FASTQ reads with abundances drawn from a normal distribution."
    )

    parser.add_argument("--input", required=True,
                        help="Input barcode model JSON file")

    parser.add_argument("--output", required=True,
                        help="Output FASTQ file")

    parser.add_argument("--n-barcodes-main", type=int, required=True,
                        help="Number of unique barcodes")

    parser.add_argument("--n-barcodes", type=int, required=True,
                        help="Total number of FASTQ reads")

    parser.add_argument("--mean", type=float, default=100,
                        help="Mean of normal distribution (default=100)")

    parser.add_argument("--sd", type=float, default=20,
                        help="Standard deviation of normal distribution (default=20)")

    parser.add_argument("--max-homopolymer", type=int, default=5,
                        help="Maximum allowed homopolymer length (default=5, 0 disables check)")

    return parser.parse_args()


def load_model(model_file):
    """
    Load barcode model JSON.
    """

    print("[INFO] Loading barcode model...")

    with open(model_file) as f:
        model = json.load(f)

    return model


def sample_base(freq_dict):
    """
    Sample base according to positional probability distribution.
    """

    bases = list(freq_dict.keys())
    weights = list(freq_dict.values())

    return random.choices(bases, weights=weights, k=1)[0]


def generate_barcode(model):
    """
    Generate barcode sequence from positional model.
    """

    seq = []

    for pos in model["positions"]:
        seq.append(sample_base(pos))

    return "".join(seq)


def check_homopolymer(seq, max_run):
    """
    Enforce homopolymer constraint.
    """

    if max_run == 0:
        return True

    run = 1

    for i in range(1, len(seq)):

        if seq[i] == seq[i - 1]:
            run += 1
            if run > max_run:
                return False
        else:
            run = 1

    return True


def simulate_quality(length):
    """
    Simulate Illumina-like quality scores (Q30-Q40).
    """

    return "".join(chr(random.randint(30, 40) + 33) for _ in range(length))


def generate_main_barcodes(model, n_main, max_homopolymer):
    """
    Generate unique main barcode sequences.
    """

    print(f"[INFO] Generating {n_main} unique barcodes...")

    barcodes = set()

    attempts = 0

    while len(barcodes) < n_main:

        attempts += 1

        seq = generate_barcode(model)

        if not check_homopolymer(seq, max_homopolymer):
            continue

        barcodes.add(seq)

    print(f"[INFO] Barcode generation attempts: {attempts}")

    return list(barcodes)


def sample_abundances(n_main, mean, sd):
    """
    Sample barcode abundances from a normal distribution.
    """

    abundances = []

    for _ in range(n_main):

        val = int(random.gauss(mean, sd))

        if val < 1:
            val = 1

        abundances.append(val)

    return abundances


def scale_abundances(abundances, total_target):
    """
    Scale abundance list so the sum equals total_target.
    """

    total_current = sum(abundances)

    scale = total_target / total_current

    scaled = [max(1, int(a * scale)) for a in abundances]

    diff = total_target - sum(scaled)

    i = 0

    while diff != 0:

        if diff > 0:
            scaled[i] += 1
            diff -= 1
        else:
            if scaled[i] > 1:
                scaled[i] -= 1
                diff += 1

        i = (i + 1) % len(scaled)

    return scaled


def main():

    args = parse_arguments()

    print("[INFO] Starting FASTQ generation with normal abundance distribution...")

    model = load_model(args.input)

    n_main = args.n_barcodes_main
    total_reads = args.n_barcodes

    mean = args.mean
    sd = args.sd

    max_homopolymer = args.max_homopolymer

    length = model["length"]

    main_barcodes = generate_main_barcodes(model, n_main, max_homopolymer)

    print("[INFO] Sampling barcode abundances...")

    abundances = sample_abundances(n_main, mean, sd)

    print("[INFO] Scaling abundances to match total reads...")

    abundances = scale_abundances(abundances, total_reads)

    print("[INFO] Writing FASTQ output...")

    read_id = 1

    with open(args.output, "w") as out:

        for barcode, count in zip(main_barcodes, abundances):

            for _ in range(count):

                qual = simulate_quality(length)

                out.write(f"@barcode_{read_id}\n")
                out.write(barcode + "\n")
                out.write("+\n")
                out.write(qual + "\n")

                read_id += 1

    print("[INFO] FASTQ generation complete.")
    print(f"[INFO] Output reads: {read_id - 1}")
    print(f"[INFO] Output file: {args.output}")


if __name__ == "__main__":
    main()
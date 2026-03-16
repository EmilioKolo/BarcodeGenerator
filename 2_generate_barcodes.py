#!/usr/bin/env python3

import argparse
import json
import random


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(description="Generate synthetic barcode FASTQ reads from a probabilistic model.")

    parser.add_argument("--input", required=True,
                        help="Input barcode model JSON file")

    parser.add_argument("--output", required=True,
                        help="Output FASTQ file")

    parser.add_argument("--n-barcodes", type=int, required=True,
                        help="Number of synthetic barcodes to generate")

    parser.add_argument("--max-homopolymer", type=int, default=5,
                        help="Maximum allowed homopolymer length (default = 5, 0 disables check)")

    return parser.parse_args()


def load_model(model_file):
    """
    Load barcode generation model.
    """

    print("[INFO] Loading barcode model...")

    with open(model_file) as f:
        model = json.load(f)

    return model


def sample_base(freq_dict):
    """
    Sample a base according to its probability distribution.
    """

    bases = list(freq_dict.keys())
    weights = list(freq_dict.values())

    return random.choices(bases, weights=weights, k=1)[0]


def generate_barcode(model):
    """
    Generate a single barcode sequence using the positional probability model.
    """

    sequence = []

    for pos in model["positions"]:
        base = sample_base(pos)
        sequence.append(base)

    return "".join(sequence)


def check_homopolymer(seq, max_run):
    """
    Check if sequence contains homopolymers longer than allowed.
    """

    if max_run == 0:
        return True

    run = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            run += 1
            if run > max_run:
                return False
        else:
            run = 1

    return True


def simulate_quality(length):
    """
    Generate a FASTQ quality string using simulated Phred scores (Q30-Q40).
    """

    qualities = []

    for _ in range(length):

        q = random.randint(30, 40)

        qualities.append(chr(q + 33))

    return "".join(qualities)


def main():

    args = parse_arguments()

    print("[INFO] Starting barcode generation...")

    model = load_model(args.input)

    length = model["length"]

    n_target = args.n_barcodes
    max_homopolymer = args.max_homopolymer

    print(f"[INFO] Target barcodes: {n_target}")
    print(f"[INFO] Barcode length: {length}")

    generated = 0
    attempts = 0

    with open(args.output, "w") as out:

        while generated < n_target:

            attempts += 1

            seq = generate_barcode(model)

            if not check_homopolymer(seq, max_homopolymer):
                continue

            qual = simulate_quality(length)

            read_id = f"@barcode_{generated+1}"

            out.write(read_id + "\n")
            out.write(seq + "\n")
            out.write("+\n")
            out.write(qual + "\n")

            generated += 1

            if generated % 10000 == 0:
                print(f"[INFO] Generated {generated} barcodes...")

    print("[INFO] Barcode generation complete.")
    print(f"[INFO] Total attempts: {attempts}")
    print(f"[INFO] Output file: {args.output}")


if __name__ == "__main__":
    main()
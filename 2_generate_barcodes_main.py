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
                        help="Total number of FASTQ reads to generate")

    parser.add_argument("--n-barcodes-main", type=int, required=True,
                        help="Number of unique 'main' barcodes")

    parser.add_argument("--main-fraction", type=float, default=0.99,
                        help="Fraction of reads that are exact copies of main barcodes (default 0.99)")

    parser.add_argument("--modeled-error", action="store_true",
                        help="Restrict mutation to bases allowed by the positional model")

    parser.add_argument("--max-homopolymer", type=int, default=5,
                        help="Maximum allowed homopolymer length (default=5, 0 disables check)")

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
    Sample a base according to positional probability distribution.
    """

    bases = list(freq_dict.keys())
    weights = list(freq_dict.values())

    return random.choices(bases, weights=weights, k=1)[0]


def generate_barcode(model):
    """
    Generate a single barcode sequence.
    """

    seq = []

    for pos in model["positions"]:
        seq.append(sample_base(pos))

    return "".join(seq)


def check_homopolymer(seq, max_run):
    """
    Check homopolymer constraint.
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
    Simulate Illumina-like quality scores (Q30-Q40).
    """

    return "".join(chr(random.randint(30, 40) + 33) for _ in range(length))


def mutate_barcode(seq, model, modeled_error):
    """
    Apply a single random mutation to a barcode.
    """

    pos = random.randrange(len(seq))
    original = seq[pos]

    if modeled_error:

        allowed = list(model["positions"][pos].keys())

        # remove original base
        allowed = [b for b in allowed if b != original]

        if not allowed:
            return seq  # no valid mutation possible

        new_base = random.choice(allowed)

    else:

        bases = ["A", "C", "G", "T"]
        bases.remove(original)

        new_base = random.choice(bases)

    mutated = list(seq)
    mutated[pos] = new_base

    return "".join(mutated)


def generate_main_barcodes(model, n_main, max_homopolymer):
    """
    Generate unique main barcode set.
    """

    print(f"[INFO] Generating {n_main} unique main barcodes...")

    main_barcodes = set()

    attempts = 0

    while len(main_barcodes) < n_main:

        attempts += 1

        seq = generate_barcode(model)

        if not check_homopolymer(seq, max_homopolymer):
            continue

        main_barcodes.add(seq)

    print(f"[INFO] Main barcode generation attempts: {attempts}")

    return list(main_barcodes)


def main():

    args = parse_arguments()

    print("[INFO] Starting barcode generation...")

    model = load_model(args.input)

    length = model["length"]

    total_reads = args.n_barcodes
    n_main = args.n_barcodes_main
    main_fraction = args.main_fraction
    modeled_error = args.modeled_error
    max_homopolymer = args.max_homopolymer

    print(f"[INFO] Total reads: {total_reads}")
    print(f"[INFO] Main barcode count: {n_main}")
    print(f"[INFO] Main fraction: {main_fraction}")

    main_barcodes = generate_main_barcodes(model, n_main, max_homopolymer)

    generated = 0

    with open(args.output, "w") as out:

        while generated < total_reads:

            # choose a main barcode
            base_seq = random.choice(main_barcodes)

            if random.random() <= main_fraction:

                seq = base_seq

            else:

                seq = mutate_barcode(base_seq, model, modeled_error)

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
                print(f"[INFO] Generated {generated} reads...")

    print("[INFO] FASTQ generation complete.")
    print(f"[INFO] Output file: {args.output}")


if __name__ == "__main__":
    main()
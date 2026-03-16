#!/usr/bin/env python3

import argparse
import json
from collections import defaultdict


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(description="Convert barcode profile into a simplified generation model.")

    parser.add_argument("--input", required=True,
                        help="Input profile file produced by barcode_profile.py")

    parser.add_argument("--output", required=True,
                        help="Output JSON model file")

    return parser.parse_args()


def round_to_step(value, step=0.05):
    """
    Round value to nearest step (default 0.05).
    """

    return round(value / step) * step


def parse_profile(profile_file):
    """
    Extract length distribution and positional frequencies
    from barcode_profile output.
    """

    print("[INFO] Reading barcode profile...")

    length_counts = {}
    position_freqs = {}

    section = None

    with open(profile_file) as f:
        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                if "Length_distribution" in line:
                    section = "length"
                elif "Position_base_frequencies" in line:
                    section = "positions"
                else:
                    section = None
                continue

            if section == "length":

                if line.startswith("Length"):
                    continue

                length, count = line.split("\t")
                length_counts[int(length)] = int(count)

            elif section == "positions":

                if line.startswith("Position"):
                    continue

                fields = line.split("\t")

                pos = int(fields[0])

                position_freqs[pos] = {
                    "A": float(fields[1]),
                    "C": float(fields[2]),
                    "G": float(fields[3]),
                    "T": float(fields[4])
                }

    return length_counts, position_freqs


def determine_length(length_counts):
    """
    Select the most probable barcode length.
    """

    total = sum(length_counts.values())

    best_length = max(length_counts, key=lambda x: length_counts[x])
    best_count = length_counts[best_length]

    proportion = best_count / total

    if proportion < 0.95:
        print(f"[WARNING] Dominant barcode length ({best_length}) represents only {proportion:.3f} of sequences.")

    print(f"[INFO] Selected barcode length: {best_length}")

    return best_length


def simplify_position(freq_dict):
    """
    Simplify positional frequencies by removing small proportions
    and rounding to 0.05 increments.
    """

    bases = ["A", "C", "G", "T"]

    sorted_freqs = sorted(freq_dict.values(), reverse=True)

    # determine threshold rule
    if sorted_freqs[0] + sorted_freqs[1] > 0.80:
        threshold = 0.10
    else:
        threshold = 0.05

    cleaned = {}

    for base in bases:

        value = freq_dict[base]

        if value < threshold:
            continue

        value = round_to_step(value)

        if value > 0:
            cleaned[base] = value

    # renormalize
    total = sum(cleaned.values())

    if total == 0:
        # fallback uniform distribution
        cleaned = {b: 0.25 for b in bases}
        total = 1.0

    for base in cleaned:
        cleaned[base] = round(cleaned[base] / total, 2)

    return cleaned


def build_model(length, position_freqs):
    """
    Build simplified positional model.
    """

    print("[INFO] Simplifying positional distributions...")

    positions = []

    for pos in sorted(position_freqs):

        simplified = simplify_position(position_freqs[pos])

        positions.append(simplified)

    return {
        "length": length,
        "positions": positions
    }


def write_model(output_file, model):
    """
    Write model to JSON file.
    """

    print("[INFO] Writing model JSON...")

    with open(output_file, "w") as f:
        json.dump(model, f, indent=2)


def main():

    args = parse_arguments()

    print("[INFO] Starting barcode model extraction...")

    length_counts, position_freqs = parse_profile(args.input)

    length = determine_length(length_counts)

    model = build_model(length, position_freqs)

    write_model(args.output, model)

    print("[INFO] Model generation complete.")


if __name__ == "__main__":
    main()
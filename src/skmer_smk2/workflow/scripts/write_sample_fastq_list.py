#!/usr/bin/env python3
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Write two-column sample FASTQ list.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--fastqs", nargs="+", required=True)
    args = parser.parse_args()

    if len(args.samples) != len(args.fastqs):
        raise ValueError("--samples and --fastqs must have the same length")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as out:
        for sample, fastq in zip(args.samples, args.fastqs):
            out.write("{}\t{}\n".format(sample, os.path.abspath(fastq)))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import os


def norm_name(value):
    name = os.path.basename(value.strip())
    for ext in (".fq.gz", ".fastq.gz", ".fq", ".fastq", ".msh"):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name


def main():
    parser = argparse.ArgumentParser(description="Convert mash dist output to relaxed PHYLIP.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    args = parser.parse_args()

    matrix = {a: {b: 0.0 for b in args.samples} for a in args.samples}
    sample_set = set(args.samples)
    with open(args.input, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            a = norm_name(parts[0])
            b = norm_name(parts[1])
            if a in sample_set and b in sample_set:
                matrix[a][b] = float(parts[2])
                matrix[b][a] = float(parts[2])

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write("{}\n".format(len(args.samples)))
        for sample in args.samples:
            values = " ".join("{:.10f}".format(matrix[sample][other]) for other in args.samples)
            out.write("{} {}\n".format(sample, values))


if __name__ == "__main__":
    main()

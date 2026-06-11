#!/usr/bin/env python3
import argparse
import csv
import os
import re


def norm_name(value):
    name = os.path.basename(value.strip())
    for ext in (".fq.gz", ".fastq.gz", ".fq", ".fastq", ".msh"):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name


def read_table(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    return [re.split(r"\s+", line) for line in lines]


def parse_matrix(rows, samples):
    matrix = {a: {b: 0.0 for b in samples} for a in samples}
    sample_set = set(samples)

    if not rows:
        return matrix

    header = [norm_name(x) for x in rows[0]]
    if len(header) >= 2 and all(x in sample_set for x in header[1:]):
        for row in rows[1:]:
            if len(row) < 2:
                continue
            src = norm_name(row[0])
            if src not in sample_set:
                continue
            for dst, value in zip(header[1:], row[1:]):
                if dst in sample_set:
                    matrix[src][dst] = float(value)
        return mirror(matrix, samples)

    if len(rows) == len(samples) and all(len(row) >= len(samples) for row in rows):
        for src, row in zip(samples, rows):
            values = row[-len(samples):]
            for dst, value in zip(samples, values):
                matrix[src][dst] = float(value)
        return mirror(matrix, samples)

    for row in rows:
        if len(row) < 3:
            continue
        a = norm_name(row[0])
        b = norm_name(row[1])
        if a in sample_set and b in sample_set:
            matrix[a][b] = float(row[2])
            matrix[b][a] = float(row[2])
    return matrix


def mirror(matrix, samples):
    for a in samples:
        for b in samples:
            if a == b:
                matrix[a][b] = 0.0
            elif matrix[a][b] == 0.0 and matrix[b][a] != 0.0:
                matrix[a][b] = matrix[b][a]
            elif matrix[b][a] == 0.0 and matrix[a][b] != 0.0:
                matrix[b][a] = matrix[a][b]
    return matrix


def write_phylip(path, samples, matrix):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as out:
        out.write("{}\n".format(len(samples)))
        for sample in samples:
            distances = " ".join("{:.10f}".format(matrix[sample][other]) for other in samples)
            out.write("{} {}\n".format(sample, distances))


def main():
    parser = argparse.ArgumentParser(description="Convert a distance table or matrix to relaxed PHYLIP.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    args = parser.parse_args()

    rows = read_table(args.input)
    matrix = parse_matrix(rows, args.samples)
    write_phylip(args.output, args.samples, matrix)


if __name__ == "__main__":
    main()

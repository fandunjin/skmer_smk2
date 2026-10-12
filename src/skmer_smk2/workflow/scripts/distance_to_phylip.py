#!/usr/bin/env python3
import argparse
import csv
import math
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
        lines = [line.rstrip("\r\n") for line in handle if line.strip()]
    return [line.split("\t") if "\t" in line else line.split() for line in lines]


def numeric_distance(value):
    text = str(value).strip()
    wrapped = re.fullmatch(r"(?:np|numpy)\.float(?:16|32|64|128)\(([^()]*)\)", text)
    if wrapped:
        text = wrapped.group(1)
    try:
        number = float(text)
    except ValueError:
        raise ValueError("Invalid or missing distance: {!r}".format(value)) from None
    if not math.isfinite(number) or number < 0:
        raise ValueError("Distance must be finite and nonnegative: {!r}".format(value))
    return number


def parse_matrix(rows, samples):
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Sample names must be nonempty and unique")
    matrix = {a: {b: None for b in samples} for a in samples}
    sample_set = set(samples)

    if not rows:
        raise ValueError("Empty distance matrix")

    header = [norm_name(x) for x in rows[0]]
    if len(header) == len(samples) + 1 and set(header[1:]) == sample_set:
        if len(set(header[1:])) != len(samples):
            raise ValueError("Duplicate matrix columns")
        seen = set()
        for row in rows[1:]:
            if len(row) != len(header):
                raise ValueError("Wrong matrix row width: {!r}".format(row))
            src = norm_name(row[0])
            if src not in sample_set or src in seen:
                raise ValueError("Unknown or duplicate matrix row: " + src)
            seen.add(src)
            for dst, value in zip(header[1:], row[1:]):
                if dst in sample_set:
                    matrix[src][dst] = numeric_distance(value)
        if seen != sample_set:
            raise ValueError("Missing matrix rows: {}".format(sorted(sample_set - seen)))
        return mirror(matrix, samples)

    if len(rows) == len(samples) and all(len(row) in (len(samples), len(samples)+1) for row in rows):
        seen = set()
        for src, row in zip(samples, rows):
            if len(row) == len(samples)+1:
                src = norm_name(row[0])
            if src not in sample_set or src in seen:
                raise ValueError("Unknown or duplicate matrix row: " + src)
            seen.add(src)
            values = row[-len(samples):]
            for dst, value in zip(samples, values):
                matrix[src][dst] = numeric_distance(value)
        return mirror(matrix, samples)

    for row in rows:
        if len(row) != 3:
            raise ValueError("Unrecognized distance table or incomplete matrix")
        a = norm_name(row[0])
        b = norm_name(row[1])
        if a in sample_set and b in sample_set:
            number = numeric_distance(row[2])
            if matrix[a][b] is not None and not math.isclose(matrix[a][b], number, rel_tol=1e-6, abs_tol=1e-10):
                raise ValueError("Conflicting pairwise distances: {}, {}".format(a, b))
            matrix[a][b] = number
            matrix[b][a] = number
        else:
            raise ValueError("Unknown sample in distance table: {}, {}".format(a, b))
    return mirror(matrix, samples)


def mirror(matrix, samples):
    for a in samples:
        for b in samples:
            if a == b:
                if matrix[a][b] not in (None, 0.0):
                    raise ValueError("Nonzero diagonal: " + a)
                matrix[a][b] = 0.0
            elif matrix[a][b] is None:
                raise ValueError("Missing distance: {}, {}".format(a, b))
            elif not math.isclose(matrix[a][b], matrix[b][a], rel_tol=1e-6, abs_tol=1e-10):
                raise ValueError("Asymmetric distances: {}, {}".format(a, b))
    if len(samples) > 1 and not any(matrix[a][b] > 0 for a in samples for b in samples):
        raise ValueError("All distances are zero; cannot infer an informative tree")
    return matrix


def write_phylip(path, samples, matrix):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
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

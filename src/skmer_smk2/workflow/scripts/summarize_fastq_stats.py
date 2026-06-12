#!/usr/bin/env python3
import argparse
import math
import os


def read_stats(path):
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        values = handle.readline().rstrip("\n").split("\t")
    row = dict(zip(header, values))
    return {
        "sample": row["sample"],
        "reads": int(row["reads"]),
        "total_bases": int(row["total_bases"]),
        "average_length": float(row["average_length"]),
    }


def main():
    parser = argparse.ArgumentParser(description="Merge FASTQ stats and optionally choose a base cutoff.")
    parser.add_argument("stats", nargs="+")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--cutoff")
    parser.add_argument("--percentile", type=float, default=75.0)
    args = parser.parse_args()

    rows = [read_stats(path) for path in args.stats]
    rows.sort(key=lambda row: (-row["total_bases"], row["sample"]))

    os.makedirs(os.path.dirname(args.summary), exist_ok=True)
    with open(args.summary, "w", encoding="utf-8") as out:
        out.write("sample\treads\ttotal_bases\taverage_length\n")
        for row in rows:
            out.write(
                "{sample}\t{reads}\t{total_bases}\t{average_length:.6f}\n".format(**row)
            )

    if args.cutoff:
        pct = max(0.0, min(100.0, args.percentile))
        if not rows:
            cutoff = 0
            index = 0
        else:
            index = max(0, min(len(rows) - 1, int(math.ceil(len(rows) * pct / 100.0) - 1)))
            cutoff = rows[index]["total_bases"]
        os.makedirs(os.path.dirname(args.cutoff), exist_ok=True)
        with open(args.cutoff, "w", encoding="utf-8") as out:
            out.write("{}\n".format(cutoff))
            out.write("# percentile={}\n".format(args.percentile))
            out.write("# sorted_position={}\n".format(index + 1))


if __name__ == "__main__":
    main()

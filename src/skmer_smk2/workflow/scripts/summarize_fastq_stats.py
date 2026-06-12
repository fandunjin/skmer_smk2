#!/usr/bin/env python3
import argparse
import math
import os


def read_stats(path):
    with open(path, "r", encoding="utf-8-sig") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        values = handle.readline().rstrip("\n").split("\t")
    row = dict(zip(header, values))
    return {
        "sample": row["sample"],
        "reads": int(row["reads"]),
        "total_bases": int(row["total_bases"]),
        "average_length": float(row["average_length"]),
    }


def percentile_index(row_count, percentile):
    pct = max(0.0, min(100.0, percentile))
    if row_count <= 0:
        return 0
    return max(0, min(row_count - 1, int(math.ceil(row_count * pct / 100.0) - 1)))


def parse_percentiles(value):
    percentiles = []
    for item in str(value).replace(",", " ").split():
        try:
            percentiles.append(float(item))
        except ValueError:
            raise ValueError("Invalid percentile value: {}".format(item))
    return percentiles


def write_candidate_report(path, rows, percentiles, selected_percentile=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    total_available = sum(row["total_bases"] for row in rows)
    sample_count = len(rows)

    with open(path, "w", encoding="utf-8") as out:
        out.write(
            "percentile\tselected\tposition\tsample\tcutoff_bases\t"
            "samples_total\tsamples_truncated\tsamples_below_cutoff\t"
            "estimated_retained_bases\testimated_retained_percent\t"
            "largest_sample_bases\tsmallest_sample_bases\t"
            "smallest_sample_percent_of_cutoff\n"
        )
        if not rows:
            return

        largest = rows[0]["total_bases"]
        smallest = rows[-1]["total_bases"]
        seen = set()
        ordered_percentiles = []
        for pct in percentiles:
            key = "{:.6f}".format(max(0.0, min(100.0, pct)))
            if key not in seen:
                seen.add(key)
                ordered_percentiles.append(max(0.0, min(100.0, pct)))

        for pct in ordered_percentiles:
            index = percentile_index(sample_count, pct)
            cutoff = rows[index]["total_bases"]
            truncated = sum(1 for row in rows if row["total_bases"] > cutoff)
            below = sum(1 for row in rows if row["total_bases"] < cutoff)
            retained = sum(min(row["total_bases"], cutoff) for row in rows)
            retained_percent = (retained / total_available * 100.0) if total_available else 0.0
            smallest_percent = (smallest / cutoff * 100.0) if cutoff else 0.0
            selected = "yes" if selected_percentile is not None and abs(pct - selected_percentile) < 1e-9 else "no"
            out.write(
                "{percentile:g}\t{selected}\t{position}\t{sample}\t{cutoff}\t"
                "{sample_count}\t{truncated}\t{below}\t{retained}\t"
                "{retained_percent:.6f}\t{largest}\t{smallest}\t"
                "{smallest_percent:.6f}\n".format(
                    percentile=pct,
                    selected=selected,
                    position=index + 1,
                    sample=rows[index]["sample"],
                    cutoff=cutoff,
                    sample_count=sample_count,
                    truncated=truncated,
                    below=below,
                    retained=retained,
                    retained_percent=retained_percent,
                    largest=largest,
                    smallest=smallest,
                    smallest_percent=smallest_percent,
                )
            )


def main():
    parser = argparse.ArgumentParser(description="Merge FASTQ stats and optionally choose a base cutoff.")
    parser.add_argument("stats", nargs="+")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--cutoff")
    parser.add_argument("--percentile", type=float, default=75.0)
    parser.add_argument("--candidates")
    parser.add_argument("--candidate-percentiles", default="50,60,70,75,80,90,95")
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
        if not rows:
            cutoff = 0
            index = 0
        else:
            index = percentile_index(len(rows), args.percentile)
            cutoff = rows[index]["total_bases"]
        os.makedirs(os.path.dirname(args.cutoff), exist_ok=True)
        with open(args.cutoff, "w", encoding="utf-8") as out:
            out.write("{}\n".format(cutoff))
            out.write("# percentile={}\n".format(args.percentile))
            out.write("# sorted_position={}\n".format(index + 1))

    if args.candidates:
        candidate_percentiles = parse_percentiles(args.candidate_percentiles)
        if args.percentile not in candidate_percentiles:
            candidate_percentiles.append(args.percentile)
        write_candidate_report(args.candidates, rows, candidate_percentiles, args.percentile)


if __name__ == "__main__":
    main()

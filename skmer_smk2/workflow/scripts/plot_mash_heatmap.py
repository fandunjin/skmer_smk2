#!/usr/bin/env python3
import argparse
import html
import os


def norm_name(value):
    name = os.path.basename(value.strip())
    for ext in (".fq.gz", ".fastq.gz", ".fq", ".fastq", ".msh"):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name


def color(value, max_value):
    if max_value <= 0:
        ratio = 0.0
    else:
        ratio = max(0.0, min(1.0, value / max_value))
    r = int(245 - ratio * 202)
    g = int(247 - ratio * 145)
    b = int(250 - ratio * 64)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def main():
    parser = argparse.ArgumentParser(description="Write a simple SVG heatmap from mash dist output.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    args = parser.parse_args()

    matrix = {a: {b: 0.0 for b in args.samples} for a in args.samples}
    sample_set = set(args.samples)
    with open(args.input, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            a = norm_name(parts[0])
            b = norm_name(parts[1])
            if a in sample_set and b in sample_set:
                matrix[a][b] = float(parts[2])
                matrix[b][a] = float(parts[2])

    max_value = max((matrix[a][b] for a in args.samples for b in args.samples), default=0.0)
    cell = 28
    label = max(120, min(320, max(len(s) for s in args.samples) * 7 + 20))
    width = label + cell * len(args.samples) + 30
    height = label + cell * len(args.samples) + 30

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write('<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}">\n'.format(width, height))
        out.write('<rect width="100%" height="100%" fill="white"/>\n')
        out.write('<style>text{font-family:Arial,sans-serif;font-size:11px}.tick{fill:#222}.value{fill:#111;font-size:8px}</style>\n')
        for i, sample in enumerate(args.samples):
            x = label + i * cell + cell / 2
            out.write('<text class="tick" transform="translate({:.1f},{:.1f}) rotate(-60)" text-anchor="start">{}</text>\n'.format(x, label - 8, html.escape(sample)))
            y = label + i * cell + cell / 2 + 4
            out.write('<text class="tick" x="{}" y="{:.1f}" text-anchor="end">{}</text>\n'.format(label - 8, y, html.escape(sample)))
        for y_i, a in enumerate(args.samples):
            for x_i, b in enumerate(args.samples):
                value = matrix[a][b]
                x = label + x_i * cell
                y = label + y_i * cell
                out.write('<rect x="{}" y="{}" width="{}" height="{}" fill="{}" stroke="#e5e7eb"/>\n'.format(x, y, cell, cell, color(value, max_value)))
                if cell >= 24:
                    out.write('<text class="value" x="{:.1f}" y="{:.1f}" text-anchor="middle">{:.3f}</text>\n'.format(x + cell / 2, y + cell / 2 + 3, value))
        out.write("</svg>\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import gzip
import os


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def iter_fastq(path):
    with open_text(path, "rt") as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError("Incomplete FASTQ record in {}".format(path))
            yield name, seq, plus, qual


def write_stats(path, sample, reads, bases):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    avg = (bases / reads) if reads else 0.0
    with open(path, "w", encoding="utf-8") as out:
        out.write("sample\treads\ttotal_bases\taverage_length\n")
        out.write("{}\t{}\t{}\t{:.6f}\n".format(sample, reads, bases, avg))


def stats(args):
    reads = 0
    bases = 0
    for _name, seq, _plus, _qual in iter_fastq(args.input):
        reads += 1
        bases += len(seq.rstrip("\n\r"))
    write_stats(args.output, args.sample, reads, bases)


def sample(args):
    with open(args.cutoff, "r", encoding="utf-8") as handle:
        first = handle.readline().strip()
    target_bases = int(float(first.split()[0]))

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    reads = 0
    bases = 0
    with open(args.output, "w", encoding="utf-8") as out:
        for name, seq, plus, qual in iter_fastq(args.input):
            seq_len = len(seq.rstrip("\n\r"))
            if reads > 0 and bases >= target_bases:
                break
            out.write(name)
            out.write(seq)
            out.write(plus)
            out.write(qual)
            reads += 1
            bases += seq_len
            if bases >= target_bases:
                break
    write_stats(args.stats, args.sample, reads, bases)


def main():
    parser = argparse.ArgumentParser(description="FASTQ statistics and base-aware head sampling.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_stats = sub.add_parser("stats")
    p_stats.add_argument("--sample", required=True)
    p_stats.add_argument("--input", required=True)
    p_stats.add_argument("--output", required=True)
    p_stats.set_defaults(func=stats)

    p_sample = sub.add_parser("sample")
    p_sample.add_argument("--sample", required=True)
    p_sample.add_argument("--input", required=True)
    p_sample.add_argument("--cutoff", required=True)
    p_sample.add_argument("--output", required=True)
    p_sample.add_argument("--stats", required=True)
    p_sample.set_defaults(func=sample)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

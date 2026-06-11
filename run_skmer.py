#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Run the Skmer/WASTER/Mash Snakemake workflow.")
    parser.add_argument("-i", "--input", required=True, help="Directory containing paired FASTQ files.")
    parser.add_argument("-ref", "--ref", default="", help="Optional plastid/reference genome FASTA for bowtie2 filtering.")
    parser.add_argument("-s", "--sample-percentile", type=float, default=75.0, help="Sorted sample base-count percentile used as the head cutoff.")
    parser.add_argument("-j", "--jobs", default="1", help="Snakemake job count/cores.")
    parser.add_argument("--snakefile", default="snakefile", help="Workflow file path.")
    parser.add_argument("--dry-run", action="store_true", help="Run Snakemake in dry-run mode.")
    parser.add_argument("snakemake_args", nargs=argparse.REMAINDER, help="Extra arguments passed to snakemake after --.")
    args = parser.parse_args()

    snakemake = shutil.which("snakemake")
    if not snakemake:
        raise SystemExit("snakemake was not found in PATH")

    input_dir = os.path.abspath(args.input)
    ref_path = os.path.abspath(args.ref) if args.ref else ""
    cmd = [
        snakemake,
        "-s",
        args.snakefile,
        "--cores",
        str(args.jobs),
        "--config",
        "input_dir={}".format(input_dir),
        "sample_percentile={}".format(args.sample_percentile),
    ]
    if ref_path:
        cmd.extend(["ref={}".format(ref_path)])
    if args.dry_run:
        cmd.append("-n")
    if args.snakemake_args:
        extra = args.snakemake_args
        if extra and extra[0] == "--":
            extra = extra[1:]
        cmd.extend(extra)

    print("Running:", " ".join(cmd), flush=True)
    return subprocess.call(cmd)


if __name__ == "__main__":
    sys.exit(main())

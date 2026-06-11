#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
from importlib.resources import files
from pathlib import Path


PACKAGE = "skmer_smk2"
WORKFLOW_CACHE = ".skmer_smk2_workflow"


def resource_path(*parts):
    return files(PACKAGE).joinpath(*parts)


def copy_tree(src, dst):
    src = Path(src)
    dst = Path(dst)
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))


def materialize_workflow(workdir):
    workdir = Path(workdir).resolve()
    target = workdir / WORKFLOW_CACHE
    copy_tree(resource_path("workflow"), target)
    return target / "Snakefile"


def run(args):
    snakemake = shutil.which("snakemake")
    if not snakemake:
        snakemake_cmd = [sys.executable, "-m", "snakemake"]
    else:
        snakemake_cmd = [snakemake]

    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    snakefile = materialize_workflow(workdir)

    input_dir = Path(args.input).resolve()
    ref_path = Path(args.ref).resolve() if args.ref else None

    cmd = snakemake_cmd + [
        "-s",
        str(snakefile),
        "--directory",
        str(workdir),
        "--cores",
        str(args.jobs),
        "--rerun-incomplete",
        "--latency-wait",
        str(args.latency_wait),
        "--config",
        "input_dir={}".format(input_dir),
        "sample_percentile={}".format(args.sample_percentile),
        "rep_n={}".format(args.bootstraps),
        "exclude_samples={}".format(args.exclude_samples or ""),
    ]
    if ref_path:
        cmd.append("ref={}".format(ref_path))
    if args.dry_run:
        cmd.append("--dry-run")
    if args.printshellcmds:
        cmd.append("--printshellcmds")
    if args.snakemake_args:
        extra = args.snakemake_args
        if extra and extra[0] == "--":
            extra = extra[1:]
        cmd.extend(extra)

    print("Running:", " ".join(cmd), flush=True)
    return subprocess.call(cmd)


def init(args):
    outdir = Path(args.output).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if args.with_workflow:
        copy_tree(resource_path("workflow"), outdir / "workflow")
    shutil.copyfile(resource_path("templates", "scan_repair_fastq.sh"), outdir / "scan_repair_fastq.sh")
    shutil.copyfile(resource_path("templates", "jsub_submit.sh"), outdir / "jsub_submit.sh")
    print("Wrote templates to {}".format(outdir))
    return 0


def repair_fastq(args):
    script = Path(args.workdir).resolve() / "scan_repair_fastq.sh"
    shutil.copyfile(resource_path("templates", "scan_repair_fastq.sh"), script)
    if args.copy_only:
        print("Wrote {}".format(script))
        return 0
    cmd = ["bash", str(script), str(Path(args.input).resolve())]
    print("Running:", " ".join(cmd), flush=True)
    return subprocess.call(cmd)


def build_parser():
    parser = argparse.ArgumentParser(prog="skmer-smk2", description="Packaged Skmer/WASTER/Mash Snakemake workflow.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="Run the workflow with local Snakemake.")
    p_run.add_argument("-i", "--input", required=True, help="Directory containing paired FASTQ files.")
    p_run.add_argument("-ref", "--ref", default="", help="Optional plastid/reference genome FASTA for bowtie2 filtering.")
    p_run.add_argument("-s", "--sample-percentile", type=float, default=75.0, help="Sorted sample base-count percentile used as the head cutoff.")
    p_run.add_argument("-j", "--jobs", default="1", help="Snakemake job count/cores.")
    p_run.add_argument("-b", "--bootstraps", type=int, default=100, help="Bootstrap replicate count.")
    p_run.add_argument("--exclude-samples", default="", help="Comma- or whitespace-separated sample names to skip.")
    p_run.add_argument("--workdir", default=".", help="Run directory for results and workflow cache.")
    p_run.add_argument("--latency-wait", default="120", help="Snakemake latency wait seconds.")
    p_run.add_argument("--dry-run", action="store_true", help="Run Snakemake in dry-run mode.")
    p_run.add_argument("--printshellcmds", action="store_true", help="Print shell commands from Snakemake.")
    p_run.add_argument("snakemake_args", nargs=argparse.REMAINDER, help="Extra arguments passed to snakemake after --.")
    p_run.set_defaults(func=run)

    p_init = sub.add_parser("init", help="Write optional helper templates into a directory.")
    p_init.add_argument("-o", "--output", default="skmer_smk2_templates")
    p_init.add_argument("--with-workflow", action="store_true", help="Also copy the packaged Snakefile and helper scripts.")
    p_init.set_defaults(func=init)

    p_repair = sub.add_parser("repair-fastq", help="Run or copy the FASTQ scan/repair helper.")
    p_repair.add_argument("-i", "--input", default=".", help="FASTQ directory to scan.")
    p_repair.add_argument("--workdir", default=".", help="Directory where the helper shell script is written.")
    p_repair.add_argument("--copy-only", action="store_true", help="Only copy the helper script.")
    p_repair.set_defaults(func=repair_fastq)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

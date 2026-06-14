#!/usr/bin/env python3
import argparse
import glob
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from .input_qc import input_check_main, input_repair_main


PACKAGE = "skmer_smk2"
WORKFLOW_CACHE = ".skmer_smk2_workflow"
REQUIRED_TOOLS = [
    "snakemake",
    "python",
    "fastp",
    "bowtie2",
    "bowtie2-build",
    "repair.sh",
    "bbmerge.sh",
    "skmer",
    "fastme",
    "raxmlHPC",
    "waster",
    "mash",
    "seqkit",
    "gzip",
]
TOOL_CANDIDATES = {
    "snakemake": ("snakemake",),
    "fastp": ("fastp",),
    "bowtie2": ("bowtie2",),
    "bowtie2-build": ("bowtie2-build",),
    "repair.sh": ("repair.sh", "repair"),
    "bbmerge.sh": ("bbmerge.sh", "bbmerge"),
    "skmer": ("skmer",),
    "fastme": ("fastme", "FastME"),
    "raxmlHPC": (
        "raxmlHPC",
        "raxmlHPC-PTHREADS",
        "raxmlHPC-PTHREADS-SSE3",
        "raxmlHPC-PTHREADS-AVX",
        "raxmlHPC-SSE3",
        "raxmlHPC-AVX",
        "raxml",
    ),
    "waster": ("waster",),
    "mash": ("mash",),
    "seqkit": ("seqkit",),
    "gzip": ("gzip",),
}
WORKFLOW_TOOL_CONFIG = {
    "python": "tool_python",
    "fastp": "tool_fastp",
    "bowtie2": "tool_bowtie2",
    "bowtie2-build": "tool_bowtie2_build",
    "repair.sh": "tool_repair",
    "bbmerge.sh": "tool_bbmerge",
    "skmer": "tool_skmer",
    "fastme": "tool_fastme",
    "raxmlHPC": "tool_raxml",
    "waster": "tool_waster",
    "mash": "tool_mash",
}
CONDA_PACKAGES = {
    "snakemake": "snakemake",
    "fastp": "fastp",
    "bowtie2": "bowtie2",
    "bowtie2-build": "bowtie2",
    "repair.sh": "bbmap",
    "bbmerge.sh": "bbmap",
    "fastme": "fastme",
    "raxmlHPC": "raxml",
    "mash": "mash",
    "seqkit": "seqkit",
    "gzip": "gzip",
}
CORE_TOOLS = {"python"}
MANUAL_TOOLS = {
    "skmer": "Install Skmer from its upstream instructions or make sure the skmer executable is in PATH.",
    "waster": "Install WASTER from its upstream instructions or make sure the waster executable is in PATH.",
}


def resource_path(*parts):
    return Path(__file__).resolve().parent.joinpath(*parts)


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


def conda_prefixes():
    prefixes = []
    for value in (os.environ.get("CONDA_PREFIX"), sys.prefix):
        if value:
            path = Path(value)
            if path not in prefixes:
                prefixes.append(path)
    return prefixes


def find_in_conda_prefixes(tool):
    for prefix in conda_prefixes():
        for subdir in ("bin", "Scripts"):
            path = prefix / subdir / tool
            if path.exists():
                return str(path)
        for pattern in ("share/bbmap*/{}", "opt/bbmap*/{}"):
            matches = sorted(glob.glob(str(prefix / pattern.format(tool))))
            if matches:
                return matches[0]
    return ""


def find_tool(tool):
    if tool == "python":
        return sys.executable
    for candidate in TOOL_CANDIDATES.get(tool, (tool,)):
        path = shutil.which(candidate)
        if path:
            return path
        path = find_in_conda_prefixes(candidate)
        if path:
            return path
    if tool == "snakemake" and importlib.util.find_spec("snakemake"):
        return "{} -m snakemake".format(sys.executable)
    return ""


def snakemake_command():
    snakemake = shutil.which("snakemake")
    if snakemake:
        return [snakemake]
    if importlib.util.find_spec("snakemake"):
        return [sys.executable, "-m", "snakemake"]
    return []


def selected_analyses(args):
    selected = []
    if args.skmer:
        selected.append("skmer")
    if args.waster:
        selected.append("waster")
    if args.mash:
        selected.append("mash")
    return selected or ["skmer", "waster", "mash"]


def run(args):
    snakemake_cmd = snakemake_command()
    if not snakemake_cmd:
        print(
            "ERROR: Snakemake was not found in PATH and cannot be imported by this Python. "
            "Activate an environment with snakemake, or run `skmer-smk2 doctor --install`.",
            file=sys.stderr,
        )
        return 127

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
        "candidate_percentiles={}".format(args.candidate_percentiles),
        "rep_n={}".format(args.bootstraps),
        "analyses={}".format(",".join(selected_analyses(args))),
        "exclude_samples={}".format(args.exclude_samples or ""),
        "bbmerge_timeout={}".format(args.bbmerge_timeout),
    ]
    for tool, config_key in WORKFLOW_TOOL_CONFIG.items():
        path = find_tool(tool)
        if path:
            cmd.append("{}={}".format(config_key, path))
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


def choose_package_manager(manager):
    if manager == "auto":
        return shutil.which("mamba") or shutil.which("conda")
    return shutil.which(manager)


def package_manager_name(path):
    if not path:
        return ""
    return Path(path).name.lower().split(".")[0]


def doctor(args):
    missing = []
    width = max(len(tool) for tool in REQUIRED_TOOLS)

    print("Checking skmer_smk2 external tools")
    print("Mode: {}".format("strict" if args.strict else "advisory"))
    print()
    print("{:<{width}}  {:<7}  {}".format("tool", "status", "resolved command/path", width=width))
    print("{:<{width}}  {:<7}  {}".format("-" * width, "-------", "---------------------", width=width))

    for tool in REQUIRED_TOOLS:
        path = find_tool(tool)
        if path:
            print("{:<{width}}  {:<7}  {}".format(tool, "OK", path, width=width))
        else:
            aliases = ", ".join(TOOL_CANDIDATES.get(tool, (tool,)))
            print("{:<{width}}  {:<7}  {} aliases: {}".format(tool, "WARN", "-", aliases, width=width))
            missing.append(tool)

    print()
    if not missing:
        print("All checked tools were found.")
        return 0

    installable_packages = []
    manual_tools = []
    for tool in missing:
        if tool in CORE_TOOLS:
            continue
        package = CONDA_PACKAGES.get(tool)
        if package:
            installable_packages.append(package)
        else:
            manual_tools.append(tool)

    installable_packages = sorted(set(installable_packages))

    if installable_packages:
        manager_path = choose_package_manager(args.manager)
        manager = package_manager_name(manager_path) or ("mamba" if args.manager == "mamba" else "conda")
        suggestion = [manager, "install", "-c", "conda-forge", "-c", "bioconda"] + installable_packages
        print("Conda-installable missing packages:")
        print("  {}".format(" ".join(installable_packages)))
        print()
        print("Suggested command for the current environment:")
        print("  {}".format(" ".join(suggestion)))
        print()
    else:
        manager_path = ""

    if manual_tools:
        print("Manual installation needed:")
        for tool in manual_tools:
            print("  {}: {}".format(tool, MANUAL_TOOLS.get(tool, "Install manually and add it to PATH.")))
        print()

    if not args.install:
        if args.strict:
            print("Strict check failed because one or more tools were not resolved from the current shell.")
        else:
            print("This is an advisory check. Missing entries may still work if your scheduler or shell initializes PATH differently.")
            print("Run `skmer-smk2 doctor --strict` when you want missing entries to return a non-zero exit code.")
        print("Run `skmer-smk2 doctor --install` to install the conda-available missing packages into the active environment.")
        return 1 if args.strict else 0

    if not installable_packages:
        print("No missing tools can be installed automatically by conda/mamba.")
        return 1 if args.strict and manual_tools else 0

    manager_path = choose_package_manager(args.manager)
    if not manager_path:
        print("ERROR: neither mamba nor conda was found in PATH.", file=sys.stderr)
        return 2

    install_cmd = [
        manager_path,
        "install",
        "-y",
        "-c",
        "conda-forge",
        "-c",
        "bioconda",
    ] + installable_packages
    print("Installing into the currently active environment:")
    print("  {}".format(" ".join(install_cmd)))
    print()
    rc = subprocess.call(install_cmd)
    if rc != 0:
        return rc
    if manual_tools:
        print()
        print("Automatic install finished, but these tools still need manual installation:")
        for tool in manual_tools:
            print("  {}".format(tool))
        return 1 if args.strict else 0
    return 0


def init(args):
    outdir = Path(args.output).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if args.with_workflow:
        copy_tree(resource_path("workflow"), outdir / "workflow")
    shutil.copyfile(resource_path("templates", "01_input_check.sh"), outdir / "01_input_check.sh")
    shutil.copyfile(resource_path("templates", "02_input_repair.sh"), outdir / "02_input_repair.sh")
    shutil.copyfile(resource_path("templates", "scan_repair_fastq.sh"), outdir / "scan_repair_fastq.sh")
    shutil.copyfile(resource_path("templates", "scan_repair_fastq_hpc.sh"), outdir / "scan_repair_fastq_hpc.sh")
    shutil.copyfile(resource_path("templates", "submit_example.sh"), outdir / "submit_example.sh")
    print("Wrote templates to {}".format(outdir))
    return 0


def repair_fastq(args):
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    script = workdir / "scan_repair_fastq.sh"
    hpc_script = workdir / "scan_repair_fastq_hpc.sh"
    input_check_script = workdir / "01_input_check.sh"
    input_repair_script = workdir / "02_input_repair.sh"
    shutil.copyfile(resource_path("templates", "01_input_check.sh"), input_check_script)
    shutil.copyfile(resource_path("templates", "02_input_repair.sh"), input_repair_script)
    shutil.copyfile(resource_path("templates", "scan_repair_fastq.sh"), script)
    shutil.copyfile(resource_path("templates", "scan_repair_fastq_hpc.sh"), hpc_script)
    if args.copy_only:
        print("Wrote {}".format(input_check_script))
        print("Wrote {}".format(input_repair_script))
        print("Wrote {}".format(script))
        print("Wrote {}".format(hpc_script))
        return 0
    cmd = ["bash", str(script), str(Path(args.input).resolve()), str(workdir)]
    print("Running:", " ".join(cmd), flush=True)
    return subprocess.call(cmd)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="skmer-smk2",
        description="Packaged Skmer/WASTER/Mash Snakemake workflow.",
    )
    parser.add_argument("--version", action="version", version="skmer-smk2 {}".format(__version__))
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser(
        "run",
        help="Run the workflow.",
        description="Run the packaged workflow. The same command works locally and inside an HPC scheduler script.",
    )
    p_run.add_argument("-i", "--input", required=True, help="Directory containing paired FASTQ files.")
    p_run.add_argument("-ref", "--ref", default="", help="Optional plastid/reference genome FASTA for bowtie2 filtering.")
    p_run.add_argument("-s", "--sample-percentile", type=float, default=75.0, help="Sorted sample base-count percentile used as the head cutoff.")
    p_run.add_argument("--candidate-percentiles", default="50,60,70,75,80,90,95", help="Comma- or whitespace-separated percentiles to report for choosing -s.")
    p_run.add_argument("-skmer", action="store_true", help="Run the Skmer branch. If no branch is selected, all branches run.")
    p_run.add_argument("-waster", action="store_true", help="Run the WASTER branch. If no branch is selected, all branches run.")
    p_run.add_argument("-mash", action="store_true", help="Run the Mash branch. If no branch is selected, all branches run.")
    p_run.add_argument("-j", "--jobs", default="1", help="Snakemake job count/cores.")
    p_run.add_argument("-b", "--bootstraps", type=int, default=100, help="Bootstrap replicate count.")
    p_run.add_argument("--exclude-samples", default="", help="Comma- or whitespace-separated sample names to skip.")
    p_run.add_argument("--bbmerge-timeout", type=int, default=14400, help="Seconds to wait for each BBMerge job before using unmerged-read fallback. Use 0 to disable.")
    p_run.add_argument("--workdir", default=".", help="Run directory for results and workflow cache.")
    p_run.add_argument("--latency-wait", default="120", help="Snakemake latency wait seconds.")
    p_run.add_argument("--dry-run", action="store_true", help="Run Snakemake in dry-run mode.")
    p_run.add_argument("--printshellcmds", action="store_true", help="Print shell commands from Snakemake.")
    p_run.add_argument("snakemake_args", nargs=argparse.REMAINDER, help="Extra arguments passed to snakemake after --.")
    p_run.set_defaults(func=run)

    p_doctor = sub.add_parser(
        "doctor",
        help="Check required external tools and optionally install conda packages.",
        description="Check external tools used by the workflow and print install advice.",
    )
    p_doctor.add_argument("--install", action="store_true", help="Install missing conda-available packages into the active environment.")
    p_doctor.add_argument("--manager", choices=("auto", "mamba", "conda"), default="auto", help="Package manager for --install. Default: auto.")
    p_doctor.add_argument("--strict", action="store_true", help="Return a non-zero exit code when tools are missing.")
    p_doctor.set_defaults(func=doctor)

    p_init = sub.add_parser("init", help="Write optional helper templates into a directory.")
    p_init.add_argument("-o", "--output", default="skmer_smk2_templates")
    p_init.add_argument("--with-workflow", action="store_true", help="Also copy the packaged Snakefile and helper scripts.")
    p_init.set_defaults(func=init)

    p_repair = sub.add_parser("repair-fastq", help="Run or copy the FASTQ scan/repair helper.")
    p_repair.add_argument("-i", "--input", default=".", help="FASTQ directory to scan.")
    p_repair.add_argument("--workdir", default=".", help="Directory where the helper shell script is written.")
    p_repair.add_argument("--copy-only", action="store_true", help="Only copy the helper script.")
    p_repair.set_defaults(func=repair_fastq)

    p_input_check = sub.add_parser("input-check", help="Check paired FASTQ inputs and write detailed reports.")
    p_input_check.add_argument("-i", "--input", required=True, help="Directory containing raw FASTQ files.")
    p_input_check.add_argument("-o", "--output", default="input_qc", help="Output directory for input QC reports.")
    p_input_check.set_defaults(func=input_check_main)

    p_input_repair = sub.add_parser("input-repair", help="Repair selected samples by sample name.")
    p_input_repair.add_argument("-i", "--input", required=True, help="Directory containing raw FASTQ files.")
    p_input_repair.add_argument("-o", "--output", default="input_qc", help="Output directory for repair reports and input_for_skmer.")
    p_input_repair.add_argument("--samples", required=True, help="Comma- or whitespace-separated sample names to repair.")
    p_input_repair.set_defaults(func=input_repair_main)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

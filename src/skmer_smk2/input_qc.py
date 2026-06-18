#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


FASTQ_EXTS = (".fq.gz", ".fastq.gz", ".fq", ".fastq")
PAIR_RE = re.compile(r"^(?P<sample>.+?)(?P<sep>[_.-])(?P<mate>R?[12]|r?[12])$")


FILE_FIELDS = [
    "sample",
    "mate",
    "file",
    "path",
    "size_bytes",
    "mtime",
    "gzip_status",
    "seqkit_status",
    "reads",
    "bases",
    "average_length",
    "note",
]

SAMPLE_FIELDS = [
    "sample",
    "r1_file",
    "r2_file",
    "r1_size_bytes",
    "r2_size_bytes",
    "r1_gzip_status",
    "r2_gzip_status",
    "r1_reads",
    "r2_reads",
    "r1_bases",
    "r2_bases",
    "r1_avg_len",
    "r2_avg_len",
    "pair_status",
    "suggested_action",
    "note",
]

POST_FIELDS = SAMPLE_FIELDS + ["repair_action", "output_r1", "output_r2"]


def split_fastq_ext(path):
    name = Path(path).name
    for ext in FASTQ_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)], ext
    return Path(path).stem, Path(path).suffix


def discover_fastqs(input_dir):
    input_dir = Path(input_dir).resolve()
    rows = []
    samples = {}
    for path in sorted(input_dir.iterdir()):
        if not path.is_file() or not path.name.endswith(FASTQ_EXTS):
            continue
        stem, _ = split_fastq_ext(path)
        match = PAIR_RE.match(stem)
        sample = ""
        mate = ""
        if match:
            sample = match.group("sample")
            mate_raw = match.group("mate").upper()
            mate = "R1" if mate_raw in ("1", "R1") else "R2"
            key = mate.lower()
            pairs = samples.setdefault(sample, {})
            if key in pairs:
                pairs.setdefault("duplicates", []).append(path.resolve())
            else:
                pairs[key] = path.resolve()
        rows.append({"path": path.resolve(), "sample": sample, "mate": mate})
    return rows, dict(sorted(samples.items()))


def run_command(cmd, log_path=None):
    if log_path:
        with open(log_path, "w", encoding="utf-8") as log:
            log.write("Running: {}\n\n".format(" ".join(str(part) for part in cmd)))
            try:
                return subprocess.call(cmd, stdout=log, stderr=log)
            except OSError as exc:
                log.write("Failed to start command: {}\n".format(exc))
                return 127
    try:
        return subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError:
        return 127


def gzip_check(path, log_dir):
    if not str(path).endswith(".gz"):
        return "NA"
    log_path = log_dir / (path.name + ".gzip.log")
    rc = run_command(["gzip", "-t", str(path)], log_path)
    return "OK" if rc == 0 else "BAD"


def seqkit_stats(path):
    try:
        proc = subprocess.run(
            ["seqkit", "stats", "-T", str(path)],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError:
        return {"seqkit_status": "ERROR", "reads": "", "bases": "", "average_length": ""}
    if proc.returncode != 0:
        return {"seqkit_status": "BAD", "reads": "", "bases": "", "average_length": ""}
    lines = [line for line in proc.stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        return {"seqkit_status": "BAD", "reads": "", "bases": "", "average_length": ""}
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    table = dict(zip(header, values))
    reads = clean_number(table.get("num_seqs", ""))
    bases = clean_number(table.get("sum_len", ""))
    avg = table.get("avg_len", "")
    return {
        "seqkit_status": "OK",
        "reads": reads,
        "bases": bases,
        "average_length": clean_number(avg),
    }


def clean_number(value):
    return str(value).replace(",", "").strip()


def stat_file(path):
    try:
        stat = path.stat()
        return str(stat.st_size), str(int(stat.st_mtime))
    except OSError:
        return "", ""


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def evaluate_sample(sample, pair, file_by_path):
    r1 = pair.get("r1")
    r2 = pair.get("r2")
    if not r1 or not r2:
        return "MISSING_PAIR", "MISSING_PAIR", "missing R1 or R2"
    if pair.get("duplicates"):
        return "PAIRED", "CHECK_MANUALLY", "multiple files were detected for the same sample mate"

    r1_row = file_by_path[str(r1)]
    r2_row = file_by_path[str(r2)]
    if r1_row["gzip_status"] == "BAD" or r2_row["gzip_status"] == "BAD":
        return "PAIRED", "REPAIR_CANDIDATE", "gzip check failed; can try repair if selected"
    if r1_row["seqkit_status"] != "OK" or r2_row["seqkit_status"] != "OK":
        return "PAIRED", "REPAIR_CANDIDATE", "seqkit stats failed; try repair if selected"
    if not r1_row["reads"] or r1_row["reads"] != r2_row["reads"]:
        return "PAIRED", "REPAIR_CANDIDATE", "R1/R2 read counts differ"
    return "PAIRED", "USE_AS_IS", "passed input checks"


def process_fastq_item(item, log_dir):
    path = item["path"]
    size, mtime = stat_file(path)
    gzip_status = gzip_check(path, log_dir)
    stats = seqkit_stats(path) if gzip_status != "BAD" else {
        "seqkit_status": "SKIPPED",
        "reads": "",
        "bases": "",
        "average_length": "",
    }
    note = ""
    if not item["sample"] or not item["mate"]:
        note = "unsupported FASTQ mate naming style"
    elif gzip_status == "BAD":
        note = "gzip -t failed"
    elif stats["seqkit_status"] != "OK":
        note = "seqkit stats failed"
    return {
        "sample": item["sample"],
        "mate": item["mate"],
        "file": path.name,
        "path": str(path),
        "size_bytes": size,
        "mtime": mtime,
        "gzip_status": gzip_status,
        "note": note,
        **stats,
    }


def build_reports(input_dir, output_dir, jobs=1):
    if not shutil.which("gzip"):
        raise RuntimeError("gzip was not found in PATH")
    if not shutil.which("seqkit"):
        raise RuntimeError("seqkit was not found in PATH")

    output_dir = Path(output_dir).resolve()
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    fastqs, samples = discover_fastqs(input_dir)

    jobs = max(1, int(jobs))
    if jobs == 1:
        file_rows = [process_fastq_item(item, log_dir) for item in fastqs]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            file_rows = list(pool.map(lambda item: process_fastq_item(item, log_dir), fastqs))

    file_by_path = {row["path"]: row for row in file_rows}
    sample_rows = []
    for sample, pair in samples.items():
        r1 = pair.get("r1")
        r2 = pair.get("r2")
        pair_status, suggested_action, note = evaluate_sample(sample, pair, file_by_path)
        r1_row = file_by_path.get(str(r1), {})
        r2_row = file_by_path.get(str(r2), {})
        sample_rows.append({
            "sample": sample,
            "r1_file": Path(r1).name if r1 else "",
            "r2_file": Path(r2).name if r2 else "",
            "r1_size_bytes": r1_row.get("size_bytes", ""),
            "r2_size_bytes": r2_row.get("size_bytes", ""),
            "r1_gzip_status": r1_row.get("gzip_status", ""),
            "r2_gzip_status": r2_row.get("gzip_status", ""),
            "r1_reads": r1_row.get("reads", ""),
            "r2_reads": r2_row.get("reads", ""),
            "r1_bases": r1_row.get("bases", ""),
            "r2_bases": r2_row.get("bases", ""),
            "r1_avg_len": r1_row.get("average_length", ""),
            "r2_avg_len": r2_row.get("average_length", ""),
            "pair_status": pair_status,
            "suggested_action": suggested_action,
            "note": note,
        })

    paired = set(samples)
    for row in file_rows:
        if row["sample"] and row["sample"] not in paired:
            sample_rows.append({
                "sample": row["sample"],
                "r1_file": row["file"] if row["mate"] == "R1" else "",
                "r2_file": row["file"] if row["mate"] == "R2" else "",
                "pair_status": "MISSING_PAIR",
                "suggested_action": "MISSING_PAIR",
                "note": "only one mate found",
            })

    write_tsv(output_dir / "input_file_report.tsv", FILE_FIELDS, file_rows)
    write_tsv(output_dir / "input_sample_report.tsv", SAMPLE_FIELDS, sorted(sample_rows, key=lambda r: r["sample"]))
    return file_rows, sample_rows


def parse_samples(value):
    return set(token for token in re.split(r"[,\s]+", value.strip()) if token)


def link_or_copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.symlink(str(Path(src).resolve()), str(dst))
    except OSError:
        shutil.copyfile(src, dst)


def run_sana(src, dst, log_dir):
    log = log_dir / (Path(src).name + ".sana.log")
    return run_command(["seqkit", "sana", str(src), "-o", str(dst)], log)


def run_repair(r1, r2, out1, out2, singletons, sample, log_dir, threads=2):
    log = log_dir / (sample + ".repair.log")
    cmd = [
        "repair.sh",
        "in1={}".format(r1),
        "in2={}".format(r2),
        "out1={}".format(out1),
        "out2={}".format(out2),
        "outs={}".format(singletons),
        "repair=t",
        "overwrite=t",
        "threads={}".format(threads),
    ]
    return run_command(cmd, log)


def repair_one_sample(sample, row, pairs, selected, repaired_dir, input_for_skmer, log_dir, repair_threads):
    row = dict(row)
    pair = pairs.get(sample, {})
    r1 = pair.get("r1")
    r2 = pair.get("r2")
    row["repair_action"] = "SKIPPED"
    row["output_r1"] = ""
    row["output_r2"] = ""

    if row["suggested_action"] == "USE_AS_IS" and sample not in selected and r1 and r2:
        out1 = input_for_skmer / Path(r1).name
        out2 = input_for_skmer / Path(r2).name
        link_or_copy(r1, out1)
        link_or_copy(r2, out2)
        row["repair_action"] = "LINKED_ORIGINAL"
        row["output_r1"] = str(out1)
        row["output_r2"] = str(out2)
        return row

    if sample not in selected:
        return row

    if not r1 or not r2:
        row["repair_action"] = "REPAIR_FAILED"
        row["note"] = "selected for repair but missing R1 or R2"
        return row

    sana_r1 = repaired_dir / "{}.R1.sana.fq.gz".format(sample)
    sana_r2 = repaired_dir / "{}.R2.sana.fq.gz".format(sample)
    repaired_r1 = repaired_dir / "{}.R1.repaired.fq.gz".format(sample)
    repaired_r2 = repaired_dir / "{}.R2.repaired.fq.gz".format(sample)
    singletons = repaired_dir / "{}.singletons.fq.gz".format(sample)

    sana_rc = run_sana(r1, sana_r1, log_dir)
    sana_rc2 = run_sana(r2, sana_r2, log_dir)
    if sana_rc != 0 or sana_rc2 != 0:
        row["repair_action"] = "SANA_FAILED"
        row["note"] = "seqkit sana failed"
        return row

    repair_rc = run_repair(sana_r1, sana_r2, repaired_r1, repaired_r2, singletons, sample, log_dir, repair_threads)
    stats1 = seqkit_stats(repaired_r1)
    stats2 = seqkit_stats(repaired_r2)
    if (
        repair_rc == 0
        and stats1["seqkit_status"] == "OK"
        and stats2["seqkit_status"] == "OK"
        and stats1["reads"]
        and stats1["reads"] == stats2["reads"]
    ):
        out1 = input_for_skmer / (Path(r1).name)
        out2 = input_for_skmer / (Path(r2).name)
        link_or_copy(repaired_r1, out1)
        link_or_copy(repaired_r2, out2)
        row["r1_reads"] = stats1["reads"]
        row["r2_reads"] = stats2["reads"]
        row["r1_bases"] = stats1["bases"]
        row["r2_bases"] = stats2["bases"]
        row["r1_avg_len"] = stats1["average_length"]
        row["r2_avg_len"] = stats2["average_length"]
        if row["r1_gzip_status"] == "BAD" or row["r2_gzip_status"] == "BAD":
            row["repair_action"] = "GZIP_BAD_REPAIRED_WARN"
            row["note"] = "gzip was BAD before repair; repaired output was generated"
        else:
            row["repair_action"] = "REPAIRED"
            row["note"] = "repaired output was generated"
        row["output_r1"] = str(out1)
        row["output_r2"] = str(out2)
    else:
        row["repair_action"] = "REPAIR_FAILED"
        row["note"] = "repair.sh failed or repaired R1/R2 counts are unequal"
    return row


def repair_inputs(input_dir, output_dir, samples_text, jobs=1, repair_threads=2):
    output_dir = Path(output_dir).resolve()
    input_dir = Path(input_dir).resolve()
    selected = parse_samples(samples_text)
    if not selected:
        print("ERROR: --samples must include at least one sample name", file=sys.stderr)
        return 2

    _, sample_rows = build_reports(input_dir, output_dir, jobs=jobs)
    sample_map = {row["sample"]: row for row in sample_rows}
    _, pairs = discover_fastqs(input_dir)
    log_dir = output_dir / "logs"
    repaired_dir = output_dir / "repaired"
    input_for_skmer = output_dir / "input_for_skmer"
    repaired_dir.mkdir(parents=True, exist_ok=True)
    input_for_skmer.mkdir(parents=True, exist_ok=True)

    jobs = max(1, int(jobs))
    repair_threads = max(1, int(repair_threads))
    tasks = [
        (sample, sample_map[sample])
        for sample in sorted(sample_map)
    ]
    if jobs == 1:
        post_rows = [
            repair_one_sample(sample, row, pairs, selected, repaired_dir, input_for_skmer, log_dir, repair_threads)
            for sample, row in tasks
        ]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = [
                pool.submit(
                    repair_one_sample,
                    sample,
                    row,
                    pairs,
                    selected,
                    repaired_dir,
                    input_for_skmer,
                    log_dir,
                    repair_threads,
                )
                for sample, row in tasks
            ]
            post_rows = [future.result() for future in futures]

    unknown = sorted(selected - set(sample_map))
    for sample in unknown:
        post_rows.append({
            "sample": sample,
            "pair_status": "MISSING_PAIR",
            "suggested_action": "MISSING_PAIR",
            "repair_action": "REPAIR_FAILED",
            "note": "sample name was selected but not found",
        })

    post_file_rows, _ = build_reports(input_for_skmer, output_dir / "post_repair_check", jobs=jobs)
    write_tsv(output_dir / "post_repair_sample_report.tsv", POST_FIELDS, post_rows)
    write_tsv(output_dir / "post_repair_file_report.tsv", FILE_FIELDS, post_file_rows)
    return 0


def input_check_main(args):
    try:
        build_reports(args.input, args.output, jobs=args.jobs)
    except RuntimeError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 127
    print("Input check finished.")
    print("File report: {}".format(Path(args.output).resolve() / "input_file_report.tsv"))
    print("Sample report: {}".format(Path(args.output).resolve() / "input_sample_report.tsv"))
    return 0


def input_repair_main(args):
    if not shutil.which("repair.sh"):
        print("ERROR: repair.sh was not found in PATH", file=sys.stderr)
        return 127
    try:
        rc = repair_inputs(args.input, args.output, args.samples, jobs=args.jobs, repair_threads=args.repair_threads)
    except RuntimeError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 127
    if rc == 0:
        print("Input repair finished.")
        print("Post-repair report: {}".format(Path(args.output).resolve() / "post_repair_sample_report.tsv"))
        print("Skmer-ready input: {}".format(Path(args.output).resolve() / "input_for_skmer"))
    return rc

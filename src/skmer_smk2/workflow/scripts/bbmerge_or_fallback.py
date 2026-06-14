#!/usr/bin/env python3
import argparse
import gzip
import os
import signal
import shutil
import subprocess
import sys


def open_input(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def copy_fastq(src, dst):
    with open_input(src) as inp, open(dst, "wb") as out:
        shutil.copyfileobj(inp, out, length=1024 * 1024)


def remove_outputs(paths):
    for path in paths:
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def outputs_exist(paths):
    return all(os.path.exists(path) for path in paths)


def log_message(log, message):
    with open(log, "a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")


def stop_process(proc, log):
    if os.name == "nt":
        try:
            subprocess.call(
                ["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                stdout=log,
                stderr=log,
            )
        except OSError:
            proc.kill()
    else:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass

    try:
        proc.wait(timeout=30)
        return
    except subprocess.TimeoutExpired:
        log.write("BBMerge did not stop after SIGTERM; killing it.\n")

    if os.name == "nt":
        try:
            subprocess.call(
                ["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                stdout=log,
                stderr=log,
            )
        except OSError:
            proc.kill()
    else:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    proc.wait()


def run_bbmerge(args):
    cmd = [
        args.bbmerge,
        "-da",
        f"t={args.threads}",
        "overwrite=t",
        f"in1={args.r1}",
        f"in2={args.r2}",
        f"out={args.merged}",
        f"outu1={args.unmerged1}",
        f"outu2={args.unmerged2}",
    ]
    with open(args.log, "w", encoding="utf-8") as log:
        log.write("Running BBMerge:\n")
        log.write(" ".join(cmd) + "\n\n")
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=log,
                stderr=log,
                start_new_session=(os.name != "nt"),
            )
        except OSError as exc:
            log.write(f"Failed to start BBMerge: {exc}\n")
            return 127
        timeout = None if args.timeout <= 0 else args.timeout
        try:
            proc.communicate(timeout=timeout)
            return proc.returncode
        except subprocess.TimeoutExpired:
            log.write(f"\nBBMerge timed out after {args.timeout} seconds.\n")
            stop_process(proc, log)
            return 124


def bbmerge_log_has_exception(log):
    if not os.path.exists(log):
        return False
    with open(log, "r", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    return "Exception in thread" in text or "java.lang." in text


def write_fallback(args, reason):
    remove_outputs([args.merged, args.unmerged1, args.unmerged2])
    open(args.merged, "wb").close()
    copy_fastq(args.r1, args.unmerged1)
    copy_fastq(args.r2, args.unmerged2)
    log_message(args.log, "")
    log_message(args.log, "BBMERGE_FALLBACK_USED")
    log_message(args.log, f"reason: {reason}")
    log_message(args.log, "sample_merged.fq was written as an empty file.")
    log_message(args.log, "R1/R2 were copied to sample_unmerged1.fq and sample_unmerged2.fq.")


def main():
    parser = argparse.ArgumentParser(description="Run bbmerge.sh, falling back to unmerged reads if it fails or times out.")
    parser.add_argument("--bbmerge", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--merged", required=True)
    parser.add_argument("--unmerged1", required=True)
    parser.add_argument("--unmerged2", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=14400, help="BBMerge timeout in seconds. Use 0 to disable.")
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.merged) or ".", exist_ok=True)
    remove_outputs([args.merged, args.unmerged1, args.unmerged2])

    rc = run_bbmerge(args)
    outputs = [args.merged, args.unmerged1, args.unmerged2]
    if rc == 124:
        write_fallback(args, f"bbmerge timed out after {args.timeout} seconds")
    elif rc != 0:
        write_fallback(args, f"bbmerge exited with status {rc}")
    elif bbmerge_log_has_exception(args.log):
        write_fallback(args, "bbmerge reported a Java exception")
    elif not outputs_exist(outputs):
        missing = ", ".join(path for path in outputs if not os.path.exists(path))
        write_fallback(args, f"bbmerge did not create all expected outputs: {missing}")
    else:
        log_message(args.log, "")
        log_message(args.log, "BBMERGE_SUCCESS")

    if not outputs_exist(outputs):
        print("ERROR: fallback did not create all expected bbmerge outputs", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

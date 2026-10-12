#!/usr/bin/env python3
"""Validate raw distances and archive invalid downstream Skmer outputs."""
import argparse
from datetime import datetime
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src/skmer_smk2/workflow/scripts"))
from distance_to_phylip import read_table, parse_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--replicates", type=int, default=100)
    args = parser.parse_args()
    root = Path(args.workdir).resolve()
    skmer = root / "results/skmer"
    raw_main = skmer / "dimtrx_main.txt"
    rows = read_table(raw_main)
    samples = rows[0][1:]
    parse_matrix(rows, samples)
    for i in range(args.replicates):
        rep = skmer / "subsample" / ("rep%d" % i)
        parse_matrix(read_table(rep / "dimtrx_rep.txt"), samples)
        for name in ("block_size.npy", "sample_size.npy"):
            path = rep / name
            if not path.is_file() or not path.stat().st_size:
                raise ValueError("Missing metadata: " + str(path))
    print("Validated raw main matrix and {} replicate matrices; {} samples".format(args.replicates, len(samples)))
    backup = root / ("skmer_invalid_backup_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
    paths = [p for p in skmer.iterdir() if p.is_file() and p.name not in ("dimtrx_main.txt", "dimtrx_main.txt.txt")]
    paths += [skmer / "subsample/_summary.csv"]
    for i in range(args.replicates):
        rep = skmer / "subsample" / ("rep%d" % i)
        paths.extend(rep / name for name in ("dimtrx_rep_cor.txt", "dimtrx_rep_cor_cons.txt", "dimtrx_rep.phy", "tree.tre"))
    count = 0
    for path in paths:
        if not path.exists():
            continue
        source = path.resolve()
        if skmer.resolve() not in source.parents:
            raise ValueError("Output path escapes results/skmer: " + str(source))
        dest = backup / path.relative_to(root)
        if root not in dest.resolve().parents:
            raise ValueError("Backup path escapes workdir")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(dest))
        count += 1
    print("Archived {} downstream files to {}".format(count, backup))
    print("Raw matrices, NumPy metadata, logs and preprocessing files were retained.")
    print("Now resume with version 0.1.33 and --rerun-incomplete --rerun-triggers mtime.")


if __name__ == "__main__":
    main()

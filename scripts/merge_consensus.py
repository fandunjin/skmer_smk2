#!/usr/bin/env python3
import argparse
import re


def read_tree(path):
    with open(path, "r", encoding="utf-8") as handle:
        return "".join(line.strip() for line in handle if line.strip())


def normalize_raxml_support(tree):
    return re.sub(r"(:[0-9.eE+-]+)\[([0-9.]+)\]", r"\2\1", tree)


def main():
    parser = argparse.ArgumentParser(description="Merge direct and bootstrap Newick trees.")
    parser.add_argument("--bootstrap", required=True)
    parser.add_argument("--direct", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    bootstrap = normalize_raxml_support(read_tree(args.bootstrap))
    direct = read_tree(args.direct)

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("[direct]\n")
        out.write(direct.rstrip(";") + ";\n")
        out.write("[bootstrap]\n")
        out.write(bootstrap.rstrip(";") + ";\n")


if __name__ == "__main__":
    main()

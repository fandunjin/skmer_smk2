#!/usr/bin/env python3
import argparse
import re
from dataclasses import dataclass, field


@dataclass
class Node:
    children: list = field(default_factory=list)
    label: str = ""
    length: str = ""


def read_tree(path):
    with open(path, "r", encoding="utf-8") as handle:
        return "".join(line.strip() for line in handle if line.strip())


def write_tree(path, tree):
    with open(path, "w", encoding="utf-8") as out:
        out.write(tree.rstrip(";") + ";\n")


def normalize_raxml_support(tree):
    # RAxML MRE can write support after branch length, e.g. :1.0[59].
    # Many viewers parse that as an invalid branch length, so use the standard
    # internal-node-label form: )59:1.0.
    return re.sub(
        r"(\))([^,():;\[\]]*)(:[0-9.eE+-]+)\[([0-9.]+)\]",
        r"\1\4\3",
        tree,
    )


def parse_newick(tree):
    text = normalize_raxml_support(tree).strip().rstrip(";")
    idx = 0

    def parse_label_length(pos):
        label_start = pos
        while pos < len(text) and text[pos] not in ",():;":
            pos += 1
        label = text[label_start:pos].strip()
        length = ""
        if pos < len(text) and text[pos] == ":":
            pos += 1
            length_start = pos
            while pos < len(text) and text[pos] not in ",();":
                pos += 1
            length = text[length_start:pos].strip()
        return label, length, pos

    def parse_node(pos):
        if pos < len(text) and text[pos] == "(":
            pos += 1
            children = []
            while True:
                child, pos = parse_node(pos)
                children.append(child)
                if pos >= len(text):
                    raise ValueError("Unexpected end of Newick tree inside a clade")
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Unexpected character in Newick tree: {text[pos]!r}")
            label, length, pos = parse_label_length(pos)
            return Node(children=children, label=label, length=length), pos

        label, length, pos = parse_label_length(pos)
        if not label:
            raise ValueError("Leaf node without a label in Newick tree")
        return Node(label=label, length=length), pos

    root, idx = parse_node(idx)
    if idx != len(text):
        raise ValueError(f"Unexpected trailing Newick text: {text[idx:]!r}")
    return root


def serialize(node):
    if node.children:
        text = "(" + ",".join(serialize(child) for child in node.children) + ")"
        if node.label:
            text += node.label
    else:
        text = node.label
    if node.length:
        text += ":" + node.length
    return text


def clade_supports(root):
    supports = {}

    def visit(node):
        if not node.children:
            return frozenset([node.label])
        tips = frozenset().union(*(visit(child) for child in node.children))
        if node.label:
            supports[tips] = node.label
        return tips

    visit(root)
    return supports


def apply_supports(root, supports):
    def visit(node):
        if not node.children:
            return frozenset([node.label])
        tips = frozenset().union(*(visit(child) for child in node.children))
        if tips in supports:
            node.label = supports[tips]
        return tips

    visit(root)
    return root


def normalize_command(args):
    write_tree(args.output, normalize_raxml_support(read_tree(args.input)))


def merge_command(args):
    bootstrap = parse_newick(read_tree(args.bootstrap))
    direct = parse_newick(read_tree(args.direct))
    supports = clade_supports(bootstrap)
    merged = apply_supports(direct, supports)
    write_tree(args.output, serialize(merged))


def main():
    parser = argparse.ArgumentParser(
        description="Normalize RAxML support labels or merge bootstrap support onto a direct Newick tree."
    )
    subparsers = parser.add_subparsers(dest="command")

    normalize_parser = subparsers.add_parser(
        "normalize",
        help="Convert RAxML :branch[support] annotations to internal node labels.",
    )
    normalize_parser.add_argument("--input", required=True)
    normalize_parser.add_argument("--output", required=True)
    normalize_parser.set_defaults(func=normalize_command)

    merge_parser = subparsers.add_parser(
        "merge",
        help="Write one direct-topology tree with matching bootstrap support labels.",
    )
    merge_parser.add_argument("--bootstrap", required=True)
    merge_parser.add_argument("--direct", required=True)
    merge_parser.add_argument("--output", required=True)
    merge_parser.set_defaults(func=merge_command)

    # Backward-compatible interface used by older Snakefiles.
    parser.add_argument("--bootstrap")
    parser.add_argument("--direct")
    parser.add_argument("--output")
    args = parser.parse_args()

    if args.command:
        args.func(args)
    elif args.bootstrap and args.direct and args.output:
        merge_command(args)
    else:
        parser.error("use a subcommand, or provide --bootstrap --direct --output")


if __name__ == "__main__":
    main()

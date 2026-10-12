#!/usr/bin/env python3
"""Run Skmer's correction algorithm with modern pandas and useful tracebacks.

Replace the removed DataFrame.append call and normalize numeric NumPy scalar
strings before correction. The installed Skmer source is never modified;
all statistical calculations remain upstream's implementation.
"""
import argparse
import ast
import inspect
import textwrap
from pathlib import Path
from distance_to_phylip import parse_matrix


class ValidatedPandas:
    """Normalize serialized NumPy scalars before upstream coerces them to NaN."""
    def __init__(self, pandas):
        self.pandas = pandas

    def __getattr__(self, name):
        return getattr(self.pandas, name)

    def read_csv(self, path, *args, **kwargs):
        frame = self.pandas.read_csv(path, *args, **kwargs)
        try:
            samples = list(frame.iloc[:, 0].astype(str))
            rows = [list(frame.columns)] + frame.astype(object).values.tolist()
            matrix = parse_matrix(rows, samples)
            for column in frame.columns[1:]:
                frame[column] = [matrix[sample][column] for sample in samples]
        except Exception as exc:
            raise ValueError("Invalid distance matrix {}: {}".format(path, exc)) from exc
        return frame


class CorrectionCompatibility(ast.NodeTransformer):
    def visit_Call(self, node):
        self.generic_visit(node)
        if (isinstance(node.func, ast.Attribute) and node.func.attr == "append"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "combo_result_df"):
            if len(node.args) != 1 or len(node.keywords) != 1 or node.keywords[0].arg != "ignore_index":
                raise RuntimeError("Unsupported Skmer append signature; inspect the installed correction()")
            return ast.copy_location(ast.Call(
                func=ast.Attribute(value=ast.Name(id="pd", ctx=ast.Load()), attr="concat", ctx=ast.Load()),
                args=[ast.List(elts=[node.func.value, node.args[0]], ctx=ast.Load())],
                keywords=node.keywords,
            ), node)
        return node

    def visit_Raise(self, node):
        # Upstream's broad exception handler otherwise hides missing metadata,
        # pandas failures, and unreadable NumPy files behind a directory error.
        if (isinstance(node.exc, ast.Call) and isinstance(node.exc.func, ast.Name)
                and node.exc.func.id == "ValueError" and node.exc.args
                and isinstance(node.exc.args[0], ast.Constant)
                and node.exc.args[0].value == "Please check subsample directory and try again"):
            return ast.copy_location(ast.Raise(exc=None, cause=None), node)
        return node


def compatible_correction(module):
    tree = ast.parse(textwrap.dedent(inspect.getsource(module.correction)))
    tree = ast.fix_missing_locations(CorrectionCompatibility().visit(tree))
    namespace = dict(vars(module))
    namespace["pd"] = ValidatedPandas(module.pd)
    exec(compile(tree, "<skmer correction compatibility>", "exec"), namespace)
    return namespace["correction"]


def run(main, sub, replicates):
    # Import in the same environment used for the workflow Python tools.
    import skmer.__main__ as skmer

    sub = Path(sub)
    expected = {"rep%d" % i for i in range(replicates)}
    observed = {p.name for p in sub.iterdir() if "rep" in p.name}
    if observed != expected:
        raise ValueError("Replicate directories differ: missing={}, unexpected={}".format(
            sorted(expected - observed), sorted(observed - expected)))
    for name in sorted(expected):
        for filename in ("dimtrx_rep.txt", "block_size.npy", "sample_size.npy"):
            path = sub / name / filename
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError("Missing or empty correction input: {}".format(path))
    compatible_correction(skmer)(argparse.Namespace(main=str(main), sub=str(sub)))
    # The upstream calculation uses errors='coerce'; never let NaN output
    # reach FastME, even when malformed metadata caused the missing values.
    reader = ValidatedPandas(skmer.pd)
    reader.read_csv(Path(main).with_name(Path(main).stem + "_cor_.txt"), sep="\t")
    for name in sorted(expected):
        path = sub / name / "dimtrx_rep_cor.txt"
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError("Skmer did not write corrected replicate: {}".format(path))
        reader.read_csv(path, sep="\t")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main", required=True)
    parser.add_argument("--sub", required=True)
    parser.add_argument("--replicates", type=int, required=True)
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be positive")
    run(args.main, args.sub, args.replicates)


if __name__ == "__main__":
    main()

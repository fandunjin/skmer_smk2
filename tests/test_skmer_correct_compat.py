import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "src/skmer_smk2/workflow/scripts/skmer_correct_compat.py"
spec = importlib.util.spec_from_file_location("skmer_correct_compat", SCRIPT)
compat = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compat)


class CorrectionTests(unittest.TestCase):
    def setUp(self):
        try:
            import skmer.__main__
        except ImportError:
            self.skipTest("Install Skmer and its Python dependencies for integration tests")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.main = self.root / "main.txt"
        self.sub = self.root / "subsample"
        self.write_matrix(self.main, 0.2)
        for i, distance in enumerate((0.1, 0.3)):
            rep = self.sub / ("rep%d" % i)
            rep.mkdir(parents=True)
            self.write_matrix(rep / "dimtrx_rep.txt", distance)
            np.save(rep / "block_size.npy", {"a": 25, "b": 25})
            np.save(rep / "sample_size.npy", {"a": 100, "b": 100})

    def write_matrix(self, path, distance):
        frame = pd.DataFrame([[0, distance], [distance, 0]], index=["a", "b"], columns=["a", "b"])
        frame.index.name = "sample"
        frame.to_csv(path, sep="\t")

    def test_modern_pandas_and_correction_formula(self):
        compat.run(self.main, self.sub, 2)
        for i, expected in enumerate((0.15, 0.25)):
            # sqrt(25/100) * (replicate - replicate_mean) + main
            frame = pd.read_csv(self.sub / ("rep%d" % i) / "dimtrx_rep_cor.txt", sep="\t", index_col=0)
            np.testing.assert_allclose(frame.to_numpy(), [[0, expected], [expected, 0]])
        # Repeating correction must not accumulate or change distances.
        before = (self.sub / "rep0/dimtrx_rep_cor.txt").read_bytes()
        compat.run(self.main, self.sub, 2)
        self.assertEqual(before, (self.sub / "rep0/dimtrx_rep_cor.txt").read_bytes())

    def test_missing_metadata_names_exact_file(self):
        (self.sub / "rep1/block_size.npy").unlink()
        with self.assertRaisesRegex(ValueError, "block_size.npy"):
            compat.run(self.main, self.sub, 2)

    def test_corrupt_metadata_preserves_real_exception(self):
        (self.sub / "rep0/block_size.npy").write_bytes(b"invalid numpy data")
        with self.assertRaises(Exception) as caught:
            compat.run(self.main, self.sub, 2)
        self.assertNotIn("Please check subsample directory", str(caught.exception))

    def test_unexpected_replicate_not_silently_combined(self):
        (self.sub / "rep2").mkdir()
        with self.assertRaisesRegex(ValueError, "unexpected=.*rep2"):
            compat.run(self.main, self.sub, 2)


if __name__ == "__main__":
    unittest.main()

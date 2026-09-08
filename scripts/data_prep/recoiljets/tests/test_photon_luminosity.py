"""Synthetic checks of luminosity arithmetic and input conventions."""
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "photon_luminosity", Path(__file__).resolve().parents[1] / "photon_luminosity.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class LuminosityTests(unittest.TestCase):
    def test_ppg12_observed_count_identity(self):
        # Independent expression in observed counts and prescales.
        expected = 950 * 100 * 1.1 / 4 / 0.0252
        result = module.photon_luminosity(100000, 8000, 2000, 0.0252,
                                          pileup=1.1, mb_coverage=950 / 1000)
        self.assertAlmostEqual(result, expected)

    def test_zero_exposure(self):
        self.assertEqual(module.photon_luminosity(100, 0, 0, 1), 0)

    def test_invalid_inputs(self):
        for values in ((1, 2, 3, 1), (1, 0, 1, 1), (-1, 1, 1, 1),
                       (1, 1, 1, 0), (math.nan, 1, 1, 1)):
            with self.subTest(values=values), self.assertRaises(ValueError):
                module.photon_luminosity(*values)

    def test_csv_units_and_duplicate_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.csv"
            header = "run,mb_live,photon_live,photon_scaled,pileup,mb_coverage\n"
            row = "1,1000000000,100,50,1,1\n"
            path.write_text(header + row)
            self.assertEqual(module.calculate_csv(path, 1, "nb^-1")["luminosity"], 0.5)
            self.assertEqual(module.calculate_csv(path, 1, "pb^-1")["luminosity"], 0.0005)
            path.write_text(header + row + row)
            with self.assertRaises(ValueError):
                module.calculate_csv(path, 1, "nb^-1")


if __name__ == "__main__":
    unittest.main()

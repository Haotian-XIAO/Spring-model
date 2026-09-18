from __future__ import annotations
from pathlib import Path
import tempfile
import unittest

from analysis.load_results import discover_seeds
from analysis.plotting import timecourses
from analysis.collaborator_package import classify


class AnalysisRobustnessTests(unittest.TestCase):
    def test_seed_discovery_follows_symlinked_batch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            batch = root / "real_batch"
            seed = batch / "seed_1001"
            (seed / "control").mkdir(parents=True)
            (seed / "softening").mkdir()
            (seed / "pair_metadata.json").write_text("{}")
            for path in (seed / "pair_timeseries.csv", seed / "control" / "timeseries.csv", seed / "softening" / "timeseries.csv"):
                path.write_text("phase\n0\n")
            linked = root / "linked_root"
            linked.mkdir()
            (linked / "batch_01").symlink_to(batch, target_is_directory=True)
            runs = discover_seeds(linked)
            self.assertEqual([run.seed for run in runs], [1001])
            self.assertEqual(runs[0].status, "success")

    def test_empty_timecourse_has_actionable_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "no successful paired time-course rows"):
                timecourses([], Path(temporary))

    def test_operational_region_boundaries_are_inclusive(self) -> None:
        import numpy as np
        np.testing.assert_array_equal(classify(np.array([0.5, 0.75, 1.0, 1.25, 1.5])), [0, 1, 1, 1, 2])


if __name__ == "__main__":
    unittest.main()

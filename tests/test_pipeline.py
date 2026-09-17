"""Fast invariants for the reproducible paired FF experiment pipeline."""

from __future__ import annotations

import csv
from dataclasses import fields
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from pulmonary_fibrosis_model import generate_initial_state, run_pair, smoke_config
from pulmonary_fibrosis_model.core import apply_softening, polygon_areas, state_fingerprint


ROOT = Path(__file__).resolve().parents[1]


class PairedPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = smoke_config()

    def test_same_seed_reproduces_identical_initial_state(self) -> None:
        first, _, _, _ = generate_initial_state(11, self.config)
        second, _, _, _ = generate_initial_state(11, self.config)
        self.assertEqual(state_fingerprint(first), state_fingerprint(second))

    def test_deep_copy_has_no_shared_mutable_state(self) -> None:
        initial, _, _, _ = generate_initial_state(12, self.config)
        control = initial.copy()
        softening = initial.copy()
        self.assertEqual(state_fingerprint(control), state_fingerprint(softening))
        for field in fields(initial):
            self.assertFalse(
                np.shares_memory(getattr(control, field.name), getattr(softening, field.name)),
                field.name,
            )

    def test_on_off_start_identically_and_only_on_softens(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = run_pair(13, self.config, root)
            self.assertTrue(result["same_initial_state"])
            self.assertEqual(
                result["common_initial_state_fingerprint"],
                self._read_json_value(root / "seed_13" / "pair_metadata.json", "initial_state_fingerprint"),
            )
            control0 = np.load(root / "seed_13" / "control" / "spatial_phase_0000.npz")
            softening0 = np.load(root / "seed_13" / "softening" / "spatial_phase_0000.npz")
            for key in ("points", "spring_stiffness", "E", "A", "D", "activation", "initial_fibrotic"):
                np.testing.assert_array_equal(control0[key], softening0[key], err_msg=key)

            control1 = np.load(root / "seed_13" / "control" / "spatial_phase_0001.npz")
            softening1 = np.load(root / "seed_13" / "softening" / "spatial_phase_0001.npz")
            self.assertEqual(int(np.count_nonzero(control1["ever_softened"])), 0)
            triggered = softening1["ever_softened"].astype(bool)
            self.assertGreater(int(np.count_nonzero(triggered)), 0)
            self.assertTrue(np.all(softening1["E"][triggered] < control1["E"][triggered]))
            self.assertTrue(
                np.all(softening1["spring_stiffness"][triggered] < control1["spring_stiffness"][triggered])
            )
            # The intervention itself never writes A; phase-one A is therefore exact-match.
            np.testing.assert_array_equal(softening1["A"], control1["A"])

    def test_a_is_not_modified_and_springs_trigger_at_most_once(self) -> None:
        state, _, _, _ = generate_initial_state(14, self.config)
        area_before = state.A.copy()
        apply_softening(state, self.config, True)
        triggered = state.ever_softened.copy()
        pre_event = state.E_pre_softening.copy()
        target = state.E_target.copy()
        self.assertGreater(int(np.count_nonzero(triggered)), 0)
        np.testing.assert_array_equal(state.A, area_before)
        # Continue through hold/recovery; a spring may never be reset and retriggered.
        for _ in range(8):
            apply_softening(state, self.config, True)
        np.testing.assert_array_equal(state.ever_softened, triggered)
        np.testing.assert_array_equal(state.E_pre_softening, pre_event)
        np.testing.assert_array_equal(state.E_target, target)

    def test_cell_area_calculation(self) -> None:
        points = np.asarray(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)))
        cells = np.asarray(((0, 1, 2, 3),))
        np.testing.assert_allclose(polygon_areas(points, cells), (1.0,))

    def test_cli_completes_tiny_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "cli_output"
            command = [
                sys.executable,
                str(ROOT / "scripts" / "run_pairs.py"),
                "15",
                "--config",
                "smoke",
                "--softening-magnitude",
                "0.5",
                "--output-root",
                str(output),
            ]
            completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertTrue((output / "seed_15" / "pair_summary.csv").exists())
            with (output / "batch_summary.csv").open(newline="", encoding="utf-8") as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row["status"], "success")

    @staticmethod
    def _read_json_value(path: Path, key: str) -> str:
        import json

        return json.loads(path.read_text(encoding="utf-8"))[key]


if __name__ == "__main__":
    unittest.main(verbosity=2)

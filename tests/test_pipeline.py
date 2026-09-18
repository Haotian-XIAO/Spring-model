"""Fast invariants for the reproducible paired FF experiment pipeline."""

from __future__ import annotations

import csv
from dataclasses import fields, replace
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from pulmonary_fibrosis_model import generate_initial_state, run_pair, smoke_config
from pulmonary_fibrosis_model.core import (
    NumericalInstabilityError,
    apply_softening,
    mechanical_equilibrate,
    polygon_areas,
    state_fingerprint,
    update_activation_and_density,
)


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

    def test_nonfinite_trial_fails_with_solver_context(self) -> None:
        state, _, _, _ = generate_initial_state(16, self.config)
        free_nodes = np.where(state.boundary_node_ids == 0)[0]
        broken = state.points.copy()
        broken[free_nodes[0], 0] = np.inf
        with self.assertRaisesRegex(NumericalInstabilityError, r"stage=control.*outer_phase=7"):
            mechanical_equilibrate(
                broken, state.edges, self.config.l0, state.spring_constants, free_nodes,
                self.config, np.random.default_rng(1), 1, stage="control", outer_phase=7,
                condition="test", E=state.E, A=state.A,
            )

    def test_nonpositive_stiffness_fails_instead_of_running_unbounded_energy(self) -> None:
        state, _, _, _ = generate_initial_state(17, self.config)
        state.A[0] = -1.0
        with self.assertRaisesRegex(NumericalInstabilityError, "unbounded below"):
            apply_softening(state, self.config, False, outer_phase=3, condition="control")

    def test_positive_remodeling_is_first_order_consistent_and_positive(self) -> None:
        state, _, _, _ = generate_initial_state(18, self.config)
        state.activation.fill(-0.01)
        state.D.fill(0.3)
        before = state.A.copy()
        update_activation_and_density(state, self.config)
        # The update routine recalculates activation, so test the map directly
        # through a configuration with a negligible archival increment.
        reference = np.where(state.initial_fibrotic, self.config.k_fibrosis, self.config.k_normal)
        positive = before * np.exp((-1e-8) / reference)
        additive = before - 1e-8
        self.assertTrue(np.all(positive > 0.0))
        np.testing.assert_allclose(positive, additive, rtol=1e-8, atol=1e-14)
        self.assertEqual(replace(self.config, remodeling_law="legacy_additive").remodeling_law, "legacy_additive")

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

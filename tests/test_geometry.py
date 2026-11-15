"""Independent regression tests; run with python -m unittest discover -s tests."""
import json
import math
import random
import unittest

from explorative_math.geometry import (
    angle_measurement, count_geometry, iter_geometry_runs, shell_fraction,
    unit_vector, validate_geometry_config,
)


class GeometryTests(unittest.TestCase):
    def runs(self, **kwargs):
        return list(iter_geometry_runs({"dimensions": [2, 8], "samples": 12,
                                        "steps": 5, **kwargs}))

    def test_exact_counts_and_blank_slots(self):
        config = {"dimensions": [2, 8], "samples": 12, "steps": 5}
        runs = list(iter_geometry_runs(config))
        counts = count_geometry(config)
        self.assertEqual(len(runs), counts["runs"])
        self.assertEqual(sum(len(r["outputs"]) for r in runs), counts["outputs"])
        self.assertEqual(len({r["id"] for r in runs}), len(runs))
        for run in runs:
            self.assertEqual(run["status"], "OBSERVED")
            for field in ("human_interest", "ai_interest", "compressibility", "divergence_note"):
                self.assertIsNone(run[field])
            for output in run["outputs"]:
                self.assertIn(output["status"], ("OBSERVED", "DERIVED", "PROPOSED"))
        json.dumps(runs, allow_nan=False)

    def test_repeatable_not_filtered(self):
        first = self.runs(seed=13)
        self.assertEqual(first, self.runs(seed=13))
        self.assertNotEqual(first, self.runs(seed=14))
        orth = first[0]
        self.assertEqual([x["index"] for x in orth["outputs"][:-1]], list(range(12)))

    def test_sphere_vectors_and_angles(self):
        rng = random.Random(7)
        for n in (1, 2, 8, 64):
            vector = unit_vector(rng, n)
            self.assertAlmostEqual(math.fsum(x*x for x in vector), 1.0, places=12)
        self.assertAlmostEqual(angle_measurement([1, 0], [0, 1])["angle_degrees"], 90)
        self.assertAlmostEqual(angle_measurement([1, 0], [-1, 0])["angle_degrees"], 180)
        self.assertEqual(angle_measurement([1.0000000000000002], [1])["angle_degrees"], 0)
        with self.assertRaises(ValueError):
            angle_measurement([1], [1, 2])

    def test_near_orthogonality_squared_dot(self):
        runs = list(iter_geometry_runs({"experiment": "orthogonality", "dimensions": [2, 64],
                                        "samples": 1500, "seed": 17}))
        low, high = [r["outputs"][-1]["mean_squared_dot"] for r in runs]
        self.assertAlmostEqual(low, 0.5, delta=0.04)
        self.assertAlmostEqual(high, 1/64, delta=0.003)
        self.assertLess(high, low)

    def test_shell_formula_and_uniform_volume(self):
        self.assertEqual(shell_fraction(8, 0), 0)
        self.assertEqual(shell_fraction(8, 1), 1)
        self.assertAlmostEqual(shell_fraction(2, 0.1), 0.19)
        self.assertGreater(shell_fraction(64, 0.1), shell_fraction(2, 0.1))
        run = list(iter_geometry_runs({"experiment": "shell", "dimensions": [8],
                                       "samples": 2000, "seed": 21}))[0]
        for item in run["outputs"][:-1]:
            self.assertAlmostEqual(item["radius"]**8, item["uniform"], places=12)
            self.assertAlmostEqual(math.sqrt(math.fsum(x*x for x in item["point"])),
                                   item["radius"], places=12)
        summary = run["outputs"][-1]
        self.assertAlmostEqual(summary["empirical_fraction"], summary["exact_fraction"], delta=0.04)

    def test_superposition_records_every_pair(self):
        run = self.runs(experiment="superposition", dimensions=[3], vectors=7)[0]
        directions = [x for x in run["outputs"] if x["kind"] == "direction"]
        pairs = [x for x in run["outputs"] if x["kind"] == "pair_angle"]
        self.assertEqual(len(directions), 7)
        self.assertEqual(len(pairs), 21)
        self.assertEqual({(p["left_index"], p["right_index"]) for p in pairs},
                         {(i,j) for i in range(7) for j in range(i+1,7)})
        self.assertEqual(run["outputs"][-1]["min_angle_degrees"],
                         min(p["angle_degrees"] for p in pairs))

    def test_paths_include_initial_and_every_transition(self):
        walk, mapping = self.runs(experiment="trajectories", dimensions=[3])
        for run in (walk, mapping):
            self.assertEqual([p["step"] for p in run["outputs"]], list(range(6)))
        for before, after in zip(walk["outputs"], walk["outputs"][1:]):
            self.assertEqual(after["state"], [x+d for x,d in zip(before["state"], after["increment"])])
        for before, after in zip(mapping["outputs"], mapping["outputs"][1:]):
            self.assertEqual(after["state"], [mapping["params"]["scale"]*x+b
                for x,b in zip(before["state"], mapping["params"]["translation"])])

    def test_zero_step_paths_are_kept(self):
        runs = self.runs(experiment="trajectories", dimensions=[1], steps=0)
        self.assertEqual([len(r["outputs"]) for r in runs], [1, 1])

    def test_invalid_bounds_rejected(self):
        bad_configs = [
            {"dimensions": []}, {"dimensions": [True]}, {"dimensions": [2,2]},
            {"dimensions": [0]}, {"samples": 0}, {"steps": -1}, {"seed": True},
            {"vectors": 2, "dimensions": [2]}, {"shell_thickness": float("nan")},
            {"shell_thickness": -0.1}, {"shell_thickness": 2}, {"map_scale": float("inf")},
            {"map_scale": 1.1}, {"experiment": "unknown"},
        ]
        for config in bad_configs:
            with self.subTest(config=config), self.assertRaises(ValueError):
                validate_geometry_config(config)


if __name__ == "__main__":
    unittest.main()

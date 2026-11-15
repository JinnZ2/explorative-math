"""Projection regression tests; stdlib-only and independent of the real catalog."""
import copy
import json
import math
import unittest
from unittest.mock import patch

from explorative_math.projections import (
    DISCLAIMER, GOLDEN_ANGLE_DEGREES, PHI, _compare, _measure,
    count_projection, default_catalog, iter_projection_runs,
    validate_projection_config,
)
from explorative_math.storage import check_run_record


def catalog_fixture():
    """The four requested signatures, with deliberately synthetic source URLs."""
    signatures = [("golden-phyllotaxis", "Golden phyllotaxis", GOLDEN_ANGLE_DEGREES),
                  ("icosahedral-capsids", "Icosahedral capsids", 72.0),
                  ("quasicrystal-fivefold", "Quasicrystal fivefold", 72.0),
                  ("hexagonal-packing", "Hexagonal packing", 60.0)]
    return [{"id": identifier, "label": label, "target_angle_degrees": target,
             "source_url": "https://example.invalid/" + identifier,
             "status": "DERIVED", "interpretation": "Angular signature only; not all physical angles."}
            for identifier, label, target in signatures]


class ProjectionTests(unittest.TestCase):
    def run_projection(self, catalog=None, **config):
        records = list(iter_projection_runs(config, catalog_fixture() if catalog is None else catalog))
        self.assertEqual(len(records), 1)
        check_run_record(records[0])
        return records[0]

    def test_defaults_and_exact_archive_counts(self):
        cfg = validate_projection_config({})
        self.assertEqual(cfg, {"dimension": 5, "shape": "crosspolytope", "frames": 12,
                               "seed": 0, "composition": "phi", "tolerance_degrees": 1.0,
                               "max_vertices": 4096})
        self.assertEqual(count_projection({}), {"runs": 1, "outputs": 13})
        run = self.run_projection()
        self.assertEqual(len(run["outputs"]), 13)
        self.assertEqual(run["family"], "projections")
        self.assertEqual(run["id"], run["run_id"])
        self.assertTrue(run["id"].startswith("projection-"))
        self.assertEqual(run["status"], "OBSERVED")
        self.assertEqual(run["outputs"][0]["kind"], "geometry_definition")
        for field in ("human_interest", "ai_interest", "compressibility", "divergence_note"):
            self.assertIsNone(run[field])
        json.dumps(run, allow_nan=False)

    def test_all_shape_counts_and_unit_circumradii(self):
        for dimension in (2, 3, 5, 8):
            for shape, expected in (("crosspolytope", 2 * dimension),
                                    ("hypercube", 2 ** dimension), ("simplex", dimension + 1)):
                with self.subTest(dimension=dimension, shape=shape):
                    run = self.run_projection(dimension=dimension, shape=shape, frames=1,
                                              composition="none")
                    vertices = run["outputs"][0]["original_vertices"]
                    self.assertEqual(len(vertices), expected)
                    self.assertEqual(run["outputs"][0]["generated_vertex_count"], expected)
                    self.assertEqual(len({tuple(v) for v in vertices}), expected)
                    for vertex in vertices:
                        self.assertEqual(len(vertex), dimension)
                        self.assertAlmostEqual(math.fsum(x * x for x in vertex), 1.0, places=12)
                    for axis in range(dimension):
                        self.assertAlmostEqual(math.fsum(v[axis] for v in vertices), 0.0, places=12)

    def test_centered_regular_simplex_pairwise_dots(self):
        for dimension in (2, 3, 5, 12):
            vertices = self.run_projection(dimension=dimension, shape="simplex", frames=1)["outputs"][0]["original_vertices"]
            for i, left in enumerate(vertices):
                for j, right in enumerate(vertices):
                    dot = math.fsum(x * y for x, y in zip(left, right))
                    self.assertAlmostEqual(dot, 1.0 if i == j else -1.0 / dimension, places=12)
                    if i != j:
                        distance2 = math.fsum((x - y) ** 2 for x, y in zip(left, right))
                        self.assertAlmostEqual(distance2, 2.0 * (dimension + 1) / dimension, places=12)

    def test_givens_rotations_preserve_norms_and_inner_products(self):
        for shape in ("crosspolytope", "hypercube", "simplex"):
            run = self.run_projection(dimension=4, shape=shape, frames=7, seed=27)
            original = run["outputs"][0]["original_vertices"]
            for frame in run["outputs"][1:]:
                self.assertTrue(frame["rotation_metadata"]["sampled_not_exhaustive"])
                rotated = frame["rotated_vertices"]
                self.assertEqual(len(rotated), len(original))
                for i, vertex in enumerate(rotated):
                    self.assertAlmostEqual(math.fsum(x * x for x in vertex), 1.0, places=12)
                    self.assertEqual(frame["projected_points"][i], vertex[:2])
                    for j in range(len(rotated)):
                        self.assertAlmostEqual(math.fsum(x * y for x, y in zip(vertex, rotated[j])),
                                               math.fsum(x * y for x, y in zip(original[i], original[j])), places=12)
                for vertex in frame["shapes"][1]["rotated_vertices"]:
                    self.assertAlmostEqual(math.fsum(x * x for x in vertex), PHI ** 2, places=12)

    def test_unrotated_first_frame_and_even_angle_schedule(self):
        run = self.run_projection(dimension=5, frames=6, seed=-5)
        frames = run["outputs"][1:]
        self.assertEqual(frames[0]["rotated_vertices"], run["outputs"][0]["original_vertices"])
        for rotation in frames[0]["rotation_metadata"]["rotations"]:
            self.assertEqual(rotation["angle_degrees"], 0.0)
        for frame_index, frame in enumerate(frames):
            metadata = frame["rotation_metadata"]
            self.assertEqual(frame["frame_index"], frame_index)
            self.assertEqual(metadata["progress"], frame_index / 6)
            touched = set()
            for rotation in metadata["rotations"]:
                i, j = rotation["plane"]
                self.assertLess(i, j)
                touched.update((i, j))
                self.assertAlmostEqual(rotation["angle_degrees"], frame_index / 6 * rotation["sweep_degrees"])
            self.assertEqual(touched, set(range(5)))
        self.assertNotEqual(frames[-1]["rotated_vertices"], run["outputs"][0]["original_vertices"])

    def test_reproducibility_and_catalog_in_identity(self):
        first = self.run_projection(frames=4, seed=19)
        self.assertEqual(first, self.run_projection(frames=4, seed=19))
        self.assertNotEqual(first["id"], self.run_projection(frames=4, seed=20)["id"])
        self.assertNotEqual(first["outputs"][2]["rotated_vertices"],
                            self.run_projection(frames=4, seed=20)["outputs"][2]["rotated_vertices"])
        catalog = catalog_fixture()
        catalog[0]["interpretation"] += " Changed metadata."
        changed = self.run_projection(catalog=catalog, frames=4, seed=19)
        self.assertNotEqual(first["id"], changed["id"])
        self.assertEqual(first["outputs"][2]["rotated_vertices"], changed["outputs"][2]["rotated_vertices"])
        self.assertEqual(self.run_projection(frames=2, tolerance_degrees=1)["id"],
                         self.run_projection(frames=2, tolerance_degrees=1.0)["id"])

    def test_every_frame_and_every_target_comparison_is_retained(self):
        run = self.run_projection(frames=9, dimension=4)
        self.assertEqual(count_projection({"frames": 9}), {"runs": 1, "outputs": 10})
        self.assertEqual([f["frame_index"] for f in run["outputs"][1:]], list(range(9)))
        comparisons = [c for frame in run["outputs"][1:] for c in frame["catalog_comparisons"]]
        self.assertEqual(len(comparisons), 36)
        self.assertEqual({(c["frame_index"], c["catalog_index"]) for c in comparisons},
                         {(frame, target) for frame in range(9) for target in range(4)})
        for frame in run["outputs"][1:]:
            self.assertEqual(frame["kind"], "projection_frame")
            self.assertEqual(frame["status"], "OBSERVED")
            self.assertEqual(frame["composed_measurements"]["status"], "DERIVED")
            gaps = frame["angular_gaps"]
            self.assertAlmostEqual(math.fsum(g["gap_degrees"] for g in gaps), 360.0, places=11)
            for comparison in frame["catalog_comparisons"]:
                target = comparison["target_angle_degrees"]
                errors = [min(abs(g["gap_degrees"] - target), 360 - abs(g["gap_degrees"] - target))
                          for g in gaps]
                self.assertEqual(comparison["measurement_errors_degrees"], errors)
                self.assertEqual(comparison["absolute_error_degrees"], min(errors))
                index = comparison["measurement_index"]
                self.assertEqual(comparison["observed_nearest_degrees"], gaps[index]["gap_degrees"])
                self.assertEqual(comparison["nearest_measurement"], gaps[index])
                self.assertEqual(comparison["nearest_measurement_indices"],
                                 [i for i, e in enumerate(errors) if e == min(errors)])
                self.assertEqual(comparison["matching_measurement_indices"],
                                 [i for i, e in enumerate(errors) if e <= 1])
                self.assertEqual(comparison["matched"], min(errors) <= 1)
                self.assertEqual(comparison["status"], "DERIVED")
                self.assertEqual(comparison["disclaimer"], DISCLAIMER)

    def test_both_matches_and_misses_retained_without_selection(self):
        catalog = catalog_fixture()
        catalog[0]["target_angle_degrees"] = 90.0
        catalog[0]["source_title"] = "Synthetic 90-degree reference"
        catalog[0]["related_source_url"] = "https://example.invalid/related"
        frame = self.run_projection(catalog=catalog, dimension=2, frames=1,
                                    composition="none")["outputs"][1]
        comparisons = frame["catalog_comparisons"]
        self.assertEqual(len(comparisons), 4)
        self.assertEqual([c["result"] for c in comparisons], ["match", "miss", "miss", "miss"])
        self.assertEqual(comparisons[0]["matching_measurement_indices"], [0, 1, 2, 3])
        for index, comparison in enumerate(comparisons):
            self.assertEqual(comparison["catalog_record"], catalog[index])
            self.assertEqual(comparison["catalog_record"]["status"], "DERIVED")

    def test_phi_scaling_alone_does_not_change_angles(self):
        run = self.run_projection(dimension=3, frames=5, composition="phi")
        for frame in run["outputs"][1:]:
            before, scaled = frame["precomposition_measurements"], frame["pure_scale_measurements"]
            self.assertEqual(before["polar_angles_degrees"], scaled["polar_angles_degrees"])
            self.assertEqual(before["angular_gaps"], scaled["angular_gaps"])
            for p, q in zip(before["polar"], scaled["polar"]):
                self.assertAlmostEqual(q["radius"], PHI * p["radius"], places=12)
            self.assertNotEqual(before["polar_angles_degrees"],
                                frame["composed_measurements"]["polar_angles_degrees"])
            for original, scaled_point in zip(frame["projected_points"], frame["pure_scale_points"]):
                self.assertEqual(scaled_point, [PHI * x for x in original])
        first = run["outputs"][1]
        before_angles = first["projected_polar_angles_degrees"]
        copy_angles = first["composed_measurements"]["polar_angles_degrees"][6:]
        for before, after in zip(before_angles, copy_angles):
            if before is None:
                self.assertIsNone(after)
            else:
                self.assertAlmostEqual(after, (before + GOLDEN_ANGLE_DEGREES) % 360, places=12)

    def test_compositions_preserve_both_complete_shapes(self):
        for composition in ("none", "phi", "second"):
            run = self.run_projection(dimension=3, shape="simplex", frames=3, composition=composition)
            copies = 1 if composition == "none" else 2
            definition = run["outputs"][0]
            self.assertEqual(len(definition["shapes"]), copies)
            self.assertEqual(definition["generated_vertex_count"], 4 * copies)
            for shape in definition["shapes"]:
                self.assertEqual(shape["original_vertices"], definition["original_vertices"])
            for frame in run["outputs"][1:]:
                self.assertEqual(len(frame["shapes"]), copies)
                self.assertEqual(len(frame["composed_points"]), 4 * copies)
                self.assertEqual(len(frame["composed_measurements"]["polar"]), 4 * copies)
                self.assertEqual(frame["composed_point_origins"],
                                 [{"shape_index": shape, "vertex_index": vertex}
                                  for shape in range(copies) for vertex in range(4)])
                if composition == "second":
                    for vertex in frame["shapes"][1]["rotated_vertices"]:
                        self.assertAlmostEqual(math.fsum(x * x for x in vertex), 1, places=12)

    def test_degenerate_zero_points_are_retained_with_null_angles(self):
        run = self.run_projection(dimension=5, frames=1)
        frame = run["outputs"][1]
        before = frame["precomposition_measurements"]
        self.assertEqual(len(frame["projected_points"]), 10)
        self.assertEqual(len(before["polar"]), 10)
        self.assertEqual(before["polar_angles_degrees"].count(None), 6)
        self.assertEqual(before["degeneracies"]["undefined_point_indices"], list(range(4, 10)))
        self.assertEqual(len(before["angular_gaps"]), 4)
        after = frame["composed_measurements"]
        self.assertEqual(len(frame["composed_points"]), 20)
        self.assertEqual(after["polar_angles_degrees"].count(None), 12)
        self.assertEqual(len(after["angular_gaps"]), 8)
        for p in after["polar"]:
            if p["angle_degrees"] is None:
                self.assertTrue(p["degenerate"])
                self.assertEqual(p["degeneracy_reason"], "zero_projection")

    def test_duplicate_directions_and_zero_gaps_are_not_deduplicated(self):
        frame = self.run_projection(dimension=4, shape="hypercube", frames=1,
                                    composition="none")["outputs"][1]
        measures = frame["precomposition_measurements"]
        self.assertEqual(len(measures["polar"]), 16)
        self.assertEqual(len(measures["angular_gaps"]), 16)
        self.assertEqual(len(measures["degeneracies"]["duplicate_direction_groups"]), 4)
        self.assertEqual(sum(g["gap_degrees"] == 0 for g in measures["angular_gaps"]), 12)
        self.assertEqual(sorted(measures["sorted_point_indices"]), list(range(16)))
        for gap in measures["angular_gaps"]:
            if gap["gap_degrees"] == 0:
                self.assertEqual(gap["degeneracy_reason"], "duplicate_direction")

    def test_empty_gap_lists_log_all_misses_and_single_direction_wraps(self):
        origins = [{"shape_index": 0, "vertex_index": 0}]
        empty = _measure([[0.0, 0.0]], origins)
        self.assertEqual(empty["angular_gaps"], [])
        self.assertEqual(empty["polar_angles_degrees"], [None])
        comparisons = _compare(empty["angular_gaps"], catalog_fixture(), 1.0, 3)
        self.assertEqual(len(comparisons), 4)
        for comparison in comparisons:
            self.assertFalse(comparison["matched"])
            self.assertEqual(comparison["result"], "miss")
            self.assertEqual(comparison["reason"], "no_defined_angular_gaps")
            self.assertIsNone(comparison["measurement_index"])
            self.assertIsNone(comparison["absolute_error_degrees"])
            self.assertEqual(comparison["measurement_errors_degrees"], [])
        single = _measure([[1.0, 0.0]], origins)
        self.assertEqual(single["angular_gaps"][0]["gap_degrees"], 360.0)
        same = _measure([[1.0, 0.0], [2.0, 0.0]], origins * 2)
        self.assertEqual([g["gap_degrees"] for g in same["angular_gaps"]], [0, 360])

    def test_circular_distance_and_tolerance_inclusivity(self):
        target = catalog_fixture()[0]
        target["target_angle_degrees"] = 10.0
        gaps = [{"measurement_index": 0, "gap_degrees": 350.0}]
        exact = _compare(gaps, [target], 20.0, 0)[0]
        self.assertEqual(exact["absolute_error_degrees"], 20.0)
        self.assertTrue(exact["matched"])
        self.assertFalse(_compare(gaps, [target], 19.99, 0)[0]["matched"])
        for tolerance, expected in ((0, False), (180, True)):
            comparison = self.run_projection(frames=1, composition="none", tolerance_degrees=tolerance)["outputs"][1]["catalog_comparisons"][0]
            self.assertEqual(comparison["matched"], expected)

    def test_catalog_is_copied_and_explicit_empty_catalog_is_supported(self):
        catalog = catalog_fixture()
        original = copy.deepcopy(catalog)
        run = self.run_projection(catalog=catalog, frames=1)
        self.assertEqual(catalog, original)
        catalog[0]["label"] = "Mutated"
        self.assertEqual(run["params"]["catalog"], original)
        run = self.run_projection(catalog=[], frames=2)
        self.assertTrue(all(f["catalog_comparisons"] == [] for f in run["outputs"][1:]))

    def test_packaged_catalog_loader_and_missing_resource(self):
        for body in (catalog_fixture(), {"entries": catalog_fixture()}):
            with patch("explorative_math.projections.resources.files") as files:
                files.return_value.joinpath.return_value.read_text.return_value = json.dumps(body)
                self.assertEqual(default_catalog(), catalog_fixture())
                files.assert_called_once_with("explorative_math")
                files.return_value.joinpath.assert_called_once_with("nature_catalog.json")
        with patch("explorative_math.projections.resources.files") as files:
            files.return_value.joinpath.return_value.read_text.side_effect = FileNotFoundError("missing")
            with self.assertRaisesRegex(ValueError, "nature_catalog.json"):
                list(iter_projection_runs({"frames": 1}))
        with patch("explorative_math.projections.resources.files", side_effect=AssertionError("must not read")):
            self.run_projection(catalog=catalog_fixture(), frames=1)
            self.assertEqual(count_projection({}), {"runs": 1, "outputs": 13})

    def test_generated_vertex_limit_refuses_not_truncates(self):
        with self.assertRaisesRegex(ValueError, "exceeds max_vertices=31; no vertices truncated"):
            validate_projection_config({"shape": "hypercube", "dimension": 5,
                                        "composition": "none", "max_vertices": 31})
        with self.assertRaisesRegex(ValueError, "exceeds max_vertices=19"):
            validate_projection_config({"dimension": 5, "max_vertices": 19})
        self.assertEqual(self.run_projection(dimension=5, max_vertices=20, frames=1)["params"]["generated_vertex_count"], 20)
        self.assertEqual(count_projection({"shape": "hypercube", "dimension": 12,
                                           "composition": "none"}), {"runs": 1, "outputs": 13})
        with self.assertRaises(ValueError):
            count_projection({"shape": "hypercube", "dimension": 12})
        with self.assertRaises(ValueError):
            count_projection({"shape": "hypercube", "dimension": 10 ** 100})

    def test_invalid_bounds_rejected(self):
        invalid = [None, [], {"dimension": 1}, {"dimension": 0}, {"dimension": True},
                   {"dimension": 2.0}, {"frames": 0}, {"frames": -1}, {"frames": True},
                   {"frames": 1.0}, {"seed": True}, {"seed": "0"}, {"seed": 0.0},
                   {"shape": "sphere"}, {"shape": None}, {"composition": "golden"},
                   {"composition": None}, {"max_vertices": 0}, {"max_vertices": True},
                   {"max_vertices": 4096.0}, {"tolerance_degrees": -0.1},
                   {"tolerance_degrees": 180.1}, {"tolerance_degrees": float("nan")},
                   {"tolerance_degrees": float("inf")}, {"tolerance_degrees": True},
                   {"tolerance_degrees": "1"}, {"tolerance_degrees": 10 ** 1000}]
        for config in invalid:
            with self.subTest(config=config), self.assertRaises(ValueError):
                validate_projection_config(config)

    def test_invalid_catalog_rejected(self):
        for field, value in (("target_angle_degrees", float("nan")),
                             ("target_angle_degrees", -1), ("target_angle_degrees", 181),
                             ("status", "OBSERVED"), ("id", ""), ("source_url", None),
                             ("interpretation", "")):
            catalog = catalog_fixture()
            catalog[0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                self.run_projection(catalog=catalog, frames=1)
        for catalog in ("bad", [None], {"unknown": []}, {"entries": [], "records": []}):
            with self.subTest(catalog=catalog), self.assertRaises(ValueError):
                list(iter_projection_runs({"frames": 1}, catalog=catalog))


if __name__ == "__main__":
    unittest.main()

"""cross -- the same inputs through both builds, logged field by field.

CC0. Co-authored: JinnZ2, Manus AI, Claude. Stdlib only.

Build A = Manus      (explorative_math/)
Build B = Claude Code (toolkit/, experiments/)

Each comparison feeds one input to the matching function in each build and
logs every field as AGREE / DIFFER / ONLY_IN_A / ONLY_IN_B. Numbers carry
absolute and relative difference. AGREE means exactly equal; there is no
tolerance and nothing is filtered, so a 1e-16 difference is a DIFFER with
its size recorded. Reading a DIFFER as "same up to rounding" is a reading,
left to the observer slots, which are null.

Some private functions of each build are called directly (for example the
Manus build's _vertices, _rotate, _measure). This reads each build as it is;
it changes neither.

  python3 -m compare            # writes compare/runs/<name>.jsonl
"""
import json
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "experiments"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from explorative_math import geometry as A_geo  # noqa: E402
from explorative_math import projections as A_proj  # noqa: E402
from toolkit import composers as B_comp  # noqa: E402
from toolkit import primitives as B_prim  # noqa: E402
from toolkit import probes as B_probe  # noqa: E402
import near_orthogonality as B_orth  # noqa: E402
import shell_concentration as B_shell  # noqa: E402

OUT = ROOT / "compare" / "runs"
AGREE, DIFFER, ONLY_A, ONLY_B = "AGREE", "DIFFER", "ONLY_IN_A", "ONLY_IN_B"
VERDICTS = (AGREE, DIFFER, ONLY_A, ONLY_B)

# union of both builds' run schemas (see docs/two_builds.md)
SCHEMA_A = ("id", "params", "outputs", "human_interest", "ai_interest",
            "compressibility", "divergence_note", "status")
SCHEMA_B = ("params", "outputs", "angles", "human_interest", "ai_interest",
            "compressibility", "nature_match", "divergence_note")
SCHEMA_UNION = tuple(dict.fromkeys(SCHEMA_A + SCHEMA_B))


# ---- field comparison -------------------------------------------------------

def _num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def field(name, a, b, present_a=True, present_b=True, note=None):
    rec = {"field": name}
    if not present_a and not present_b:
        raise ValueError("a field must exist in at least one build")
    if not present_b:
        rec.update(verdict=ONLY_A, a=a, b=None)
    elif not present_a:
        rec.update(verdict=ONLY_B, a=None, b=b)
    elif _num(a) and _num(b):
        diff = abs(a - b)
        scale = max(abs(a), abs(b))
        rec.update(verdict=AGREE if a == b else DIFFER, a=a, b=b, abs_diff=diff,
                   rel_diff=(diff / scale) if scale else 0.0)
    else:
        rec.update(verdict=AGREE if a == b else DIFFER, a=a, b=b)
    if note:
        rec["note"] = note
    return rec


def vector_fields(prefix, va, vb):
    """One field per element; lengths that differ give ONLY_IN_A / ONLY_IN_B tails."""
    out = []
    for i in range(max(len(va), len(vb))):
        out.append(field(f"{prefix}[{i}]", va[i] if i < len(va) else None,
                         vb[i] if i < len(vb) else None, i < len(va), i < len(vb)))
    return out


def tally(fields):
    t = {v: 0 for v in VERDICTS}
    for f in fields:
        t[f["verdict"]] += 1
    nums = [f["abs_diff"] for f in fields if "abs_diff" in f]
    t["max_abs_diff"] = max(nums) if nums else None
    return t


def record(name, params, fields, extra=None):
    """extra: keys placed in outputs beside fields and tally. `bounds` holds the
    instrument bound each probe reading was taken at (docs/resolution_leak.md);
    it is not a compared field and does not enter the tally."""
    rec = {k: None for k in SCHEMA_UNION}
    rec.update(id=f"compare-{name}", params=params,
               outputs={"fields": fields, "tally": tally(fields), **(extra or {})},
               status="OBSERVED")
    return rec


# ---- comparisons --------------------------------------------------------------

def cmp_constants():
    a = A_proj.GOLDEN_ANGLE_DEGREES
    b = 360.0 * (1.0 - 1.0 / B_prim.PHI)
    fields = [field("phi", A_proj.PHI, B_prim.PHI),
              field("golden_angle_deg", a, b,
                    note="A: 360/phi^2; B: 360*(1-1/phi); equal in exact arithmetic")]
    return [record("constants", {"what": "phi and golden angle"}, fields)]


def cmp_shell_formula():
    recs = []
    for d in (1, 2, 3, 5, 10, 20, 50, 100, 500):
        for t in (0.01, 0.05, 0.1, 0.5):
            f = [field("shell_fraction", A_geo.shell_fraction(d, t), B_shell.exact(d, t),
                       note="A: -expm1(n*log1p(-t)); B: 1-(1-t)**n")]
            recs.append(record(f"shell_formula_d{d}_t{t}", {"dimension": d, "thickness": t}, f))
    return recs


def cmp_seed_alignment():
    """Does the same seed give the same random stream in both builds?"""
    recs = []
    for d in (2, 8, 32):
        cfg = A_geo.validate_geometry_config({"dimensions": [d], "samples": 4, "seed": 0,
                                              "experiment": "orthogonality"})
        run_a = A_geo._orthogonality(d, cfg)
        pairs_a = [o for o in run_a["outputs"] if o["kind"] == "pair"]
        rng_b = random.Random(0)
        pairs_b = []
        for _ in range(4):
            u, v = B_prim.random_unit(rng_b, d), B_prim.random_unit(rng_b, d)
            pairs_b.append((u, v))
        fields = []
        for k in range(4):
            fields += vector_fields(f"pair{k}.left", pairs_a[k]["left"], pairs_b[k][0])
        recs.append(record(f"seed_alignment_d{d}",
                           {"dimension": d, "seed": 0, "pairs": 4,
                            "A_stream": "random.Random(sha256(canonical params))",
                            "B_stream": "random.Random(seed)"}, fields))
    return recs


def cmp_shared_stream_sampler():
    """Same explicit Random object state into both samplers: isolates the sampler."""
    recs = []
    for d in (2, 8, 32, 128):
        for s in (0, 1):
            va = A_geo.unit_vector(random.Random(s), d)
            vb = B_prim.random_unit(random.Random(s), d)
            recs.append(record(f"shared_stream_unit_vector_d{d}_s{s}",
                               {"dimension": d, "explicit_seed": s,
                                "A_norm": "math.fsum", "B_norm": "sum"},
                               vector_fields("unit_vector", va, vb)))
    return recs


def cmp_orthogonality_stats():
    """Each build's own experiment at its own defaults where dims overlap: statistics only."""
    recs = []
    per_b, _ = B_orth.run(seed=0)
    for d in (2, 20, 50, 100, 200):
        cfg = A_geo.validate_geometry_config({"dimensions": [d], "samples": B_orth.PAIRS,
                                              "seed": 0, "experiment": "orthogonality"})
        summ = [o for o in A_geo._orthogonality(d, cfg)["outputs"] if o["kind"] == "summary"][0]
        b = per_b[str(d)]
        fields = [field("mean_angle_degrees", summ["mean_angle_degrees"], b["mean_angle"]),
                  field("mean_squared_dot", summ["mean_squared_dot"], None, True, False),
                  field("sd_cos", None, b["sd_cos"], False, True),
                  field("expected_squared_dot", summ["expected_squared_dot"], None, True, False),
                  field("frac_within_band", None, b["frac_within_band"], False, True)]
        recs.append(record(f"orthogonality_stats_d{d}",
                           {"dimension": d, "pairs": B_orth.PAIRS, "seed": 0,
                            "note": "different random streams; statistics, not samples"}, fields))
    return recs


def _circumradius(vs):
    return max(B_prim.norm(v) for v in vs)


def _sorted_normalized_distances(vs):
    r = _circumradius(vs)
    return sorted(round(math.dist(a, b) / r, 12) for i, a in enumerate(vs) for b in vs[i + 1:])


SHAPE_PAIRS = {"simplex": ("simplex", B_prim.simplex),
               "hypercube": ("hypercube", B_prim.hypercube),
               "crosspolytope": ("crosspolytope", B_prim.cross_polytope)}


def cmp_polytopes():
    recs = []
    for name, (a_name, b_fn) in SHAPE_PAIRS.items():
        for n in (3, 4, 5):
            va, vb = A_proj._vertices(a_name, n), b_fn(n)
            fields = [field("vertex_count", len(va), len(vb)),
                      field("ambient_dimension", len(va[0]), len(vb[0])),
                      field("circumradius", _circumradius(va), _circumradius(vb)),
                      field("vertex_list_identical", va, vb),
                      field("normalized_pairwise_distances", _sorted_normalized_distances(va),
                            _sorted_normalized_distances(vb),
                            note="rounded to 12 places; invariant under rotation and scale")]
            recs.append(record(f"polytope_{name}_n{n}", {"shape": name, "dimension": n}, fields))
    a_only = [field(f"shape:{s}", True, None, True, False) for s in ()]
    b_only = [field(f"shape:{s}", None, True, False, True) for s in ("cell24", "icosahedron", "regular_polygon")]
    recs.append(record("polytope_catalogue", {"what": "shapes each build can generate"}, a_only + b_only))
    return recs


def _shared_rotations(n, seed):
    rng = random.Random(seed)
    planes = B_prim.all_planes(n)
    return [(pl, rng.uniform(0, 360)) for pl in planes]


def cmp_rotate_project():
    """Shared explicit rotations, same vertex set (A's), applied by each build's rotator."""
    recs = []
    for name in ("simplex", "hypercube", "crosspolytope"):
        n = 4
        verts = A_proj._vertices(name, n)
        rots = _shared_rotations(n, seed=7)
        pa = [v[:2] for v in A_proj._rotate(verts, [{"plane": list(pl), "angle_degrees": deg}
                                                   for pl, deg in rots])]
        pb = B_comp.rotate_project(verts, [(pl, math.radians(deg)) for pl, deg in rots])
        fields = []
        for i, (x, y) in enumerate(zip(pa, pb)):
            fields += vector_fields(f"point{i}", x, y)
        recs.append(record(f"rotate_project_{name}_n{n}",
                           {"shape": name, "dimension": n, "vertices": "A's _vertices for both",
                            "rotations_deg": [[list(pl), deg] for pl, deg in rots],
                            "A_angle_unit": "degrees", "B_angle_unit": "radians"}, fields))
    return recs


def cmp_angle_measure():
    """Polar angles of one projected point set: A about the origin, B about the centroid."""
    recs = []
    for name in ("simplex", "hypercube", "crosspolytope"):
        n = 4
        verts = A_proj._vertices(name, n)
        rots = _shared_rotations(n, seed=7)
        pts = B_comp.rotate_project(verts, [(pl, math.radians(deg)) for pl, deg in rots])
        meas = A_proj._measure(pts, [{"shape_index": 0, "vertex_index": i} for i in range(len(pts))])
        a_sorted = sorted(x for x in meas["polar_angles_degrees"] if x is not None)
        b_sorted = B_probe.polar_angles(pts)
        fields = vector_fields("sorted_polar_angle_deg", a_sorted, b_sorted)
        fields.append(field("angular_gaps", [g["gap_degrees"] for g in meas["angular_gaps"]],
                            None, True, False))
        fields.append(field("rotational_symmetry_k", None, B_probe.rotational_symmetry(pts),
                            False, True))
        fields.append(field("mirror_axes_deg", None, B_probe.mirror_symmetry(pts), False, True))
        fields.append(field("catalog_comparisons", "A compares gaps to its nature catalog",
                            None, True, False))
        recs.append(record(f"angle_measure_{name}_n{n}",
                           {"shape": name, "dimension": n,
                            "A_origin": "coordinate origin", "B_origin": "centroid"}, fields,
                           extra={"bounds": {
                               "rotational_symmetry_k": {
                                   "b": B_probe.rotational_symmetry_bound(pts)},
                               "mirror_axes_deg": {
                                   "b": B_probe.mirror_symmetry_bound(pts)}}}))
    return recs


def cmp_phi_scale_vs_phi_scale_rotate():
    """A's own projection pipeline (composition=phi); B's phi_stack on A's projected points.

    A's frames come from A's seeded tour. Each frame's projected points (before
    composition) are handed to B's phi_stack, so both compositions start from
    identical 2D points and differ only in what "phi composition" means.

    "phi composition" is a welded term: one label, two operations. In this
    comparison's output the two are named phi_scale (B: scale only) and
    phi_scale_rotate (A: scale and rotate by the golden angle). Neither
    build's own code or naming is changed (docs/two_builds.md)."""
    recs = []
    for name in ("simplex", "crosspolytope"):
        n = 4
        run = next(A_proj.iter_projection_runs(
            {"dimension": n, "shape": name, "frames": 4, "seed": 0, "composition": "phi"}))
        for frame in [o for o in run["outputs"] if o["kind"] == "projection_frame"]:
            pts = frame["projected_points"]
            a_union = frame["composed_points"]
            b_union = [r["point"] for r in B_comp.phi_stack(pts, (0, 1))]
            a_ang = B_probe.polar_angles(a_union)
            b_ang = B_probe.polar_angles(b_union)
            fields = [field("point_count", len(a_union), len(b_union)),
                      field("distinct_directions", len({round(x, 9) for x in a_ang}),
                            len({round(x, 9) for x in b_ang})),
                      field("rotational_symmetry_k", B_probe.rotational_symmetry(a_union),
                            B_probe.rotational_symmetry(b_union),
                            note="B's probe applied to both unions; A has no symmetry probe"),
                      field("catalog_matches",
                            [c["catalog_id"] for c in frame["catalog_comparisons"] if c["matched"]],
                            None, True, False)]
            fields += vector_fields("sorted_polar_angle_deg", a_ang, b_ang)
            recs.append(record(f"phi_scale_vs_phi_scale_rotate_{name}_n{n}_frame{frame['frame_index']}",
                               {"shape": name, "dimension": n, "frame": frame["frame_index"],
                                "A_source": "explorative_math.projections.iter_projection_runs",
                                "label_in_both_builds": "phi composition",
                                "A_operation": "phi_scale_rotate",
                                "A_phi": "copy scaled by phi AND rotated in-plane by the golden angle",
                                "B_operation": "phi_scale",
                                "B_phi": "copy scaled by phi, not rotated (phi_stack ks 0, 1)"},
                               fields,
                               extra={"bounds": {
                                   "rotational_symmetry_k": {
                                       "a": B_probe.rotational_symmetry_bound(a_union),
                                       "b": B_probe.rotational_symmetry_bound(b_union)},
                                   "distinct_directions": {"rounding_decimals": 9}}}))
    return recs


def cmp_capabilities():
    a_only = ["superposition (k>n random packing)", "random walk trajectory",
              "iterated affine map", "elementary CA / substitution / affine enumeration",
              "nature-catalog comparison with match and miss", "append-only observer log",
              "divergence events", "manifest + validation replay", "zipapp phone build"]
    b_only = ["cut-and-project (Z^2, slope, window)", "Fibonacci chain + period profile",
              "golden-angle phyllotaxis packing", "rotational / mirror symmetry probe",
              "LZ76 + zlib compressibility", "random experiment generator",
              "open_questions.md puzzle layer", "cell24, icosahedron, lattices"]
    fields = ([field(f"capability:{c}", True, None, True, False) for c in a_only] +
              [field(f"capability:{c}", None, True, False, True) for c in b_only])
    return [record("capabilities", {"what": "features present in one build only, listed by hand"},
                   fields, {"provenance": "hand-listed from each build's code and README"})]


def cmp_schema():
    fields = []
    for k in SCHEMA_UNION:
        fields.append(field(f"schema_key:{k}", k if k in SCHEMA_A else None,
                            k if k in SCHEMA_B else None, k in SCHEMA_A, k in SCHEMA_B))
    fields += [
        field("ai_interest_rule",
              "always null in records.jsonl; filled only by an explicit observer event",
              "null in automatic runs; an AI reading or 'not interpretable atm' in hand-made runs"),
        field("interest_value_type", "boolean or null", "free text or null"),
        field("compressibility_rule",
              "null unless an observer reports a number with a named method",
              "computed automatically (LZ76 phrases + zlib ratio) with a definition"),
        field("storage", "dataset dir: manifest + records/observers/divergence/notes jsonl",
              "one JSON file per run in runs/"),
        field("status_field", "per run and per output (OBSERVED/DERIVED/PROPOSED)",
              "no status key; tags live in ai_interest text and README"),
    ]
    return [record("schema", {"what": "run schema of each build"}, fields)]


COMPARISONS = (cmp_constants, cmp_shell_formula, cmp_seed_alignment, cmp_shared_stream_sampler,
               cmp_orthogonality_stats, cmp_polytopes, cmp_rotate_project, cmp_angle_measure,
               cmp_phi_scale_vs_phi_scale_rotate, cmp_capabilities, cmp_schema)


def run_all():
    return {fn.__name__[4:]: fn() for fn in COMPARISONS}


def write_all(out_dir=OUT):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, recs in run_all().items():
        p = out_dir / f"{name}.jsonl"
        with p.open("w", encoding="utf-8") as fh:
            for r in recs:
                fh.write(json.dumps(r, sort_keys=False, allow_nan=False) + "\n")
        paths.append(p)
    return paths


def summary(results=None):
    results = results or run_all()
    lines = []
    for name, recs in results.items():
        t = {v: 0 for v in VERDICTS}
        mx = None
        for r in recs:
            for v in VERDICTS:
                t[v] += r["outputs"]["tally"][v]
            m = r["outputs"]["tally"]["max_abs_diff"]
            if m is not None:
                mx = m if mx is None else max(mx, m)
        lines.append(f"{name:24s} " + "  ".join(f"{v}={t[v]}" for v in VERDICTS)
                     + (f"  max_abs_diff={mx:.3e}" if mx is not None else ""))
    return "\n".join(lines)

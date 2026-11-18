"""adversarial_probes -- attack the probes and the run writer with inputs chosen to break them.
CC0. Co-authored. Stdlib only.
Each case has a SPEC (a docstring, a repo rule, or a ground truth computed here by
an independent method) and a PREDICTION written before the run (docs/prereg/
prereg_adversarial_q5slopes_2026-11-18.json, sha256 e31a31a1...2bdd). Observed class:
  HOLDS      the probe agrees with its spec
  DEFECT     the probe contradicts its docstring, a repo rule, or an independent truth, on in-scope input
  LIMIT      misbehaves only outside its stated scope, or in a regime its bound line names
  UNDEFINED  the spec is silent; behaviour recorded
Nothing is dropped; every case and its observed value is logged. No probe is changed.
  python3 experiments/adversarial_probes.py [seed] [--small]
--small: fewer trials in A09, A10, A14, A17 and shorter timing in A15 (own run file).
"""
import json
import math
import random
import sys
import tempfile
import time
from decimal import Decimal, getcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import emlib  # noqa: E402
from toolkit import loghook  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "adversarial_probes"
PREREG = ROOT / "docs" / "prereg" / "prereg_adversarial_q5slopes_2026-11-18.json"
NAN, INF = float("nan"), float("inf")


def finite(x):
    """Replace NaN / inf by the strings "NaN" / "Infinity" / "-Infinity" so the run file is strict JSON
    (the writer itself would emit bare tokens; that is case A19). Nothing else is changed."""
    if isinstance(x, float) and not math.isfinite(x):
        return "NaN" if x != x else ("Infinity" if x > 0 else "-Infinity")
    if isinstance(x, dict):
        return {k: finite(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [finite(v) for v in x]
    return x


def predicted_class(case_id):
    text = json.loads(PREREG.read_text())["cases"][case_id]["predict"]
    return text.split(":")[0].split(" ")[0].strip()


def rec(case_id, observed, cls, detail=None):
    pred = predicted_class(case_id)
    return {"id": case_id, "observed": finite(observed), "class": cls, "predicted": pred,
            "prediction_conforms": (cls == pred) or (pred.startswith("UNDEFINED") and cls == "UNDEFINED")
                                   or (pred == "UNDEFINED/LIMIT" and cls in ("UNDEFINED", "LIMIT")),
            "detail": detail}


def hexagon(t=0.0):
    return [[t + math.cos(math.pi * i / 3), t + math.sin(math.pi * i / 3)] for i in range(6)]


def call(fn, *a, **k):
    try:
        return {"returns": fn(*a, **k)}
    except Exception as e:  # recorded, never hidden
        return {"raises": type(e).__name__, "message": str(e)[:120]}


# ---------------------------------------------------------------- three_gap_sum
def a01():
    r = Q.three_gap_sum((1e-9, 2e-9, 7e-9))
    return rec("A01", r, "DEFECT" if r["verdict"] == "pass" else "HOLDS",
               "7e-9 is not 3e-9; verdict compares |residual| to an absolute 1e-6")


def a02():
    getcontext().prec = 60
    phi = (1 + Decimal(5).sqrt()) / 2
    scale = Decimal(10) ** 12
    a, b, c = float(scale), float(scale * phi), float(scale * phi * phi)
    r = Q.three_gap_sum((a, b, c))
    return rec("A02", r, "DEFECT" if r["verdict"] == "fail" else "HOLDS",
               "a true golden triple (1, phi, phi^2) x 1e12, each rounded to a float on its own")


def a03():
    obs = {"nan": Q.three_gap_sum((1.0, NAN, 3.0)),
           "duplicates_1_1_2": Q.three_gap_sum((1.0, 1.0, 2.0)),
           "zeros": Q.three_gap_sum((0.0, 0.0, 0.0))}
    return rec("A03", obs, "UNDEFINED", "docstring asks for distinct lengths; none of these is refused")


# ---------------------------------------------------------------- angle_histogram, polar_angles
def a04():
    h = Q.angle_histogram([-5.0, 185.0, 360.0, 90.0], width=10.0, lo=0.0, hi=180.0)
    total = sum(h["counts"])
    dropped_reported = any(k not in ("lo", "hi", "width", "counts") for k in h)
    return rec("A04", {"counts_sum": total, "angles_in": 4, "keys": sorted(h)},
               "DEFECT" if total < 4 and not dropped_reported else "HOLDS",
               "three angles outside [lo, hi] vanish with no count of them")


def a05():
    h = Q.angle_histogram([181.0], width=7.0, lo=0.0, hi=180.0)
    return rec("A05", {"n_bins": len(h["counts"]), "counted": sum(h["counts"])},
               "DEFECT" if sum(h["counts"]) == 1 else "HOLDS",
               "181 > hi = 180 but the last bin of width 7 reaches 182")


def a06():
    centre = Q.polar_angles([[-1.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    t = math.radians(-1e-15)
    tiny = Q.polar_angles([[math.cos(t), math.sin(t)], [-math.cos(t), -math.sin(t)]])
    bad = (360.0 in tiny) or centre.count(0.0) == 2
    return rec("A06", {"point_at_centroid": centre, "tiny_negative_angle": tiny},
               "DEFECT" if bad else "HOLDS",
               "the centroid point is given direction 0.0; an angle just below 0 is returned as 360.0")


# ---------------------------------------------------------------- rotational_symmetry
def abs_sym(points, k, tol):
    """Independent symmetry truth: exact bipartite matching of the k-th rotation within tol."""
    c = P.centroid(points)
    pts = [P.sub(p, c) for p in points]
    img = P.rotate(pts, 0, 1, 2 * math.pi / k)
    return has_matching(img, pts, tol)


def has_matching(a, b, tol):
    n = len(a)
    adj = [[j for j in range(n) if math.dist(a[i], b[j]) <= tol] for i in range(n)]
    match = [-1] * n

    def aug(i, seen):
        for j in adj[i]:
            if j in seen:
                continue
            seen.add(j)
            if match[j] < 0 or aug(match[j], seen):
                match[j] = i
                return True
        return False
    return len(a) == len(b) and all(aug(i, set()) for i in range(n))


def a07():
    cluster = [[0.3, 0.0], [-0.1, 0.25], [-0.15, -0.2]]
    ring = [[1e6 * math.cos(2 * math.pi * i / 12), 1e6 * math.sin(2 * math.pi * i / 12)] for i in range(12)]
    pts = cluster + ring
    got = Q.rotational_symmetry(pts, 12)
    truth = [k for k in range(2, 13) if abs_sym(pts, k, 1e-3)]
    return rec("A07", {"probe": got, "truth_abs_tol_1e-3": truth, "bound": Q.rotational_symmetry_bound(pts)},
               "DEFECT" if got != truth else "HOLDS",
               "tolerance is rel_tol times the largest radius, so a cluster far smaller than that is invisible")


def a08():
    base = [[math.cos(2 * math.pi * i / 12), math.sin(2 * math.pi * i / 12)] for i in range(12)]
    full = Q.rotational_symmetry(base, 12)

    def intact(eps):
        pts = [list(p) for p in base]
        pts[0] = [(1 + eps) * pts[0][0], (1 + eps) * pts[0][1]]
        return Q.rotational_symmetry(pts, 12) == full
    lo, hi = 1e-9, 1e-3
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if intact(mid):
            lo = mid
        else:
            hi = mid
    return rec("A08", {"full_list": full, "last_eps_still_symmetric": lo, "first_eps_broken": hi},
               "HOLDS" if 0.8e-6 <= hi <= 1.25e-6 else "DEFECT",
               "flip point of the tolerance, vertex displaced radially, scale = 1")


def a09(trials):
    rng = random.Random(9)
    disagree, bad_pos, n_sym = 0, 0, 0
    example = None
    for _ in range(trials):
        pts = []
        for _j in range(3):
            c = [rng.uniform(-1, 1), rng.uniform(-1, 1)]
            for sgn in (1, -1):
                for _m in range(rng.randint(1, 3)):
                    r, t = 1.5e-6 * math.sqrt(rng.random()), rng.uniform(0, 2 * math.pi)
                    pts.append([sgn * c[0] + r * math.cos(t), sgn * c[1] + r * math.sin(t)])
        rng.shuffle(pts)
        got = 2 in Q.rotational_symmetry(pts, 2)
        cen = P.centroid(pts)
        rel = [P.sub(p, cen) for p in pts]
        tol = 1e-6 * Q._scale_of(pts)
        img = P.rotate(rel, 0, 1, math.pi)
        truth = has_matching(img, rel, tol)
        n_sym += truth
        if got != truth:
            disagree += 1
            bad_pos += got and not truth
            example = example or pts
    obs = {"trials": trials, "exact_matching_says_symmetric": n_sym, "greedy_disagrees": disagree,
           "greedy_says_symmetric_but_none_exists": bad_pos, "example_points": example}
    return rec("A09", obs, "LIMIT" if disagree >= 1 and bad_pos == 0 else ("DEFECT" if bad_pos else "HOLDS"),
               "greedy first-fit assignment in _matches vs exact bipartite matching, near-coincident clusters")


# ---------------------------------------------------------------- mirror grid
def a10():
    rng = random.Random(10)
    n_grid, n_pairs, n_any = 0, 0, 0
    n = 500
    for _ in range(n):
        th = rng.uniform(0, 72.0)
        pent = P.regular_polygon(5, phase=math.radians(th))
        g = Q.mirror_symmetry(pent, 360)
        n_grid += len(g) == 5
        n_any += len(g) > 0
        n_pairs += len(Q.mirror_symmetry_from_pairs(pent)) == 5
    control = [len(Q.mirror_symmetry(P.regular_polygon(5, phase=math.radians(0.5 * j)), 360)) == 5
               for j in range(0, 144, 7)]
    obs = {"trials": n, "grid_finds_5_axes": n_grid, "grid_finds_any_axis": n_any, "pairs_finds_5_axes": n_pairs,
           "on_grid_controls": len(control), "on_grid_controls_found_5": sum(control),
           "bound": Q.mirror_symmetry_bound(P.regular_polygon(5))}
    ok = n_grid <= 2 and n_pairs == n and all(control)
    return rec("A10", obs, "LIMIT" if ok else "DEFECT",
               "a shape's mirror lines are found on the fixed axis grid only; the bound line states the step")


def a11():
    out, first_bad = {}, None
    for t in (1e3, 1e6, 1e8, 1e9, 1e10, 1e12):
        got = Q.rotational_symmetry(hexagon(t), 12)
        out[f"{t:g}"] = got
        if got != [2, 3, 6] and first_bad is None:
            first_bad = t
    return rec("A11", {"by_translation": out, "first_wrong_translation": first_bad,
                       "bound_at_1e10": Q.rotational_symmetry_bound(hexagon(1e10))},
               "LIMIT" if first_bad is not None else "HOLDS",
               "float spacing at the translation vs tolerance 1e-6; the bound line logs abs_tol, not float precision")


def a12():
    cube = [[x, y, z] for x in (-1.0, 1.0) for y in (-1.0, 1.0) for z in (-1.0, 1.0)]
    a, b = 0.7, 0.4
    ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
    tilt = []
    for x, y, z in cube:
        y1, z1 = ca * y - sa * z, sa * y + ca * z
        x2, z2 = cb * x + sb * z1, -sb * x + cb * z1
        tilt.append([x2, y1, z2])
    cases = {"empty": [], "one_point": [[1.0, 1.0]], "identical_points": [[1.0, 1.0]] * 3,
             "nan_coordinate": [[NAN, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 1.0]],
             "inf_coordinate": [[INF, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 1.0]],
             "cube_3d_axes_aligned": cube, "cube_3d_tilted": tilt}
    obs = {}
    for name, pts in cases.items():
        obs[name] = {"rotational_symmetry": call(Q.rotational_symmetry, pts),
                     "mirror_symmetry": call(lambda p=pts: len(Q.mirror_symmetry(p, 36))),
                     "degeneracy": call(lambda p=pts: Q.degeneracy(p)["tags"]["1e-06"]["degenerate"])}
    return rec("A12", obs, "UNDEFINED", "docstrings say 2D point set; out-of-scope and invalid input are not refused")


def a13():
    pts = hexagon() + [hexagon()[0]]
    got = Q.rotational_symmetry(pts, 12)
    return rec("A13", {"with_one_vertex_listed_twice": got, "as_a_set": Q.rotational_symmetry(hexagon(), 12)},
               "HOLDS" if got == [2, 3, 6] else "UNDEFINED",
               "the centroid is taken over the list, so a repeated point moves it; no note or refusal")


# ---------------------------------------------------------------- lz76, compressibility, periods
def ks_complexity(s):
    """Kaspar and Schuster (1987) LZ76 complexity, written independently of emlib.lz76."""
    n = len(s)
    if n == 0:
        return 0
    if n == 1:
        return 1
    c, l, i, k, kmax = 1, 1, 0, 1, 1
    while True:
        if s[i + k - 1] == s[l + k - 1]:
            k += 1
            if l + k > n:
                c += 1
                break
        else:
            kmax = max(kmax, k)
            i += 1
            if i == l:
                c += 1
                l += kmax
                if l + 1 > n:
                    break
                i, k, kmax = 0, 1, 1
            else:
                k = 1
    return c


def a14(trials):
    rng = random.Random(14)
    known = {"0001101001000101": 6, "": 0, "a": 1, "aaaaaaaa": 2, "ab": 2}
    known_obs = {s: {"expected": v, "emlib": emlib.lz76(s), "ks": ks_complexity(s)} for s, v in known.items()}
    mism = []
    for _ in range(trials):
        n = rng.randint(1, 120)
        alpha = "abcd"[:rng.choice((2, 3, 4))]
        s = "".join(rng.choice(alpha) for _ in range(n))
        if emlib.lz76(s) != ks_complexity(s):
            mism.append(s)
    structured = ["ab" * 50, "abc" * 40, "0" * 100 + "1", "".join("01"[(i * i) % 2] for i in range(200))]
    mism += [s for s in structured if emlib.lz76(s) != ks_complexity(s)]
    known_ok = all(v["expected"] == v["emlib"] for v in known_obs.values())
    return rec("A14", {"trials": trials, "mismatches": len(mism), "first_mismatch": mism[:1], "known": known_obs,
                       "ks_matches_known": all(v["expected"] == v["ks"] for v in known_obs.values())},
               "HOLDS" if not mism and known_ok else "DEFECT", "emlib.lz76 vs independent Kaspar-Schuster")


def a15(sizes):
    rng = random.Random(15)
    times = {}
    for n in sizes:
        s = "".join(rng.choice("01") for _ in range(n))
        t0 = time.perf_counter()
        emlib.lz76(s)
        times[str(n)] = time.perf_counter() - t0
    ks = list(times)
    exp = (math.log(times[ks[-1]] / times[ks[0]]) / math.log(int(ks[-1]) / int(ks[0]))
           if len(ks) > 1 and times[ks[0]] > 0 else None)
    slow = times[ks[-1]] > 2.0
    return rec("A15", {"seconds_by_n": times, "growth_exponent": exp},
               "LIMIT" if slow else "HOLDS", "repo rule: scripts run in seconds; this machine only, not a phone")


def a16():
    x = emlib.compressibility(["1", "10"], "x")
    y = emlib.compressibility(["11", "0"], "x")
    return rec("A16", {"a": x, "b": y, "same": x == y}, "UNDEFINED" if x == y else "HOLDS",
               "two different token sequences give one string, so one reading")


def border_period(s):
    n = len(s)
    if n == 0:
        return 0
    pi = [0] * n
    for i in range(1, n):
        j = pi[i - 1]
        while j and s[i] != s[j]:
            j = pi[j - 1]
        pi[i] = j + (s[i] == s[j])
    return n - pi[-1]


def a17(trials):
    rng = random.Random(17)
    known = {"abcabcab": 3, "aaaa": 1, "abcd": 4, "a": 1, "": 0}
    kobs = {s: {"expected": v, "probe": Q.smallest_period(s), "border": border_period(s)} for s, v in known.items()}
    mism = 0
    for _ in range(trials):
        n = rng.randint(0, 80)
        s = "".join(rng.choice("ab" if rng.random() < 0.7 else "abc") for _ in range(n))
        mism += Q.smallest_period(s) != border_period(s)
    ok = mism == 0 and all(v["expected"] == v["probe"] for v in kobs.values())
    return rec("A17", {"trials": trials, "mismatches": mism, "known": kobs}, "HOLDS" if ok else "DEFECT",
               "smallest_period vs n minus longest proper border")


def a18():
    obs = {str(f): Q.min_stretch_tag(100, f) for f in (-0.5, 5.0, NAN)}
    ok = all(v["assessed"] for v in obs.values())
    return rec("A18", obs, "LIMIT" if ok else "HOLDS", "frequency outside (0, 1] and NaN are accepted and assessed")


def a19():
    with tempfile.TemporaryDirectory() as d:
        p = emlib.write_run("adv_nan", 0, {"k": 1}, {"x": NAN, "y": INF, "z": -INF}, None, None, None, None, out_dir=d)
        text = Path(p).read_text()

        def refuse(name):
            raise ValueError(name)
        try:
            json.loads(text, parse_constant=refuse)
            strict = True
        except ValueError as e:
            strict = False
            tok = str(e)
    return rec("A19", {"file_has_bare_tokens": any(t in text for t in ("NaN", "Infinity")), "strict_json_parses": strict},
               "DEFECT" if not strict else "HOLDS",
               "json.dumps default allow_nan=True writes NaN / Infinity, which are not JSON")


def a20():
    out = {}
    base = dict(params={}, outputs={}, angles=None, ai_interest=None, compressibility=None, nature_match=None)
    with tempfile.TemporaryDirectory() as d:
        r1 = call(emlib.check_run, {"params": {}, "outputs": {}, "angles": None, "human_interest": "yes",
                                    "ai_interest": None, "compressibility": None, "nature_match": None,
                                    "divergence_note": None})
        r2 = call(emlib.check_run, {"params": {}})
        r3 = call(emlib.check_run, {"params": {}, "outputs": {}, "angles": None, "human_interest": None,
                                    "ai_interest": None, "compressibility": None, "nature_match": None,
                                    "divergence_note": None, "extra": 1})
        out.update(human_interest_nonnull=r1, missing_keys=r2, extra_key=r3)
        pa = emlib.write_run("adv_det", 0, {"k": 1}, {"v": [1, 2, 3]}, None, None, None, None, out_dir=d + "/a", gz=True)
        pb = emlib.write_run("adv_det", 0, {"k": 1}, {"v": [1, 2, 3]}, None, None, None, None, out_dir=d + "/b", gz=True)
        same = Path(pa).read_bytes() == Path(pb).read_bytes()
    refused = all("raises" in out[k] for k in out)
    return rec("A20", dict(out, same_bytes_twice=same), "HOLDS" if refused and same else "DEFECT",
               "check_run refuses bad records; same record gives identical gz bytes")


def run(seed=0, small=False):
    t_a09, t_a14, t_a17 = (2000, 300, 300) if small else (20000, 3000, 3000)
    sizes = (500, 1000, 2000) if small else (2000, 5000, 10000, 20000)
    cases = [a01(), a02(), a03(), a04(), a05(), a06(), a07(), a08(), a09(t_a09), a10(), a11(), a12(), a13(),
             a14(t_a14), a15(sizes), a16(), a17(t_a17), a18(), a19(), a20()]
    by_class, conf = {}, 0
    for c in cases:
        by_class[c["class"]] = by_class.get(c["class"], 0) + 1
        conf += c["prediction_conforms"]
    return {"cases": cases, "class_counts": by_class, "predictions_conforming": conf, "n_cases": len(cases),
            "defects": [c["id"] for c in cases if c["class"] == "DEFECT"],
            "limits": [c["id"] for c in cases if c["class"] == "LIMIT"],
            "undefined": [c["id"] for c in cases if c["class"] == "UNDEFINED"],
            "nonconforming_predictions": [c["id"] for c in cases if not c["prediction_conforms"]],
            "small": small}


def main(argv):
    small = "--small" in argv
    args = [a for a in argv if not a.startswith("--")]
    seed = int(args[0]) if args else 0
    out = run(seed, small)
    name = NAME + ("_small" if small else "")
    p = loghook.log_run(name, seed, {"prereg": str(PREREG.relative_to(ROOT)), "small": small,
                                     "repo_commit_under_test": "3162c9950136a4a7f14e1e064f8d335db33968c8"}, out)
    print(p)
    print({k: out[k] for k in ("class_counts", "predictions_conforming", "n_cases", "defects", "limits",
                               "undefined", "nonconforming_predictions")})


if __name__ == "__main__":
    main(sys.argv[1:])

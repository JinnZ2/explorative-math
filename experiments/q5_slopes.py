"""q5_slopes -- Q5 (how does the number of gap lengths depend on the window width?) across slopes.
CC0. Co-authored. Stdlib only.
Predictions P1..P5: docs/prereg/prereg_adversarial_q5slopes_2026-11-18.json (sha256 e31a31a1...2bdd),
written before this ran. Mode: CONFIRM for P1..P5; the record-classification table is EXPLORE (PROPOSED).

Set-up (same pipeline as window_extent.py): Z^2 cut along a line of slope s, hard window |q| <= h centred
on the line (q = (n - s m)/sqrt(1+s^2)), points kept via toolkit.composers.cut_and_project, gaps clustered
at SEP = 1e-6. Width is given in the internal coordinate n - s m as w; the canonical strip is w = 1 + |s|
(r = w / (1 + |s|) = 1). Stretch: about TARGET points; a count that is not 2 at the canonical width, or that
changes with a longer stretch, is re-read at LONG points and both readings are logged.
Rational slopes: the window edge can land exactly on lattice points, so the count is logged at
w*(1-1e-9), w, w*(1+1e-9).
  python3 experiments/q5_slopes.py [seed] [--small]
"""
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from toolkit import loghook  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402
import window_extent as WE  # noqa: E402

NAME = "q5_slopes"
SEP = 1e-6
TARGET = 1200
LONG = 20000
PHI = P.PHI


def lattice(s, h, x_max):
    nrm = math.hypot(1.0, s)
    M = int((x_max + abs(s) * h) / nrm) + 3
    out = []
    for m in range(-M, M + 1):
        c = s * m
        for n in range(math.floor(c - h * nrm) - 1, math.ceil(c + h * nrm) + 2):
            out.append((float(m), float(n)))
    return out


def gaps_at(s, w, target=TARGET):
    """distinct gap lengths for window width w (internal coordinate), about `target` points."""
    nrm = math.hypot(1.0, s)
    h = w / (2.0 * nrm)
    x_max = target / (2.0 * 2.0 * h)
    par, perp = C.line_basis_2d(s)
    pts = C.cut_and_project(lattice(s, h, x_max), par, perp, lambda q: abs(q[0]) <= h)
    xs = sorted(p[0] for p in pts if abs(p[0]) <= x_max)
    cl = WE.clusters([b - a for a, b in zip(xs, xs[1:])], SEP)
    total = sum(c["count"] for c in cl) or 1
    counts = sorted(c["count"] for c in cl)
    lengths = [c["length"] for c in cl]
    sep = min((b - a for a, b in zip(sorted(lengths), sorted(lengths)[1:])), default=None)
    return {"n_points": len(xs), "n_gap_lengths": len(cl), "lengths": lengths,
            "min_gap_frequency": counts[0] / total if counts else None,
            "min_separation_between_lengths": sep,
            "three_gap_sum": Q.three_gap_sum(lengths) if len(cl) == 3 else None}


def reading(s, w, target=TARGET, long=LONG):
    r = gaps_at(s, w, target)
    out = {"n_gap_lengths": r["n_gap_lengths"], "min_gap_frequency": r["min_gap_frequency"],
           "n_points": r["n_points"], "three_gap_sum": r["three_gap_sum"],
           "min_separation_between_lengths": r["min_separation_between_lengths"]}
    out["stretch_tag"] = Q.min_stretch_tag(r["n_points"], r["min_gap_frequency"])
    if r["n_gap_lengths"] != 2 or out["stretch_tag"]["tag"]:
        rl = gaps_at(s, w, long)
        out["long"] = {"n_points": rl["n_points"], "n_gap_lengths": rl["n_gap_lengths"],
                       "min_gap_frequency": rl["min_gap_frequency"], "three_gap_sum": rl["three_gap_sum"]}
    return out


def cf_errors(s, q_max=6):
    """|q*|s| - p| for the continued-fraction convergents p/q of |s| with q <= q_max,
    starting with the trivial terms (q=0: 1) and (q=1)."""
    x = abs(s)
    errs = [(0, 1, 1.0)]
    p0, q0, p1, q1 = 1, 0, math.floor(x), 1
    errs.append((q1, p1, abs(q1 * x - p1)))
    frac = x - math.floor(x)
    for _ in range(12):
        if frac < 1e-9:
            break
        x = 1.0 / frac
        a = math.floor(x)
        frac = x - a
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
        if q1 > q_max:
            break
        errs.append((q1, p1, abs(q1 * abs(s) - p1)))
    return errs


def records(s, m_max=6):
    """One-sided record approximations of |s|: m >= 0 where frac(m|s|) (from above) or 1 - frac (from below)
    sets a new minimum, starting from 1 at m = 0."""
    x = abs(s)
    best_up, best_dn = 1.0, 1.0
    rec = {}
    for m in range(1, m_max + 1):
        f = m * x - math.floor(m * x)
        if f < best_up - 1e-12:
            best_up = f
            rec[(m, math.floor(m * x))] = ("record_above", f)
        g = 1.0 - f
        if g < best_dn - 1e-12:
            best_dn = g
            rec[(m, math.floor(m * x) + 1)] = ("record_below", g)
    rec[(0, 1)] = ("start", 1.0)
    return rec


def scan(s, tag, rmin=0.02, rmax=2.05):
    x = abs(s)
    rec = records(s)
    seen, rows = {}, []
    for m in range(0, 7):
        for n in range(-6, 7):
            if m == 0 and n <= 0:
                continue
            d = abs(m * x - n)
            r = d / (1.0 + x)
            if not (rmin < r < rmax):
                continue
            key = round(r, 10)
            if key in seen:
                continue
            seen[key] = (m, n)
    for key, (m, n) in sorted(seen.items()):
        w = key * (1.0 + x)
        rd = reading(s, w)
        cls = rec.get((m, n), ("neither", None))[0]
        rows.append({"r": key, "m": m, "n": n, "record_class": cls, "reading": rd})
    return rows


def summarise_scan(rows):
    two = [r for r in rows if r["reading"]["n_gap_lengths"] == 2]
    rec = [r for r in rows if r["record_class"] != "neither"]
    tp = [r for r in two if r["record_class"] != "neither"]
    fp = [r for r in two if r["record_class"] == "neither"]
    fn = [r for r in rec if r["reading"]["n_gap_lengths"] != 2]
    three = [r for r in rows if r["reading"]["n_gap_lengths"] == 3]
    return {"n_widths": len(rows), "n_two": len(two), "n_record": len(rec), "two_and_record": len(tp),
            "two_not_record": [[r["r"], r["m"], r["n"]] for r in fp],
            "record_not_two": [[r["r"], r["m"], r["n"], r["reading"]["n_gap_lengths"]] for r in fn],
            "n_three": len(three),
            "three_gap_sum_all_pass": all(r["reading"]["three_gap_sum"]["verdict"] == "pass" for r in three),
            "two_widths_r": [r["r"] for r in two]}


def main(argv):
    small = "--small" in argv
    args = [a for a in argv if not a.startswith("--")]
    seed = int(args[0]) if args else 0
    named = {"1/phi": 1 / PHI, "sqrt2-1": math.sqrt(2) - 1, "sqrt3-1": math.sqrt(3) - 1, "1/sqrt2": 1 / math.sqrt(2),
             "e-2": math.e - 2, "pi-3": math.pi - 3, "-1/phi": -1 / PHI, "phi": PHI}
    rational = {"2/5": 0.4, "1/2": 0.5, "1": 1.0, "0": 0.0, "3": 3.0, "100": 100.0, "1e-6": 1e-6}
    rng = random.Random(seed + 5)
    n_rand = 20 if small else 150
    rand = []
    while len(rand) < n_rand:
        s = rng.uniform(-3, 3)
        if abs(s) > 0.02:
            rand.append(s)

    canonical = {}
    for name, s in list(named.items()) + [(f"rand{i}", s) for i, s in enumerate(rand)]:
        canonical[name] = {"slope": s, **reading(s, 1.0 + abs(s))}
    rational_rows = {}
    for name, s in rational.items():
        w = 1.0 + abs(s)
        rational_rows[name] = {"slope": s, "rational": True,
                               "below": reading(s, w * (1 - 1e-9)), "at": reading(s, w),
                               "above": reading(s, w * (1 + 1e-9))}
    scans = {}
    for name, s in named.items():
        rows = scan(s, name)
        scans[name] = {"slope": s, "cf_errors_q_le_6": [list(e) for e in cf_errors(s)],
                       "summary": summarise_scan(rows), "rows": rows}
    n_two_rand = sum(1 for k, v in canonical.items() if k.startswith("rand") and v["n_gap_lengths"] == 2)
    not_two_irr = {k: v["n_gap_lengths"] for k, v in canonical.items() if v["n_gap_lengths"] != 2}
    long_fix = {k: v["long"]["n_gap_lengths"] for k, v in canonical.items() if "long" in v}
    # P3: convergent errors are two-length widths
    p3 = {}
    for name, sc in scans.items():
        s = sc["slope"]
        by_r = {r["r"]: r["reading"]["n_gap_lengths"] for r in sc["rows"]}
        miss = []
        for q, p, d in cf_errors(s):
            rr = round(d / (1.0 + abs(s)), 10)
            if 0.02 < rr < 2.05 and by_r.get(rr) != 2:
                miss.append([q, p, rr, by_r.get(rr)])
        p3[name] = {"convergent_errors_not_two": miss}
    # P4: golden widths
    def is_phi_power(r):
        k = round(math.log(r) / math.log(PHI))
        return abs(r - PHI ** k) < 1e-8, k
    p4 = {name: {"two_widths_r": sc["summary"]["two_widths_r"],
                 "phi_powers": [[r, is_phi_power(r)[1]] for r in sc["summary"]["two_widths_r"] if is_phi_power(r)[0]],
                 "not_phi_powers": [r for r in sc["summary"]["two_widths_r"] if not is_phi_power(r)[0]]}
          for name, sc in scans.items()}
    outputs = {
        "canonical": canonical, "rational": rational_rows, "scans": scans, "P3_convergent_errors": p3, "P4_phi_powers": p4,
        "P1_two_at_canonical": {"named_irrational": {k: canonical[k]["n_gap_lengths"] for k in named},
                                "random_two": n_two_rand, "random_total": len(rand),
                                "not_two_at_canonical": not_two_irr, "not_two_after_long_stretch": long_fix},
        "P2_rational_counts": {k: [v["below"]["n_gap_lengths"], v["at"]["n_gap_lengths"], v["above"]["n_gap_lengths"]]
                               for k, v in rational_rows.items()},
        "P5_three_gap_sum_all_pass": all(sc["summary"]["three_gap_sum_all_pass"] for sc in scans.values()),
        "mode": "CONFIRM P1-P5; record classification EXPLORE, status PROPOSED", "small": small}
    name = NAME + ("_small" if small else "")
    p = loghook.log_run(name, seed, {"prereg": "docs/prereg/prereg_adversarial_q5slopes_2026-11-18.json",
                                     "target_points": TARGET, "long_points": LONG, "sep": SEP,
                                     "n_random_slopes": len(rand), "small": small,
                                     "repo_commit_under_test": "3162c9950136a4a7f14e1e064f8d335db33968c8"}, outputs)
    print(p)
    print("P1 named:", outputs["P1_two_at_canonical"]["named_irrational"])
    print("P1 random two:", n_two_rand, "/", len(rand), "| not two:", not_two_irr, "| after long:", long_fix)
    print("P2 rational [below,at,above]:", outputs["P2_rational_counts"])
    for k, sc in scans.items():
        s = sc["summary"]
        print("scan", k, "widths", s["n_widths"], "two", s["n_two"], "two&record", s["two_and_record"], "two_not_record", s["two_not_record"],
              "record_not_two", s["record_not_two"], "three_ok", s["three_gap_sum_all_pass"])
    print("P3 misses:", {k: v["convergent_errors_not_two"] for k, v in p3.items()})
    print("P4:", {k: (len(v["phi_powers"]), len(v["not_phi_powers"])) for k, v in p4.items()})
    print("P5:", outputs["P5_three_gap_sum_all_pass"])


if __name__ == "__main__":
    main(sys.argv[1:])

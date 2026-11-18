"""q5_slopes_followup -- EXPLORE follow-up to q5_slopes (not pre-registered; every reading is PROPOSED).
CC0. Co-authored. Stdlib only.
Expectations written before the run (they can be wrong; the run decides):
  E1  record approximations with m up to 20 (one-sided records of |s|, see q5_slopes.records) still give
      exactly 2 gap lengths, for the 5 irrational slopes outside Q(sqrt5), on a centred window.
  E2  moving the window off-centre (centre c = u * w / 2 in the internal coordinate, u in (-1, 1)) does not
      change the count (2 stays 2, 3 stays 3) at: the canonical width, five record widths, five non-record widths.
      Where it does change, the cases are logged.
  E3  the golden-slope convergent |5*phi - 8| (n = 8, outside the q5_slopes scan) gives 2 lengths.
Stretch is TARGET points; counts that are not 2 or 3 are logged as read. --small uses 2 slopes and m <= 10.
  python3 experiments/q5_slopes_followup.py [seed] [--small]
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
import window_extent as WE  # noqa: E402
import q5_slopes as S  # noqa: E402

NAME = "q5_slopes_followup"


def count_at(s, w, c_frac=0.0, target=S.TARGET):
    """gap-length count for window q in [c - h, c + h], c = c_frac * h, w in n - s m units."""
    nrm = math.hypot(1.0, s)
    h = w / (2.0 * nrm)
    c = c_frac * h
    x_max = target / (2.0 * 2.0 * h)
    par, perp = C.line_basis_2d(s)
    pts = C.cut_and_project(S.lattice(s, h + abs(c), x_max), par, perp, lambda q: abs(q[0] - c) <= h)
    xs = sorted(p[0] for p in pts if abs(p[0]) <= x_max)
    cl = WE.clusters([b - a for a, b in zip(xs, xs[1:])], S.SEP)
    counts = sorted(k["count"] for k in cl)
    return {"n_gap_lengths": len(cl), "n_points": len(xs),
            "min_gap_frequency": (counts[0] / sum(counts)) if counts else None}


def main(argv):
    small = "--small" in argv
    args = [a for a in argv if not a.startswith("--")]
    seed = int(args[0]) if args else 0
    rng = random.Random(seed + 77)
    slopes = {"sqrt2-1": math.sqrt(2) - 1, "sqrt3-1": math.sqrt(3) - 1, "1/sqrt2": 1 / math.sqrt(2),
              "e-2": math.e - 2, "pi-3": math.pi - 3}
    if small:
        slopes = dict(list(slopes.items())[:2])
    m_max = 10 if small else 20
    e1, e2 = {}, {}
    for name, s in slopes.items():
        x = abs(s)
        recs = S.records(s, m_max)
        rows = []
        for (m, n), (cls, d) in sorted(recs.items()):
            w = abs(m * x - n)
            if w < 0.02:
                rows.append({"m": m, "n": n, "w": w, "class": cls, "skipped": "w < 0.02 (stretch cost)"})
                continue
            rows.append({"m": m, "n": n, "w": w, "class": cls, **count_at(s, w)})
        e1[name] = rows
        # E2: offsets
        rec_w = [abs(m * x - n) for (m, n) in recs if abs(m * x - n) >= 0.02][:5]
        non_rec = []
        for m in range(1, 7):
            for n in range(-6, 7):
                w = abs(m * x - n)
                if (m, n) not in recs and 0.1 < w < 2.0 and len(non_rec) < 5 and (m + n) % 2 == 0:
                    non_rec.append(w)
        cases = [("canonical", 1.0 + x)] + [("record", w) for w in rec_w] + [("non_record", w) for w in non_rec]
        out = []
        for kind, w in cases:
            base = count_at(s, w, 0.0)["n_gap_lengths"]
            offs = {f"{u:+.2f}": count_at(s, w, u)["n_gap_lengths"] for u in (-0.9, -0.4, 0.3, 0.8)}
            out.append({"kind": kind, "w": w, "centred": base, "offsets": offs,
                        "changed": any(v != base for v in offs.values())})
        e2[name] = out
    phi = P.PHI
    w3 = abs(5 * (1 / phi) - 3)  # |5 s - 3| with s = 1/phi: the n = 3 form of the same error
    w8 = abs(5 * phi - 8)        # slope phi, n = 8
    e3 = {"slope_phi_w_|5phi-8|": {"w": w8, **count_at(phi, w8)},
          "slope_1/phi_w_|5/phi-3|": {"w": w3, **count_at(1 / phi, w3)}}
    summary = {"E1_all_two": all(r.get("n_gap_lengths") == 2 for rows in e1.values() for r in rows if "skipped" not in r),
               "E1_rows_not_two": [[k, r["m"], r["n"], r["n_gap_lengths"]] for k, rows in e1.items() for r in rows
                                   if "skipped" not in r and r["n_gap_lengths"] != 2],
               "E2_changed": [[k, o["kind"], o["w"], o["centred"], o["offsets"]] for k, v in e2.items() for o in v if o["changed"]],
               "E3": {k: v["n_gap_lengths"] for k, v in e3.items()}}
    outputs = {"E1_records": e1, "E2_offsets": e2, "E3_convergent_n8": e3, "summary": summary,
               "mode": "EXPLORE", "status": "PROPOSED", "small": small}
    name = NAME + ("_small" if small else "")
    p = loghook.log_run(name, seed, {"mode": "EXPLORE", "target_points": S.TARGET, "m_max": m_max,
                                     "repo_commit_under_test": "3162c9950136a4a7f14e1e064f8d335db33968c8"}, outputs)
    print(p)
    print(summary)


if __name__ == "__main__":
    main(sys.argv[1:])

"""am. simplicial_helix -- tan(n theta) = n tan(theta) on every logged closed chain.

CC0. Co-authored. Stdlib only. Status: PROPOSED.

Equation carried from REFERENCES/related_studies.md R5 (keywords only; the
paper was not read here). For each of the 17 closed chains (corkscrew:
tetrahedra in R^4; corkscrew6: m-simplices, m >= 4, in R^6):

  n      the simplex's vertex count, m + 1
  theta  HALF of each plane's turn: theta = 180 * p / N degrees. (The full
         turn is logged beside it. The half turn is the one for which the
         R^3 tetrahedral helix, turn acos(-2/3), satisfies the equation with
         n = 4; that row is logged as the reference.)

Logged per plane and per chain: tan(n theta) - n tan(theta), and the bounded
form sin(n theta) cos(theta) - n sin(theta) cos(n theta), which is zero for
the same theta and does not blow up. The four chains that matched no
Coxeter group are flagged. Observer slots null.

  python3 experiments/simplicial_helix.py [seed]
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402

NAME = "simplicial_helix"
SOURCE = "runs/coxeter_exponents_seed0.json"
TOL = 1e-9


def residuals(theta_deg, n):
    t = math.radians(theta_deg)
    bounded = math.sin(n * t) * math.cos(t) - n * math.sin(t) * math.cos(n * t)
    c1, cn = math.cos(t), math.cos(n * t)
    raw = (math.tan(n * t) - n * math.tan(t)) if abs(c1) > 1e-12 and abs(cn) > 1e-12 else None
    return {"theta_deg": theta_deg, "n": n, "tan_n_theta": math.tan(n * t) if abs(cn) > 1e-12 else None,
            "n_tan_theta": n * math.tan(t) if abs(c1) > 1e-12 else None,
            "residual_tan_form": raw, "residual_bounded_form": bounded,
            "holds": abs(bounded) < TOL}


def run(seed=0):
    rec = emlib.load_run(emlib.ROOT / SOURCE)
    rows = []
    for c in rec["outputs"]["chains"]:
        n = c["simplex_m"] + 1
        planes = []
        for p in c["turns"]:
            half = residuals(180.0 * p / c["n"], n)
            full = residuals(360.0 * p / c["n"], n)
            planes.append({"turn": p, "turn_deg": 360.0 * p / c["n"], "half_turn": half,
                           "full_turn": full,
                           "theta_is_90_deg": 2 * p == c["n"],
                           "other_n_2_to_12_that_hold_at_half_turn":
                               [k for k in range(2, 13)
                                if residuals(180.0 * p / c["n"], k)["holds"]]})
        rows.append({"dimension": c["dimension"], "simplex_m": c["simplex_m"],
                     "simplex_vertices_n": n, "n_to_close": c["n"], "turns": c["turns"],
                     "coxeter_verdict": c["verdict"], "coxeter_matches": c["matches"],
                     "survivor_unmatched": c["verdict"] == "no match",
                     "planes": planes,
                     "largest_abs_residual_half_turn":
                         max(abs(x["half_turn"]["residual_bounded_form"]) for x in planes),
                     "largest_abs_residual_full_turn":
                         max(abs(x["full_turn"]["residual_bounded_form"]) for x in planes),
                     "holds_on_every_plane_half_turn": all(x["half_turn"]["holds"] for x in planes),
                     "planes_holding_half_turn": sum(1 for x in planes if x["half_turn"]["holds"]),
                     "planes_holding_only_at_theta_90":
                         sum(1 for x in planes if x["half_turn"]["holds"] and x["theta_is_90_deg"]),
                     "holds_on_every_plane_full_turn": all(x["full_turn"]["holds"] for x in planes)})
    twist = math.degrees(math.acos(-2.0 / 3.0))
    ref = {"helix": "regular tetrahedra in R^3, turn acos(-2/3)", "turn_deg": twist,
           "half_turn": residuals(twist / 2.0, 4), "full_turn": residuals(twist, 4),
           "half_turn_with_n_3": residuals(twist / 2.0, 3)}
    surv = [r for r in rows if r["survivor_unmatched"]]
    return {"equation": "tan(n theta) = n tan(theta)",
            "source": "REFERENCES/related_studies.md R5 (carried; keywords only)",
            "bound": {"seed": seed, "tol": TOL, "theta": "half of each plane's turn",
                      "n": "simplex vertex count, m + 1",
                      "bounded_form": "sin(n theta) cos(theta) - n sin(theta) cos(n theta)",
                      "theta_choice_not_from_the_paper": "the paper was not read; the half "
                          "turn is the reading under which the R^3 reference row holds"},
            "reference_r3": ref, "chains": rows,
            "summary": {"chains": len(rows),
                        "holding_on_every_plane": sum(1 for r in rows
                                                      if r["holds_on_every_plane_half_turn"]),
                        "holding_on_no_plane": sum(1 for r in rows
                                                   if r["planes_holding_half_turn"] == 0),
                        "planes_holding": sum(r["planes_holding_half_turn"] for r in rows),
                        "planes_holding_that_are_theta_90":
                            sum(r["planes_holding_only_at_theta_90"] for r in rows),
                        "theta_90_note": "at theta = 90 degrees both tangents are unbounded; "
                                         "the bounded form is 0 there for odd n",
                        "survivors": len(surv),
                        "survivors_holding_on_every_plane":
                            sum(1 for r in surv if r["holds_on_every_plane_half_turn"]),
                        "matched_holding_on_every_plane":
                            sum(1 for r in rows if not r["survivor_unmatched"]
                                and r["holds_on_every_plane_half_turn"]),
                        "full_turn_reading_holding": sum(1 for r in rows
                                                         if r["holds_on_every_plane_full_turn"])}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(
        ["h" if r["holds_on_every_plane_half_turn"] else "n" for r in out["chains"]],
        "one letter per chain: h the equation holds on every plane (half turn), n it does not")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"source": SOURCE, "tol": TOL,
                                     "status_before_run": "PROPOSED"},
                        out, None, None, comp, None)
    print(p)
    r = out["reference_r3"]
    print("R3 reference: half turn, n=4:", r["half_turn"]["residual_bounded_form"],
          "| full turn:", r["full_turn"]["residual_bounded_form"],
          "| half turn, n=3:", r["half_turn_with_n_3"]["residual_bounded_form"])
    for c in out["chains"]:
        print(f"  {c['dimension']}D m={c['simplex_m']} n={c['simplex_vertices_n']} N={c['n_to_close']:3d} "
              f"turns {c['turns']}: residuals "
              f"{[float(f'{x['half_turn']['residual_bounded_form']:.2e}') for x in c['planes']]} "
              f"holds {c['holds_on_every_plane_half_turn']}"
              f"{'   <-- unmatched' if c['survivor_unmatched'] else ''}")
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])

"""au. h_order -- one guess about which conjugate turn sets are chains, set
beside the verdicts conjugate_map logged.

CC0. Co-authored. Stdlib only. Status: PROPOSED.

H-ORDER (PROPOSED, a guess from structure, not a known result):
  u P is REAL  <=>  multiplying by u keeps the ORDER of the reduced turn
  angles: with P = (p_1 < ... < p_k) and q_i = min(u p_i mod N, N - u p_i mod N),
  order is kept when q_1 < q_2 < ... < q_k.

Input: runs/conjugate_map_seed0.json.gz, verdicts as logged (exact for the
logged and R^5 / R^6 chains, float for the R^7 list). Nothing is solved here.

  1  cross-tab   one row per (chain P, orbit member S). A member is REAL when
                 it is a found chain, else its logged verdict. Several units
                 reach the same S; order_kept is read three ways and all
                 three are logged:
                   first   the smallest unit that reaches S
                   any     some unit reaching S keeps the order
                   every   every unit reaching S keeps the order
                 The row P -> P (u = 1) keeps the order and is REAL by
                 construction; tables are given with and without those rows.
  2  counterexamples   order kept and a share < 0; order broken and REAL.
  3  variant     the order of cos(2 pi q_i / N) in place of the order of q_i.
                 DERIVED: q_i lies in [0, N / 2], where the cosine falls
                 strictly, so the cosines are in the reverse order exactly
                 when the angles are in order. The variant is the same test;
                 it is run and its table logged beside the first.

No further fitting. Nothing is interpreted. Observer slots null.

  python3 experiments/h_order.py [seed]
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402

NAME = "h_order"
SOURCE = "runs/conjugate_map_seed0.json.gz"


def images(n, turns, u):
    return [min((u * p) % n, n - (u * p) % n) for p in turns]


def order_kept(n, turns, u):
    q = images(n, sorted(turns), u)
    return all(a < b for a, b in zip(q, q[1:]))


def cosine_order_kept(n, turns, u):
    c = [math.cos(2 * math.pi * q / n) for q in images(n, sorted(turns), u)]
    return all(a > b for a, b in zip(c, c[1:]))


def units_reaching(n, turns, target):
    t = tuple(sorted(target))
    return [u for u in range(1, n) if math.gcd(u, n) == 1
            and tuple(sorted(images(n, turns, u))) == t]


def read(test):
    """One row per (chain, orbit member)."""
    rec = emlib.load_run(emlib.ROOT / SOURCE)["outputs"]
    rows = []
    for r in rec["orbits"]:
        n, p = r["n"], r["turns"]
        for x in r["orbit"]:
            us = units_reaching(n, p, x["turns"])
            kept = [test(n, p, u) for u in us]
            rows.append({"source": r["source"], "simplex_m": r["simplex_m"], "n": n,
                         "chain": p, "member": x["turns"], "is_the_chain_itself": x["turns"] == p,
                         "verdict_logged": x["verdict"],
                         "real": x["verdict"] in ("FOUND", "REAL"),
                         "units_reaching": us, "order_kept_by_unit": kept,
                         "first": kept[0], "any": any(kept), "every": all(kept)})
    return rows


def table(rows, how):
    t = {"kept_real": 0, "kept_negative": 0, "broken_real": 0, "broken_negative": 0}
    for r in rows:
        t[("kept" if r[how] else "broken") + ("_real" if r["real"] else "_negative")] += 1
    t["rows"] = len(rows)
    t["counterexamples"] = t["kept_negative"] + t["broken_real"]
    t["h_order_holds_on_every_row"] = t["counterexamples"] == 0
    return t


def first_rows(rows, how, kept, real, k=10):
    out = [r for r in rows if r[how] == kept and r["real"] == real]
    return [{"simplex_m": r["simplex_m"], "n": r["n"], "chain": r["chain"],
             "member": r["member"], "verdict_logged": r["verdict_logged"],
             "units_reaching": r["units_reaching"],
             "order_kept_by_unit": r["order_kept_by_unit"]} for r in out[:k]]


def cross(rows):
    others = [r for r in rows if not r["is_the_chain_itself"]]
    out = {}
    for how in ("first", "any", "every"):
        out[how] = {"all_rows": table(rows, how), "without_the_chain_itself": table(others, how),
                    "order_kept_and_a_negative_share_first_10": first_rows(rows, how, True, False),
                    "order_broken_and_real_first_10": first_rows(rows, how, False, True)}
    return out


def run(seed=0):
    rows = read(order_kept)
    rows_cos = read(cosine_order_kept)
    same = all(a["order_kept_by_unit"] == b["order_kept_by_unit"] for a, b in zip(rows, rows_cos))
    not_r7 = [r for r in rows if not r["source"].endswith("R^7")]
    r7 = [r for r in rows if r["source"].endswith("R^7")]
    c = cross(rows)
    disagree = sum(1 for r in rows if r["any"] != r["every"])
    return {"bound": {"seed": seed, "source": SOURCE,
                      "verdicts": "as logged by conjugate_map: exact for the logged and "
                                  "R^5 / R^6 chains, float for the R^7 list; none solved here",
                      "row": "one per (chain, orbit member); the chain itself is a row (u = 1)",
                      "units": "every unit u mod N, 1 <= u < N; u and N - u give the same images",
                      "variant": "cosines compared as floats; ties cannot occur for distinct "
                                 "reduced turns below N / 2"},
            "hypothesis": "H-ORDER, PROPOSED: u P is REAL <=> u keeps the order of the reduced "
                          "turn angles",
            "rows": rows,
            "cross_tab": c,
            "cross_tab_logged_and_r5_r6_only": cross(not_r7),
            "cross_tab_r7_only": cross(r7),
            "variant_cosine_order": {"cross_tab": cross(rows_cos),
                                     "same_reading_as_angle_order_on_every_unit": same,
                                     "why": "DERIVED: cos is strictly decreasing on [0, pi]"},
            "summary": {
                "rows": len(rows), "rows_from_the_r7_list": len(r7),
                "rows_real": sum(1 for r in rows if r["real"]),
                "rows_negative": sum(1 for r in rows if not r["real"]),
                "rows_where_any_and_every_differ": disagree,
                "counterexamples": {how: c[how]["all_rows"]["counterexamples"]
                                    for how in ("first", "any", "every")},
                "h_order_holds": {how: c[how]["all_rows"]["h_order_holds_on_every_row"]
                                  for how in ("first", "any", "every")},
                "variant_is_the_same_test": same,
                "named_case_28_1_8_13_to_5_9_12": next(
                    {"units_reaching": r["units_reaching"],
                     "order_kept_by_unit": r["order_kept_by_unit"], "real": r["real"]}
                    for r in rows if r["n"] == 28 and r["chain"] == [1, 8, 13]
                    and r["member"] == [5, 9, 12])}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(
        [("k" if r["any"] else "b") + ("R" if r["real"] else "n") for r in out["rows"]],
        "two letters per (chain, member) row: k order kept by some unit / b broken by all; "
        "R real / n negative share")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"sources": [SOURCE], "status_before_run": "PROPOSED"},
                        out, None, None, comp, None)
    print(p)
    for name, c in (("ALL", out["cross_tab"]), ("LOGGED+R5/R6", out["cross_tab_logged_and_r5_r6_only"]),
                    ("R7", out["cross_tab_r7_only"]),
                    ("VARIANT cosine order, ALL", out["variant_cosine_order"]["cross_tab"])):
        for how in ("first", "any", "every"):
            for part in ("all_rows", "without_the_chain_itself"):
                t = c[how][part]
                print(f"{name:26s} {how:5s} {part:24s} kept&REAL {t['kept_real']:4d}  kept&NEG "
                      f"{t['kept_negative']:4d}  broken&REAL {t['broken_real']:4d}  broken&NEG "
                      f"{t['broken_negative']:4d}  counterexamples {t['counterexamples']}")
    for how in ("any", "every"):
        for key in ("order_kept_and_a_negative_share_first_10", "order_broken_and_real_first_10"):
            print(how, key)
            for r in out["cross_tab"][how][key]:
                print("   ", r)
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])

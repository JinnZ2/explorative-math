"""al. coxeter_products -- the unmatched chains against reducible Coxeter groups.

CC0. Co-authored. Stdlib only. Status: PROPOSED.

Four closed chains logged by corkscrew6 matched no irreducible group in
coxeter_exponents. Here they are set beside products of two or three
irreducible groups of total rank <= MAX_RANK, with factors A_n (n >= 1),
B_n, D_n, G2, F4 and I2(m).

For a product, the Coxeter element is the product of the factors' Coxeter
elements: its order is the lcm of their Coxeter numbers, and its eigenvalue
angles are the fractions m_j / h_k of every factor. Each factor's h and
exponents are computed from its diagram (coxeter_exponents.exponents).

A chain with turns p / N matches a product when lcm = N and every turn
fraction, or its complement, is one of the product's fractions. Logged per
chain: every match, or none. Observer slots null.

  python3 experiments/coxeter_products.py [seed]
"""
import itertools
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import coxeter_exponents as X  # noqa: E402
import emlib  # noqa: E402

NAME = "coxeter_products"
MAX_RANK = 6
I2_MAX = 60


def factor(kind, n, m=None):
    """(name, rank, h, exponents), computed from the diagram."""
    if kind == "A" and n == 1:
        edges = {}
    elif kind == "G":
        edges = {(0, 1): 6}
    elif kind == "F":
        edges = {(0, 1): 3, (1, 2): 4, (2, 3): 3}
    elif kind == "I":
        edges = {(0, 1): m}
    else:
        edges = X.diagram(kind, n)
    c = X.coxeter_element(edges, n)
    ident = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    import cmath
    traces, p, h = [float(n)], ident, None
    for k in range(1, 400):
        p = X.matmul(p, c)
        if max(abs(p[i][j] - ident[i][j]) for i in range(n) for j in range(n)) < 1e-8:
            h = k
            break
        traces.append(sum(p[i][i] for i in range(n)))
    exps = []
    for e in range(1, h):
        mult = sum(traces[j] * cmath.exp(-2j * math.pi * j * e / h) for j in range(h)) / h
        exps += [e] * int(round(mult.real))
    name = {"A": f"A{n}", "B": f"B{n}", "D": f"D{n}", "G": "G2", "F": "F4",
            "I": f"I2({m})"}[kind]
    return {"name": name, "rank": n, "h": h, "exponents": exps}


def factors():
    out = [factor("A", n) for n in range(1, MAX_RANK)]
    out += [factor("B", n) for n in range(2, MAX_RANK)]
    out += [factor("D", n) for n in range(4, MAX_RANK)]
    out += [factor("G", 2), factor("F", 4)]
    out += [factor("I", 2, m) for m in range(5, I2_MAX + 1) if m != 6]   # I2(6) is G2
    return out


def lcm(values):
    out = 1
    for v in values:
        out = out * v // math.gcd(out, v)
    return out


def products(fs):
    for size in (2, 3):
        for combo in itertools.combinations_with_replacement(range(len(fs)), size):
            parts = [fs[i] for i in combo]
            if sum(p["rank"] for p in parts) <= MAX_RANK:
                yield parts


def unmatched():
    rec = emlib.load_run(emlib.RUNS / "coxeter_exponents_seed0.json")
    return [r for r in rec["outputs"]["chains"] if r["verdict"] == "no match"]


def run(seed=0):
    fs = factors()
    prods = list(products(fs))
    rows = []
    for c in unmatched():
        n, turns = c["n"], c["turns"]
        want = [Fraction(t, n) for t in turns]
        present, full = [], []
        for parts in prods:
            fr = {Fraction(e, p["h"]) for p in parts for e in p["exponents"]}
            if all(w in fr or (1 - w) in fr for w in want):
                name = " x ".join(p["name"] for p in parts)
                h = lcm([p["h"] for p in parts])
                present.append({"product": name, "rank": sum(p["rank"] for p in parts), "h": h})
                if h == n:
                    full.append(name)
        # which single turn fractions any factor offers at all
        per_turn = {str(w): sorted({f["name"] for f in fs
                                    if any(Fraction(e, f["h"]) in (w, 1 - w)
                                           for e in f["exponents"])})[:12] for w in want}
        rows.append({"dimension": c["dimension"], "simplex_m": c["simplex_m"], "n": n,
                     "turns": turns, "turn_fractions": [str(w) for w in want],
                     "matches": full, "verdict": "match" if full else "none",
                     "products_holding_all_turn_fractions_whatever_h":
                         [p["product"] for p in present][:20],
                     "factors_offering_each_turn_fraction": per_turn,
                     "turn_fractions_no_factor_offers":
                         [k for k, v in per_turn.items() if not v]})
    return {"bound": {"seed": seed, "max_total_rank": MAX_RANK,
                      "factors": "A1..A5, B2..B5, D4, D5, G2, F4, I2(m) for m = 5..60 (m != 6)",
                      "product_sizes": [2, 3], "n_factors": len(fs), "n_products": len(prods),
                      "match": "lcm of the factors' h equals N, and every turn fraction p/N or "
                               "its complement is m_j / h_k of some factor",
                      "not_tested": ["H3, H4, E6-E8 as factors", "total rank above 6",
                                     "four or more factors"]},
            "factors": [f for f in fs if not f["name"].startswith("I2")],
            "i2_factors": "I2(m): h = m, exponents 1 and m - 1 (computed for each m)",
            "chains": rows,
            "summary": {"chains": len(rows), "match": sum(1 for r in rows if r["matches"]),
                        "none": sum(1 for r in rows if not r["matches"])}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(["m" if r["matches"] else "n" for r in out["chains"]],
                                 "one letter per chain: m match, n none")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"max_total_rank": MAX_RANK, "i2_max": I2_MAX,
                                     "source": "runs/coxeter_exponents_seed0.json",
                                     "status_before_run": "PROPOSED"},
                        out, None, None, comp, None)
    print(p)
    print("factors:", [(f["name"], f["h"], f["exponents"]) for f in out["factors"]])
    print("products tried:", out["bound"]["n_products"])
    for r in out["chains"]:
        print(f"  N={r['n']} turns {r['turns']} = {r['turn_fractions']}: {r['verdict']} "
              f"{r['matches']} | any h: {r['products_holding_all_turn_fractions_whatever_h'][:4]}"
              f" | no factor offers {r['turn_fractions_no_factor_offers']}")
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])

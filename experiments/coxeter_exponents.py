"""ak. coxeter_exponents -- are the closed chains' turns Coxeter exponents over h?

CC0. Co-authored. Stdlib only. Status: PROPOSED.

Every closed chain logged by corkscrew (tetrahedra in R^4) and corkscrew6
(m-simplices, m >= 4, in R^6) has turns 360 * p / N per plane. Here the
exponents m_j and the Coxeter number h of A_n, B_n, D_n (rank <= 8), H3,
H4, E6, E7, E8 are COMPUTED, not carried: the Coxeter element is built from
the Coxeter diagram as a matrix, h is its order, and the exponents are the
k with eigenvalue exp(2 pi i k / h), counted from the traces of its powers.

A chain matches a group when N = h and every turn p (or N - p) is an
exponent. Logged per chain: every match, the matches whose rank equals the
chain's dimension, and the matches whose exponents are exactly the chain's
turns and their complements. Observer slots null.

  python3 experiments/coxeter_exponents.py [seed]
"""
import cmath
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402

NAME = "coxeter_exponents"
MAX_RANK = 8
SOURCES = {"corkscrew": "runs/corkscrew_seed0.json", "corkscrew6": "runs/corkscrew6_seed0.json"}


def diagram(kind, n):
    """Coxeter matrix entries m_ij > 2 as {(i, j): m}."""
    e = {(i, i + 1): 3 for i in range(n - 1)}
    if kind == "B":
        e[(n - 2, n - 1)] = 4
    elif kind == "D":
        e = {(i, i + 1): 3 for i in range(n - 2)}
        e[(n - 3, n - 1)] = 3
    elif kind == "H":
        e[(0, 1)] = 5
    elif kind == "E":
        e = {(i, i + 1): 3 for i in range(n - 2)}
        e[(2, n - 1)] = 3
    return e


def coxeter_element(edges, n):
    """Product of the simple reflections, in the simple-root basis:
    s_i(a_j) = a_j - 2 B_ij a_i with B_ii = 1, B_ij = -cos(pi / m_ij)."""
    b = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for (i, j), m in edges.items():
        b[i][j] = b[j][i] = -math.cos(math.pi / m)
    c = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for i in range(n):
        s = [[(1.0 if r == k else 0.0) - (2.0 * b[i][k] if r == i else 0.0) for k in range(n)]
             for r in range(n)]
        c = matmul(c, s)
    return c


def matmul(a, b):
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def exponents(kind, n, max_order=200):
    c = coxeter_element(diagram(kind, n), n)
    ident = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    traces, p, h = [float(n)], ident, None
    for k in range(1, max_order + 1):
        p = matmul(p, c)
        if max(abs(p[i][j] - ident[i][j]) for i in range(n) for j in range(n)) < 1e-8:
            h = k
            break
        traces.append(sum(p[i][i] for i in range(n)))
    exps = []
    for m in range(1, h):
        mult = sum(traces[j] * cmath.exp(-2j * math.pi * j * m / h) for j in range(h)) / h
        exps += [m] * int(round(mult.real))
    return h, exps


def groups():
    out = []
    for kind, ranks in (("A", range(2, MAX_RANK + 1)), ("B", range(2, MAX_RANK + 1)),
                        ("D", range(4, MAX_RANK + 1)), ("H", (3, 4)), ("E", (6, 7, 8))):
        for n in ranks:
            h, e = exponents(kind, n)
            out.append({"group": f"{kind}{n}", "rank": n, "h": h, "exponents": e,
                        "n_exponents_equals_rank": len(e) == n,
                        "exponents_sum_to_rank_h_over_2": sum(e) * 2 == n * h})
    return out


def chains():
    out = []
    r4 = emlib.load_run(emlib.ROOT / SOURCES["corkscrew"])["outputs"]["tetrahedra_chain"]["r4"]
    for c in r4["closed_chains"]:
        out.append({"source": "corkscrew", "dimension": 4, "simplex_m": 3,
                    "n": c["n_tetrahedra_to_close"], "turns": [c["p"], c["q"]]})
    r6 = emlib.load_run(emlib.ROOT / SOURCES["corkscrew6"])["outputs"]["chains"]
    for m in ("4", "5", "6"):
        for c in r6["chains_using_all_three_planes"][m]["listed"]:
            out.append({"source": "corkscrew6", "dimension": 6, "simplex_m": int(m),
                        "n": c["n"], "turns": [c["p"], c["q"], c["r"]]})
    return out


def run(seed=0):
    gs = groups()
    rows = []
    for c in chains():
        n, turns = c["n"], c["turns"]
        want = sorted(set(turns) | {n - t for t in turns})
        matches = []
        for g in gs:
            if g["h"] != n:
                continue
            es = set(g["exponents"])
            if all(t in es or (n - t) in es for t in turns):
                matches.append({"group": g["group"], "rank": g["rank"],
                                "rank_equals_dimension": g["rank"] == c["dimension"],
                                "exponents_are_exactly_the_turns_and_complements":
                                    sorted(es) == want and len(g["exponents"]) == len(want)})
        rows.append(dict(c, turns_over_n=[f"{t}/{n}" for t in turns],
                         groups_with_h_equal_n=[g["group"] for g in gs if g["h"] == n],
                         matches=[m["group"] for m in matches],
                         matches_with_rank_equal_dimension=[m["group"] for m in matches
                                                            if m["rank_equals_dimension"]],
                         matches_exactly=[m["group"] for m in matches
                                          if m["exponents_are_exactly_the_turns_and_complements"]],
                         verdict="match" if matches else "no match"))
    return {"bound": {"seed": seed, "groups": "A2..A8, B2..B8, D4..D8, H3, H4, E6, E7, E8",
                      "not_in_the_list": ["F4", "G2", "I2(m)", "ranks above 8"],
                      "a_n_matches_anything": "A_(N-1) has every exponent 1..N-1 and h = N, so "
                                              "it matches any chain with N <= 9; that match "
                                              "is logged like the others",
                      "exponents": "computed from the Coxeter element's traces; h is its order "
                                   "(tolerance 1e-8)",
                      "match": "N = h and every turn p, or N - p, is an exponent"},
            "groups": gs, "chains": rows,
            "summary": {"chains": len(rows),
                        "match": sum(1 for r in rows if r["matches"]),
                        "no_match": [f"N={r['n']} {r['turns']} (m={r['simplex_m']})"
                                     for r in rows if not r["matches"]],
                        "match_with_rank_equal_dimension":
                            sum(1 for r in rows if r["matches_with_rank_equal_dimension"]),
                        "match_exactly": sum(1 for r in rows if r["matches_exactly"])}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(["m" if r["matches"] else "n" for r in out["chains"]],
                                 "one letter per chain: m match, n no match")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"sources": SOURCES, "max_rank": MAX_RANK,
                                     "status_before_run": "PROPOSED"},
                        out, None, None, comp, None)
    print(p)
    print("GROUPS", [(g["group"], g["h"], g["exponents"]) for g in out["groups"]])
    for r in out["chains"]:
        print(f"  {r['dimension']}D m={r['simplex_m']} N={r['n']:3d} turns {r['turns']}: "
              f"{r['verdict']}  all {r['matches']}  rank=dim "
              f"{r['matches_with_rank_equal_dimension']}  exact {r['matches_exactly']}")
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])

"""CC0. Co-authored. EXPLORE follow-up to Q5 (not pre-registered). Widths where the window width equals the perpendicular
coordinate of a lattice vector: r = |n - m/phi| / phi, |m|,|n| <= 4, 0.05 < r < 2.05.
Each tested at r*(1-1e-7), r, r*(1+1e-7). Status of any reading: PROPOSED.
  python3 experiments/q5_special_widths.py [--small]    (--small: |m|,|n| <= 2; own run file)"""
import sys, math
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / 'experiments'):
    sys.path.insert(0, str(_p))
import q5_q9_sweeps as R
from toolkit import loghook, composers as C, primitives as P
par, perp = C.line_basis_2d(R.SLOPE); lo, hi = C.strip_bounds_1d(perp[0]); canon = hi - lo
phi = P.PHI
cands = {}
SMALL = '--small' in sys.argv
MN = 2 if SMALL else 4
for m in range(-MN,MN+1):
    for n in range(-MN,MN+1):
        r = abs(n - m/phi)/phi
        if 0.05 < r < 2.05:
            cands.setdefault(round(r,12), (m,n))
rows = []
for r,(m,n) in sorted(cands.items()):
    rec = {"r": r, "m": m, "n": n}
    for tag,f in (("below",1-1e-7),("at",1.0),("above",1+1e-7)):
        npts, cl = R.gap_lengths(r*f*canon/2.0)
        tot = sum(c["count"] for c in cl); cnt = sorted(c["count"] for c in cl)
        rec[tag] = {"n_gap_lengths": len(cl), "min_gap_frequency": cnt[0]/tot}
    rows.append(rec)
two_at = [x["r"] for x in rows if x["at"]["n_gap_lengths"]==2]
two_any = [x["r"] for x in rows if any(x[t]["n_gap_lengths"]==2 for t in ("below","at","above"))]
outputs = {"canonical_width_q": canon, "rows": rows, "r_with_two_lengths_at": two_at,
           "r_with_two_lengths_on_any_side": two_any, "mode": "EXPLORE", "status": "PROPOSED"}
params = {"experiment": "q5_special_widths", "r_formula": "|n - m/phi|/phi", "mn_range": MN, "perturbation": 1e-7,
          "target_points": R.TARGET, "repo_commit": "3162c9950136a4a7f14e1e064f8d335db33968c8", "mode": "EXPLORE"}
print(loghook.log_run("q5_special_widths" + ("_small" if SMALL else ""), 0, params, outputs, out_dir=None))
print("tested", len(rows)); print("two at exactly r:", two_at); print("two on any side:", two_any)
for x in rows[:40]: print(x["r"], (x["m"],x["n"]), x["below"]["n_gap_lengths"], x["at"]["n_gap_lengths"], x["above"]["n_gap_lengths"], "minfreq_at %.4f"%x["at"]["min_gap_frequency"])

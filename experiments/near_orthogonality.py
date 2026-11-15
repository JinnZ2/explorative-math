"""a. near_orthogonality -- angle between random unit vectors vs dimension.

CC0. Co-authored. Stdlib only.

Draw pairs of independent random unit vectors (normalised Gaussians) in d
dimensions and record the angle between them. As d grows the angles pile up
near 90 degrees.

  python3 experiments/near_orthogonality.py [seed]
"""
import math
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402

NAME = "near_orthogonality"
DIMS = (2, 3, 5, 10, 20, 50, 100, 200)
PAIRS = 400
BAND = 5.0  # degrees either side of 90


def run(seed=0, dims=DIMS, pairs=PAIRS):
    rng = random.Random(seed)
    per_dim, angles = {}, {}
    for d in dims:
        a = []
        cosines = []
        for _ in range(pairs):
            u, v = emlib.random_unit(rng, d), emlib.random_unit(rng, d)
            c = emlib.dot(u, v)
            cosines.append(c)
            a.append(math.degrees(math.acos(max(-1.0, min(1.0, c)))))
        sd_cos = statistics.pstdev(cosines)
        per_dim[str(d)] = {
            "mean_angle": statistics.fmean(a),
            "sd_angle": statistics.pstdev(a),
            "mean_cos": statistics.fmean(cosines),
            "sd_cos": sd_cos,
            "sd_cos_times_sqrt_d": sd_cos * math.sqrt(d),
            "frac_within_band": sum(1 for x in a if abs(x - 90.0) <= BAND) / pairs,
        }
        angles[str(d)] = [round(x, 6) for x in a]
    return per_dim, angles


def main(argv):
    seed = int(argv[0]) if argv else 0
    per_dim, angles = run(seed)
    symbols = [chr(ord("a") + min(17, int(x // 10)))
               for d in DIMS for x in angles[str(d)]]
    comp = emlib.compressibility(
        symbols, "angles in 10-degree bins as letters a..r, dims in ascending order")
    ratios = [per_dim[str(d)]["sd_cos_times_sqrt_d"] for d in DIMS]
    ai = (f"OBSERVED: fraction within 90+/-{BAND:g} deg goes "
          f"{per_dim[str(DIMS[0])]['frac_within_band']:.3f} at d={DIMS[0]} to "
          f"{per_dim[str(DIMS[-1])]['frac_within_band']:.3f} at d={DIMS[-1]}; "
          f"sd(cos)*sqrt(d) stays in [{min(ratios):.3f}, {max(ratios):.3f}]. "
          "DERIVED: consistent with sd(cos) ~ 1/sqrt(d). AI reading: the d=2 and d=3 "
          "rows are the ones with structure (flat and sine-shaped densities); above "
          "d~20 the rows look alike.")
    p = emlib.write_run(NAME, seed,
                        {"dims": list(DIMS), "pairs": PAIRS, "band_deg": BAND,
                         "sampler": "normalised Gaussian"},
                        {"per_dim": per_dim}, angles, ai, comp, None)
    print(p)
    for d in DIMS:
        r = per_dim[str(d)]
        print(f"d={d:4d}  mean={r['mean_angle']:7.3f}  sd={r['sd_angle']:7.3f}  "
              f"sd_cos*sqrt(d)={r['sd_cos_times_sqrt_d']:.3f}  within={r['frac_within_band']:.3f}")


if __name__ == "__main__":
    main(sys.argv[1:])

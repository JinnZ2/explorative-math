"""ac. complement_fit -- how well a probe shape fits the complement of a baseline.

CC0. Co-authored. Stdlib only. Continues experiments/collapse_gap.py (its
steps 6, 7, 8). Status: PROPOSED. One concrete reading of a carried method;
the choices are below and in the run's bound.

  6 fit score      the baseline's pattern sits round the collapse point. A
                   probe shape is placed in what is left (its size and its
                   turn chosen to leave the least gap), and the RESIDUAL GAP
                   is read: the mean distance from a sample to the nearest
                   point of baseline + probe. Lower = the probe fills more
                   of the complement. Probes are ranked by it; the baseline
                   alone is the reference.
  7 combinations   pairs and triples of probe shapes placed together, one
                   after another (each placed to leave the least gap given
                   what is already there).
  8 environment    the same under a field/ gradient across the plane: the
                   local unit of length is the spacing a cell prefers at the
                   field value there (field/hexsheet.preferred_spacing), so
                   a gap counts by how many local units it spans. A uniform
                   zero field is the control. Fit is relative to the
                   environment; how the rankings shift is logged.

  sphere           the residual gap in full R^4: the deepest void left after
                   the probe's vertex directions are added to the
                   baseline's, in the frame the shapes are built in.

A probe with more points fills more whatever its shape; the gap closed per
probe point is logged beside the residual gap. Observer slots null.

  python3 experiments/complement_fit.py [seed]
  python3 experiments/complement_fit.py --small [seed]   # phone bound; own run file
"""
import itertools
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import collapse_gap as G  # noqa: E402
import emlib  # noqa: E402
from field import hexsheet as H  # noqa: E402
from toolkit import complement as CP  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "complement_fit"
N_RINGS, N_THETA = 6, 1680
ANGLE_STRIDE = 20            # samples every 20th grid angle: 84 per ring
TURN_STEP = 7                # probe turns every 7 grid steps: 1.5 degrees
SCALES = (0.5, 1.0)          # probe sizes tried (the baseline has radius 1)
BASELINES = ("5-cell", "16-cell")
COMBO_PROBES = ("3-gon", "5-gon", "6-gon", "8-gon", "12-gon", "tesseract")
ALPHA, SLOPE = 0.5, 1.0      # gradient: field value = SLOPE * x; spacing 1 + ALPHA * value
SPHERE_STARTS, SPHERE_STEPS = 80, 300
SMALL = {"n_theta": 840, "angle_stride": 20, "turn_step": 7, "combo_probes": ("3-gon", "5-gon", "8-gon"),
         "sphere_starts": 30, "sphere_steps": 200}
ENVIRONMENTS = ("none", "uniform_zero", "gradient")
METHOD = ("PROPOSED, carried: place a probe shape in the complement of the baseline and read "
          "the residual gap; rank the fits; combine probes; rerun under a field gradient and "
          "log how the rankings shift.")


def weights(radii, n_theta, stride, environment):
    """One weight per sample: 1 / (the local unit of length there)."""
    out = []
    for r in radii:
        row = []
        for j in range(0, n_theta, stride):
            x = r * math.cos(2 * math.pi * j / n_theta)
            if environment == "gradient":
                unit = H.preferred_spacing(SLOPE * x, ALPHA)
            else:
                unit = H.preferred_spacing(0.0, ALPHA)
            row.append(1.0 / unit)
        out.append(row)
    return out


def sample(field, n_theta, stride, shift=0):
    return [[row[(j + shift) % n_theta] for j in range(0, n_theta, stride)] for row in field]


def gap(samples, w):
    tot = cnt = 0.0
    for srow, wrow in zip(samples, w):
        for s, x in zip(srow, wrow):
            tot += s * x
            cnt += 1
    return tot / cnt


def deepest(samples, w):
    return max(s * x for srow, wrow in zip(samples, w) for s, x in zip(srow, wrow))


def place(current, probe_fields, n_theta, stride, turn_step, w):
    """Place one probe on what is already there: the size and turn that leave
    the smallest mean gap. Returns (gap, scale, turn in degrees, new samples)."""
    best = None
    for scale, pf in probe_fields.items():
        for shift in range(0, n_theta, turn_step):
            tot = cnt = 0.0
            for crow, prow, wrow in zip(current, pf, w):
                for k, j in enumerate(range(0, n_theta, stride)):
                    tot += min(crow[k], prow[(j + shift) % n_theta]) * wrow[k]
                    cnt += 1
            g = tot / cnt
            if best is None or g < best[0] - 1e-15:
                best = (g, scale, shift)
    g, scale, shift = best
    pf = probe_fields[scale]
    new = [[min(crow[k], prow[(j + shift) % n_theta])
            for k, j in enumerate(range(0, n_theta, stride))]
           for crow, prow in zip(current, pf)]
    return g, scale, 360.0 * shift / n_theta, new


def ranking(rows, key="residual_gap"):
    rows = sorted(rows, key=lambda r: (r[key], r["name"]))
    for i, r in enumerate(rows):
        r["rank"] = i + 1
    return rows


def plane_part(seed, n_theta, stride, turn_step, baselines, combo_probes):
    shapes = {n: CP.turned_to_zero(CP.centred_unit(p)) for n, _, p in G.plane_shapes(seed)}
    radii, _ = CP.polar_grid(N_RINGS, n_theta)
    fields = {n: {s: CP.distance_field([[s * x, s * y] for x, y in pat], radii, n_theta)
                  for s in SCALES} for n, pat in shapes.items()}
    out = []
    for b in baselines:
        per_env = {}
        base_samples = sample(fields[b][1.0], n_theta, stride)
        for env in ENVIRONMENTS:
            w = weights(radii, n_theta, stride, env)
            g0, d0 = gap(base_samples, w), deepest(base_samples, w)
            singles, cache = [], {}
            for p in shapes:                      # the baseline itself is a probe too ("self")
                g, scale, turn, new = place(base_samples, fields[p], n_theta, stride,
                                            turn_step, w)
                cache[(p,)] = (g, new)
                singles.append({"name": "self" if p == b else p, "shape": p,
                                "n_points": len(shapes[p]), "residual_gap": g,
                                "residual_deepest_void": deepest(new, w),
                                "scale": scale, "turn_deg": turn,
                                "gap_closed": g0 - g,
                                "gap_closed_per_probe_point": (g0 - g) / len(shapes[p])})
            singles = ranking(singles)
            entry = {"baseline_gap": g0, "baseline_deepest_void": d0, "singles": singles,
                     "ranking": [r["name"] for r in singles],
                     "ranking_by_gap_closed_per_point":
                         [r["name"] for r in sorted(singles,
                                                    key=lambda r: -r["gap_closed_per_probe_point"])]}
            if env != "uniform_zero":
                order = [r["shape"] for r in singles if r["shape"] in combo_probes]
                combos = {2: [], 3: []}
                for size in (2, 3):
                    for combo in itertools.combinations(order, size):   # best single first
                        steps = []
                        for k in range(1, size + 1):
                            pre = combo[:k]
                            if pre not in cache:
                                g, scale, turn, new = place(cache[pre[:-1]][1], fields[pre[-1]],
                                                            n_theta, stride, turn_step, w)
                                cache[pre] = (g, new, scale, turn)
                            c = cache[pre]
                            steps.append({"probe": pre[-1], "gap_after": c[0],
                                          "scale": c[2] if len(c) > 2 else next(
                                              r["scale"] for r in singles if r["shape"] == pre[0]),
                                          "turn_deg": c[3] if len(c) > 3 else next(
                                              r["turn_deg"] for r in singles
                                              if r["shape"] == pre[0])})
                        g, new = cache[combo][0], cache[combo][1]
                        npts = sum(len(shapes[p]) for p in combo)
                        combos[size].append({"name": " + ".join(combo), "probes": list(combo),
                                             "n_points": npts, "residual_gap": g,
                                             "residual_deepest_void": deepest(new, w),
                                             "gap_closed": g0 - g,
                                             "gap_closed_per_probe_point": (g0 - g) / npts,
                                             "placement": steps})
                entry["pairs"] = ranking(combos[2])
                entry["triples"] = ranking(combos[3])
                entry["pair_ranking"] = [r["name"] for r in entry["pairs"]]
                entry["triple_ranking"] = [r["name"] for r in entry["triples"]]
                entry["combination_probes"] = list(order)
            per_env[env] = entry

        def shift(a, c, key):
            ra, rc = per_env[a][key], per_env[c][key]
            pos = {n: i for i, n in enumerate(rc)}
            moved = [{"name": n, "rank_before": i + 1, "rank_after": pos[n] + 1,
                      "moved_by": i - pos[n]} for i, n in enumerate(ra) if n in pos]
            return {"from": a, "to": c, "same_order": ra == rc,
                    "n_changed_rank": sum(1 for m in moved if m["moved_by"]),
                    "largest_move": max((abs(m["moved_by"]) for m in moved), default=0),
                    "best_before": ra[0] if ra else None, "best_after": rc[0] if rc else None,
                    "moves": moved}
        out.append({"baseline": b, "environments": per_env,
                    "shift_none_to_uniform_zero": shift("none", "uniform_zero", "ranking"),
                    "shift_none_to_gradient": {
                        "singles": shift("none", "gradient", "ranking"),
                        "singles_per_point": shift("none", "gradient",
                                                   "ranking_by_gap_closed_per_point"),
                        "pairs": shift("none", "gradient", "pair_ranking"),
                        "triples": shift("none", "gradient", "triple_ranking")},
                    "gradient_placements_differ": [
                        {"name": a["name"], "turn_none": a["turn_deg"], "turn_gradient": c["turn_deg"],
                         "scale_none": a["scale"], "scale_gradient": c["scale"]}
                        for a in per_env["none"]["singles"]
                        for c in per_env["gradient"]["singles"]
                        if a["name"] == c["name"]
                        and (a["turn_deg"] != c["turn_deg"] or a["scale"] != c["scale"])]})
    return out


def sphere_part(seed, starts, steps, baselines):
    shapes = dict(G.sphere_shapes(seed))
    out = []
    for b in baselines:
        base = CP.sphere_voids(shapes[b], seed, starts, steps)
        rows = []
        for p, verts in shapes.items():
            u = CP.sphere_voids(shapes[b] + verts, seed, starts, steps)
            rows.append({"name": "self" if p == b else p, "n_points": len(verts),
                         "residual_gap": u["deepest_depth_deg"],
                         "gap_closed": base["deepest_depth_deg"] - u["deepest_depth_deg"],
                         "gap_closed_per_probe_point":
                             (base["deepest_depth_deg"] - u["deepest_depth_deg"]) / len(verts)})
        rows = ranking(rows)
        out.append({"baseline": b, "baseline_deepest_void_deg": base["deepest_depth_deg"],
                    "probes": rows, "ranking": [r["name"] for r in rows],
                    "bound": base["bound"]})
    return {"residual_gap": "deepest void (degrees) left after adding the probe's directions",
            "frame": "each shape as toolkit.primitives builds it; no rotation is searched",
            "baselines": out}


def run(seed=0, n_theta=N_THETA, stride=ANGLE_STRIDE, turn_step=TURN_STEP,
        baselines=BASELINES, combo_probes=COMBO_PROBES, sphere_starts=SPHERE_STARTS,
        sphere_steps=SPHERE_STEPS):
    return {"method": METHOD,
            "bound": {"seed": seed, "n_rings": N_RINGS, "n_theta": n_theta,
                      "samples_per_ring": n_theta // stride,
                      "probe_turn_step_deg": 360.0 * turn_step / n_theta,
                      "probe_scales": list(SCALES), "baselines": list(baselines),
                      "combination_probes": list(combo_probes),
                      "combinations": "placed one after another, best single first; not a "
                                      "joint search over all turns",
                      "residual_gap": "mean over samples of (distance to the nearest point of "
                                      "baseline + probes) / (local unit of length)",
                      "environment": {"gradient": f"field value = {SLOPE} * x; local unit = "
                                                  f"1 + {ALPHA} * value "
                                                  "(field/hexsheet.preferred_spacing, linear)",
                                      "uniform_zero": "local unit 1 everywhere (control)",
                                      "none": "no weighting"},
                      "patterns_do_not_move": "the field changes how a gap is counted, not "
                                              "where the points are",
                      "more_points_fill_more": "gap closed per probe point is logged beside "
                                               "the residual gap",
                      "sphere_starts": sphere_starts, "sphere_steps": sphere_steps},
            "plane": plane_part(seed, n_theta, stride, turn_step, baselines, combo_probes),
            "sphere": sphere_part(seed, sphere_starts, sphere_steps, baselines)}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["n_theta"], SMALL["angle_stride"], SMALL["turn_step"], BASELINES,
               SMALL["combo_probes"], SMALL["sphere_starts"], SMALL["sphere_steps"])
           if small else run(seed))
    symbols = []
    for b in out["plane"]:
        symbols += [chr(ord("a") + min(25, m["rank_after"] - 1))
                    for m in b["shift_none_to_gradient"]["singles"]["moves"]]
    comp = emlib.compressibility(
        symbols, "rank under the gradient (a = 1st) of each probe, listed in its no-field "
                 "rank order, baseline by baseline")
    bd = out["bound"]
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "n_theta": bd["n_theta"],
                         "samples_per_ring": bd["samples_per_ring"],
                         "probe_turn_step_deg": bd["probe_turn_step_deg"],
                         "probe_scales": list(SCALES), "baselines": bd["baselines"],
                         "combination_probes": bd["combination_probes"],
                         "alpha": ALPHA, "slope": SLOPE, "environments": list(ENVIRONMENTS),
                         "method": METHOD, "status_before_run": "PROPOSED"},
                        out, None, None, comp, None, gz=True)
    print(p)
    for b in out["plane"]:
        print("PLANE baseline", b["baseline"])
        for env in ("none", "gradient"):
            e = b["environments"][env]
            print(f"  [{env}] baseline gap {e['baseline_gap']:.4f}")
            print("    singles:", [(r["name"], round(r["residual_gap"], 4), r["scale"],
                                    round(r["turn_deg"], 1)) for r in e["singles"]])
            print("    per point:", e["ranking_by_gap_closed_per_point"])
            print("    pairs (first 4):", [(r["name"], round(r["residual_gap"], 4))
                                           for r in e["pairs"][:4]])
            print("    triples (first 3):", [(r["name"], round(r["residual_gap"], 4))
                                             for r in e["triples"][:3]])
        s = b["shift_none_to_gradient"]
        print("  uniform-zero control same order:", b["shift_none_to_uniform_zero"]["same_order"])
        for k in ("singles", "singles_per_point", "pairs", "triples"):
            print(f"  shift {k}: changed rank {s[k]['n_changed_rank']}, largest move "
                  f"{s[k]['largest_move']}, best {s[k]['best_before']} -> {s[k]['best_after']}")
        print("  placements that differ under the gradient:", len(b["gradient_placements_differ"]))
    for b in out["sphere"]["baselines"]:
        print("SPHERE baseline", b["baseline"], f"deepest void {b['baseline_deepest_void_deg']:.2f}")
        print("   ", [(r["name"], round(r["residual_gap"], 2)) for r in b["probes"]])


if __name__ == "__main__":
    main(sys.argv[1:])

"""t. field_response -- an ideal hexagonal sheet answering an environmental field.

CC0. Co-authored. Stdlib only. CONSTRUCTED toy (field/hexsheet.py): its
numbers are readings of this model, not findings about any material.

Status of the claims: PROPOSED before the run.
  P1  each cell's deviation from the ideal hexagon follows the field value
      where the cell sits: the deviation is a RESPONSE, not an error.
  P2  the mean over cells approximates the ideal hexagon, although single
      cells do not.
  P3  control: the same run with a uniform field gives deviations of about 0.

  gradient   field value = slope * x across the sheet (a heat gradient),
             mean zero. Each cell prefers a spacing 1 + alpha * value.
  control    the same sheet, rule, sweeps and rate, field uniformly 0.

Logged per cell: the field value and the response (size, shape, and each
neighbour offset minus the ideal offset). The probe correlates the size
response with the field value over interior cells and reads the mean cell.

  python3 experiments/field_response.py [seed]
  python3 experiments/field_response.py --small [seed]   # phone bound: smaller sheet, fewer
                                                         # sweeps, every output kept, own file
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from field import hexsheet as H  # noqa: E402

NAME = "field_response"
RADIUS = 6                # sheet of 127 cells
SLOPE = 0.5               # field units per ideal cell spacing
ALPHA = 0.01              # preferred spacing = 1 + ALPHA * field value
SWEEPS = 3000
RATE = 0.15
SMALL = {"radius": 3, "sweeps": 800}
CLAIMS = {"P1": "PROPOSED before the run: a cell's deviation from the ideal follows the field "
                "value where it sits (response, not error).",
          "P2": "PROPOSED before the run: the mean over cells approximates the ideal hexagon.",
          "P3": "PROPOSED before the run: with a uniform field the deviations are about 0."}


def one(sheet, field, alpha, sweeps, rate):
    pos, bound = H.relax(sheet, field, alpha, sweeps, rate)
    recs = H.responses(sheet, field, pos)
    inner = [r for r in recs if r["interior"]]
    fv = [r["field_value"] for r in inner]
    size = [r["size_response"] for r in inner]
    shape = [r["shape_response"] for r in inner]
    probe = {
        "over": "interior cells (six neighbours)", "n_interior_cells": len(inner),
        "n_edge_cells": len(recs) - len(inner),
        "correlation_size_response_vs_field": H.pearson(fv, size),
        "slope_size_response_vs_field": H.slope(fv, size),
        "alpha": alpha,
        "correlation_shape_response_vs_abs_field": H.pearson([abs(v) for v in fv], shape),
        "max_abs_size_response": max(abs(v) for v in size),
        "mean_size_response": sum(size) / len(size),
        "max_shape_response": max(shape),
        "max_abs_size_response_all_cells": max(abs(r["size_response"]) for r in recs),
        "mean_cell": H.mean_cell(recs)}
    return {"relaxation_bound": bound,
            "field_range": [min(field.values()), max(field.values())],
            "probe": probe, "cells": recs}


def run(seed=0, radius=RADIUS, sweeps=SWEEPS, slope=SLOPE, alpha=ALPHA, rate=RATE):
    sheet = H.cells(radius)
    grad = one(sheet, H.gradient_field(sheet, slope), alpha, sweeps, rate)
    ctrl = one(sheet, H.uniform_field(sheet, 0.0), alpha, sweeps, rate)
    return {"claims": CLAIMS, "constructed": True,
            "bound": {"radius": radius, "n_cells": len(sheet), "sweeps": sweeps, "rate": rate,
                      "slope": slope, "alpha": alpha,
                      "ideal": "axial integers; x exact Fraction, y exact Fraction * sqrt(3)/2",
                      "response_arithmetic": "binary float (64-bit)",
                      "boundary": "free edge; edge cells logged, kept out of the probe"},
            "ideal_offsets": H.ideal_offsets(),
            "gradient": grad, "control_uniform": ctrl}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed, SMALL["radius"], SMALL["sweeps"]) if small else run(seed)
    g, c = out["gradient"]["probe"], out["control_uniform"]["probe"]
    gb, cb = out["gradient"]["relaxation_bound"], out["control_uniform"]["relaxation_bound"]
    inner = [r for r in out["gradient"]["cells"] if r["interior"]]
    lo, hi = min(r["size_response"] for r in inner), max(r["size_response"] for r in inner)
    span = (hi - lo) or 1.0
    symbols = [chr(ord("a") + min(9, int(10 * (r["size_response"] - lo) / span))) for r in inner]
    comp = emlib.compressibility(
        symbols, "size response of each interior cell under the gradient, in ten equal "
                 "bins a..j, cells in sheet order")
    ai = (f"CONSTRUCTED toy; readings of the model only. OBSERVED, gradient field "
          f"({out['bound']['n_cells']} cells, {g['n_interior_cells']} interior, field from "
          f"{out['gradient']['field_range'][0]:.3g} to {out['gradient']['field_range'][1]:.3g}, "
          f"alpha {ALPHA}): correlation of size response with field value "
          f"{g['correlation_size_response_vs_field']:.6f}, slope "
          f"{g['slope_size_response_vs_field']:.6f} (alpha {g['alpha']}); largest size response "
          f"{g['max_abs_size_response']:.3e}, mean {g['mean_size_response']:.3e}; largest shape "
          f"response {g['max_shape_response']:.3e}. Mean cell: largest mean offset response "
          f"{g['mean_cell']['max_abs_mean_offset_response']:.3e}, against largest single-cell "
          f"offset response {g['mean_cell']['largest_single_cell_offset_response']:.3e}; mean "
          f"side lengths {[round(x, 6) for x in g['mean_cell']['mean_side_lengths']]}, mean of "
          f"each opposite pair "
          f"{[round(x, 6) for x in g['mean_cell']['mean_opposite_pair_lengths']]}. "
          f"DERIVED: the offset to the neighbour up the gradient crosses higher field than "
          f"the offset to the neighbour down it, by slope * alpha / 2 = {SLOPE * ALPHA / 2:g} "
          f"each way, so the six mean sides differ by that much while each opposite pair "
          f"averages to the ideal. "
          f"OBSERVED, control with uniform field: largest size response "
          f"{c['max_abs_size_response']:.3e}, largest shape response "
          f"{c['max_shape_response']:.3e}; correlation "
          f"{c['correlation_size_response_vs_field']} (void: the field has no spread). "
          f"Relaxation bound: {gb['sweeps']} sweeps at rate {gb['rate']}; net pull left "
          f"{gb['max_net_pull_left']:.3e} (gradient), {cb['max_net_pull_left']:.3e} (control). "
          f"Claims logged before the reading: {CLAIMS}. "
          "AI reading: not interpretable atm beyond these values.")
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "radius": out["bound"]["radius"],
                         "sweeps": out["bound"]["sweeps"], "rate": RATE, "slope": SLOPE,
                         "alpha": ALPHA, "field": "gradient: slope * x; control: uniform 0",
                         "claims": CLAIMS, "status_before_run": "PROPOSED",
                         "constructed": True},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])

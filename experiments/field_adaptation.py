"""v. field_adaptation -- inverse readout, and which sheet sits better in a field.

CC0. Co-authored. Stdlib only. CONSTRUCTED toy (field/hexsheet.py): its
numbers are readings of this model, not findings about any material.

Words: "ideal" is the reference sheet, not a target. A cell's difference
from it is its response, not an error. The optimum is relative to the
environment.

Claims, PROPOSED before the run:
  C1  the field can be read back from the cell offsets alone.
  C2  under a field the adapted sheet carries lower strain, net pull and
      energy than the ideal reference sheet; the ideal reference sheet is
      the better one only at zero field.
  C3  control: at uniform zero field the two sheets are equal.
  C4  with a rule that saturates at high field, the inverse readout fails
      somewhere; where is logged.

  inverse readout   input: the relaxed positions (every cell's neighbour
                    offsets) and the rule. The field is not passed in.
                    output: a field value per cell, set beside the true one.
                    Read twice under the saturating rule: with the true
                    rule inverted, and with the linear rule assumed.
  performance       the same field (the same rest lengths) applied to two
                    sheets: the ideal reference (rigid, regular) and the
                    adapted one (relaxed under that field). Also over a
                    sweep of field slopes, 0 included.
  rules             linear, and saturating (1 + alpha * s * tanh(value / s)).

  python3 experiments/field_adaptation.py [seed]
  python3 experiments/field_adaptation.py --small [seed]   # phone bound: smaller sheet, fewer
                                                           # sweeps and slopes; own run file
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from field import hexsheet as H  # noqa: E402

NAME = "field_adaptation"
RADIUS = 6
SLOPE = 0.5               # the main field: value = SLOPE * x
SLOPES = (0.0, 0.125, 0.25, 0.5, 1.0)     # performance sweep; SLOPE and 0 must be in it
ALPHA = 0.01
SATURATION = 1.5          # the saturating rule levels off at 1 +- ALPHA * SATURATION
SWEEPS = 2000
RATE = 0.15
FAIL_FRACTION = 0.10      # a readout "fails" at a cell if it is off by more than this
                          # fraction of the true field's full range, or is void
SMALL = {"radius": 3, "sweeps": 600, "slopes": (0.0, 0.5)}
CLAIMS = {
    "C1": "PROPOSED before the run: the field can be read back from the cell offsets alone.",
    "C2": "PROPOSED before the run: under a field the adapted sheet has lower strain, net pull "
          "and energy than the ideal reference sheet; the ideal reference sheet is the better "
          "one only at zero field.",
    "C3": "PROPOSED before the run: at uniform zero field the two sheets are equal.",
    "C4": "PROPOSED before the run: with a saturating rule the inverse readout fails "
          "somewhere; where is logged."}


def readout(sheet, field, pos, alpha, rule, saturation, assumed_rule=None):
    """The inverse readout beside the true field. The field is used only AFTER
    the reconstruction, to score it."""
    use = assumed_rule or rule
    rec = H.inverse_readout(sheet, pos, alpha, use, saturation if use == "saturating" else None)
    s = set(sheet)
    span = (max(field.values()) - min(field.values()))
    cells = []
    for c in sheet:
        interior = all((c[0] + dq, c[1] + dr) in s for dq, dr in H.NEIGHBOURS)
        r = rec[c]
        diff = (r - field[c]) if r is not None else None
        fails = (r is None) or (span > 0 and abs(diff) > FAIL_FRACTION * span)
        cells.append({"cell": list(c), "interior": interior, "true_field": field[c],
                      "recovered_field": r, "difference": diff,
                      "in_rule_range": r is not None, "fails": fails})

    def part(rows):
        ok = [x for x in rows if x["recovered_field"] is not None]
        t = [x["true_field"] for x in ok]
        v = [x["recovered_field"] for x in ok]
        return {"n_cells": len(rows), "n_recovered": len(ok),
                "n_void_outside_rule_range": len(rows) - len(ok),
                "correlation": H.pearson(t, v), "slope_recovered_vs_true": H.slope(t, v),
                "max_abs_difference": max((abs(x["difference"]) for x in ok), default=None),
                "mean_abs_difference": (sum(abs(x["difference"]) for x in ok) / len(ok)
                                        if ok else None),
                "n_failing": sum(1 for x in rows if x["fails"])}
    failing = [x for x in cells if x["fails"]]
    return {"rule_used_for_the_inversion": use, "rule_that_made_the_sheet": rule,
            "fail_threshold": FAIL_FRACTION * span, "true_field_range": span,
            "interior": part([x for x in cells if x["interior"]]),
            "edge": part([x for x in cells if not x["interior"]]),
            "all": part(cells),
            "where_it_fails": {
                "n": len(failing),
                "cells": [x["cell"] for x in failing],
                "interior": sum(1 for x in failing if x["interior"]),
                "edge": sum(1 for x in failing if not x["interior"]),
                "min_abs_true_field_among_failing":
                    min((abs(x["true_field"]) for x in failing), default=None),
                "max_abs_true_field_among_passing":
                    max((abs(x["true_field"]) for x in cells if not x["fails"]), default=None)},
            "cells": cells}


def two_sheets(sheet, field, alpha, sweeps, rate, rule, saturation):
    es = H.rest_lengths(sheet, field, alpha, rule, saturation)
    ideal_pos = {c: H.ideal_xy(c) for c in sheet}
    pos, bound = H.relax(sheet, field, alpha, sweeps, rate, rule, saturation)
    pi, pa = H.performance(sheet, es, ideal_pos), H.performance(sheet, es, pos)
    keys = ("energy", "max_abs_strain", "rms_strain", "max_net_pull", "mean_net_pull")
    return pos, {"relaxation_bound": bound,
                 "ideal_reference_sheet": pi, "adapted_sheet": pa,
                 "adapted_lower": {k: pa[k] < pi[k] for k in keys},
                 "equal": {k: pa[k] == pi[k] for k in keys},
                 "adapted_minus_ideal": {k: pa[k] - pi[k] for k in keys}}


def one_rule(sheet, rule, saturation, slopes, main_slope, alpha, sweeps, rate):
    sweep, main = [], None
    for sl in slopes:
        field = H.gradient_field(sheet, sl) if sl else H.uniform_field(sheet, 0.0)
        pos, perf = two_sheets(sheet, field, alpha, sweeps, rate, rule, saturation)
        sweep.append({"slope": sl, "field_range": [min(field.values()), max(field.values())],
                      "performance": perf})
        if sl == main_slope:
            main = {"slope": sl, "field": "gradient: slope * x",
                    "performance": perf,
                    "inverse_readout": readout(sheet, field, pos, alpha, rule, saturation)}
            if rule == "saturating":
                main["inverse_readout_linear_rule_assumed"] = readout(
                    sheet, field, pos, alpha, rule, saturation, assumed_rule="linear")
        if sl == 0:
            control = {"slope": 0.0, "field": "uniform zero", "performance": perf,
                       "inverse_readout": readout(sheet, field, pos, alpha, rule, saturation)}
    return {"rule": rule, "saturation": saturation, "main": main, "control_zero_field": control,
            "slope_sweep": sweep,
            "slopes_where_ideal_reference_is_not_worse":
                [r["slope"] for r in sweep
                 if r["performance"]["ideal_reference_sheet"]["energy"]
                 <= r["performance"]["adapted_sheet"]["energy"]]}


def run(seed=0, radius=RADIUS, sweeps=SWEEPS, slopes=SLOPES, main_slope=SLOPE, alpha=ALPHA,
        saturation=SATURATION, rate=RATE):
    if main_slope not in slopes or 0.0 not in slopes:
        raise ValueError("the slope sweep must contain the main slope and 0")
    sheet = H.cells(radius)
    return {"claims": CLAIMS, "constructed": True,
            "naming": {"ideal": "the reference sheet, not a target",
                       "response": "a cell's difference from the reference; not an error",
                       "optimum": "relative to the environment"},
            "bound": {"radius": radius, "n_cells": len(sheet), "sweeps": sweeps, "rate": rate,
                      "alpha": alpha, "saturation": saturation, "slopes": list(slopes),
                      "main_slope": main_slope, "fail_fraction": FAIL_FRACTION,
                      "arithmetic": "binary float (64-bit); the reference is exact",
                      "readout_input": "relaxed positions and the rule; the field is not "
                                       "passed to the readout",
                      "readout_estimator": "mean distance to the neighbours a cell has, "
                                           "inverted through the rule",
                      "boundary": "free edge"},
            "rules": {"linear": one_rule(sheet, "linear", None, slopes, main_slope, alpha,
                                         sweeps, rate),
                      "saturating": one_rule(sheet, "saturating", saturation, slopes,
                                             main_slope, alpha, sweeps, rate)}}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = (run(seed, SMALL["radius"], SMALL["sweeps"], SMALL["slopes"]) if small else run(seed))
    lin, sat = out["rules"]["linear"], out["rules"]["saturating"]
    symbols = ["f" if c["fails"] else "p"
               for c in sat["main"]["inverse_readout"]["cells"]]
    comp = emlib.compressibility(
        symbols, "saturating rule, true rule inverted: one letter per cell in sheet order, "
                 "f readout fails, p passes")

    def rd(r):
        return (f"interior r = {_f(r['interior']['correlation'])} "
                f"(max |diff| {_f(r['interior']['max_abs_difference'])}), edge r = "
                f"{_f(r['edge']['correlation'])} (max |diff| "
                f"{_f(r['edge']['max_abs_difference'])}), void {r['all']['n_void_outside_rule_range']}"
                f", failing {r['where_it_fails']['n']} ({r['where_it_fails']['interior']} "
                f"interior, {r['where_it_fails']['edge']} edge)")

    def pf(p):
        return (f"energy ideal reference {p['ideal_reference_sheet']['energy']:.3e} / adapted "
                f"{p['adapted_sheet']['energy']:.3e}; max |strain| "
                f"{p['ideal_reference_sheet']['max_abs_strain']:.3e} / "
                f"{p['adapted_sheet']['max_abs_strain']:.3e}; max net pull "
                f"{p['ideal_reference_sheet']['max_net_pull']:.3e} / "
                f"{p['adapted_sheet']['max_net_pull']:.3e}")
    ai = (f"CONSTRUCTED toy; readings of the model only. {out['bound']['n_cells']} cells, "
          f"field slope {out['bound']['main_slope']}, alpha {ALPHA}, saturation {SATURATION}. "
          f"OBSERVED, linear rule -- inverse readout: {rd(lin['main']['inverse_readout'])}. "
          f"Performance: {pf(lin['main']['performance'])}. Zero field: "
          f"{pf(lin['control_zero_field']['performance'])}; equal "
          f"{lin['control_zero_field']['performance']['equal']}. Slopes where the ideal "
          f"reference sheet is not worse: {lin['slopes_where_ideal_reference_is_not_worse']}. "
          f"OBSERVED, saturating rule -- inverse readout with the true rule: "
          f"{rd(sat['main']['inverse_readout'])}; with the linear rule assumed: "
          f"{rd(sat['main']['inverse_readout_linear_rule_assumed'])}. Performance: "
          f"{pf(sat['main']['performance'])}. Zero field: "
          f"{pf(sat['control_zero_field']['performance'])}. Slopes where the ideal reference "
          f"sheet is not worse: {sat['slopes_where_ideal_reference_is_not_worse']}. "
          "No reading is given: the run was asked for as a log, not an interpretation.")
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "radius": out["bound"]["radius"],
                         "sweeps": out["bound"]["sweeps"], "rate": RATE,
                         "slopes": out["bound"]["slopes"], "main_slope": SLOPE, "alpha": ALPHA,
                         "saturation": SATURATION, "fail_fraction": FAIL_FRACTION,
                         "claims": CLAIMS, "status_before_run": "PROPOSED",
                         "constructed": True},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


def _f(x):
    return "void" if x is None else f"{x:.6g}"


if __name__ == "__main__":
    main(sys.argv[1:])

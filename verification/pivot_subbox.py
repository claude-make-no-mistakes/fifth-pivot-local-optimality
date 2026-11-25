"""Direct rational bound on a neighborhood of the cubic witness."""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pivot"))
from region import BOX, Interval, Q, X, Y, Z, enclose, stages_and_denominators
from verify import BASE_INTERVAL, benchmark, require


def main():
    box = {
        X: Interval(Q("-0.6176001"), Q("-0.6175999")),
        Y: Interval(Q("-0.779073"), Q("-0.779072")),
        Z: Interval(Q("0.45322499"), Q("0.45322501")),
    }
    stages, denominators = stages_and_denominators()
    for var, span in box.items():
        require(BOX[var].low <= span.low <= span.high <= BOX[var].high, "inside reconstruction domain")
    for name, expr in denominators.items():
        span = enclose(expr, box)
        require(span.high < 0 or span.low > 0, f"nonzero {name} throughout sub-box")
    h = stages[-1]
    # H12 = H11, H21 = -H11 at every feasible point of the retained family.
    # Therefore its signed last pivot is H22 - H21*H12/H11 = H22 + H11.
    growth = enclose(h[0, 0] + h[1, 1], box)
    benchmark()
    require(0 < Q("4.132464") < growth.low <= growth.high < Q("4.132467") < BASE_INTERVAL[0],
            "strict lower and upper rational comparisons")
    witness = json.loads((HERE.parent / "pivot/witness.json").read_bytes())["witness"]
    for var, ends in ((X, [witness["x"]] * 2), (Y, witness["y_isolating_interval"]),
                      (Z, [witness["z"]] * 2)):
        require(box[var].low <= Q(ends[0]) <= Q(ends[1]) <= box[var].high, "contains witness")
    print("Every feasible point in the witness sub-box has 4.132464 < p5 < 4.132467 < g*.")


if __name__ == "__main__":
    main()

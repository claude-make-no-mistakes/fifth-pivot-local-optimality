"""Rational interval proof for a box of relaxed pivot matrices."""

import json
from dataclasses import dataclass

import sympy as sp

from verify import Q, SIGNED_ENTRIES, TIES, raw_matrix, require, schur_stages


@dataclass(frozen=True)
class Interval:
    low: sp.Rational
    high: sp.Rational

    def __add__(self, other):
        return Interval(self.low + other.low, self.high + other.high)

    def __mul__(self, other):
        products = (self.low * other.low, self.low * other.high, self.high * other.low, self.high * other.high)
        return Interval(min(products), max(products))

    def __pow__(self, power):
        if power < 0:
            require(self.low > 0 or self.high < 0, f"denominator contains zero: {self}")
            return Interval(1 / self.high, 1 / self.low) ** (-power)
        if power == 0:
            return Interval(Q(1), Q(1))
        ends = (self.low**power, self.high**power)
        low = 0 if power % 2 == 0 and self.low <= 0 <= self.high else min(ends)
        return Interval(low, max(ends))

    def as_strings(self):
        return list(map(str, (self.low, self.high)))


X, Y, Z = sp.symbols("x y z")
VARIABLES = (X, Y, Z)
BOX = {
    X: Interval(Q("-0.6177"), Q("-0.6174")),
    Y: Interval(Q("-0.7794"), Q("-0.7789")),
    Z: Interval(Q("0.45322490"), Q("0.45322501")),
}
LEFT_Y = Interval(Q("-0.778957"), Q("-0.778955"))
RIGHT_Y = Interval(Q("-0.779307"), Q("-0.779304"))


def enclose(expr, box):
    expr = sp.sympify(expr)
    if expr in box:
        return box[expr]
    if expr.is_Rational:
        return Interval(expr, expr)
    if expr.is_Add:
        result = Interval(Q(0), Q(0))
        for term in expr.args:
            result += enclose(term, box)
        return result
    if expr.is_Mul:
        result = Interval(Q(1), Q(1))
        for factor in expr.args:
            result *= enclose(factor, box)
        return result
    if expr.is_Pow and expr.exp.is_Integer:
        return enclose(expr.base, box) ** int(expr.exp)
    raise ValueError(f"Cannot enclose {expr}")


def point(value):
    return Interval(value, value)


def mean_value(expr, box):
    midpoint = {var: (span.low + span.high) / 2 for var, span in box.items()}
    result = point(expr.subs(midpoint))
    for var, span in box.items():
        displacement = Interval(span.low - midpoint[var], span.high - midpoint[var])
        result += enclose(sp.diff(expr, var), box) * displacement
    return result


def stages_and_denominators():
    x, y, z = VARIABLES
    a, b = 1 + x, 1 + z
    q = 2 - z * (x + z)
    t = -z / y
    v = z * (x + z) / a
    d, e = -1 + x * z, 1 - x * v
    k = (q - z * a - y * (z + v)) / (e - d)
    f, u = b * k, (b - q) / k
    matrix = sp.Matrix([
        [1, 1, t, -z, v],
        [x, a + z, u + x * t, -1, 1],
        [-1, z, q + u - t, -1, -1],
        [y, f + y, 1, 1, 1],
        [z, -1, 1, -1, 1],
    ])
    f_stage = sp.Matrix([
        [b, u, d, e],
        [b, q + u, -b, -1 + v],
        [f, b, 1 + y * z, 1 - y * v],
        [-b, 1 - z * t, -1 + z**2, 1 - z * v],
    ])
    g_stage = sp.Matrix([
        [q, -z * a, -q],
        [q, 1 + y * z - k * d, 1 - y * v - k * e],
        [1 - z * t + u, -q, 2 - (x + z) * v],
    ])
    g01, g20 = g_stage[0, 1], g_stage[2, 0]
    h = sp.Matrix([
        [g_stage[1, 1] - g01, g_stage[1, 2] + q],
        [-q - g20 * g01 / q, g_stage[2, 2] + g20],
    ])
    denominators = {"1+x": a, "1+z": b, "y": y, "e-d": e - d, "k": k, "q": q, "H11": h[0, 0]}
    return (matrix, f_stage, g_stage, h), denominators


def prove():
    stages, denominators = stages_and_denominators()
    matrix, _, _, h = stages
    original = raw_matrix(*VARIABLES)
    require((original - matrix).applyfunc(sp.cancel).is_zero_matrix, "same tie reconstruction")
    for simple, full in zip(stages, schur_stages(original, sp.cancel)):
        require((simple - full).applyfunc(sp.cancel).is_zero_matrix, "compact stage identity")
    for (i, j), value in SIGNED_ENTRIES.items():
        require(matrix[i, j] == value, "prescribed initial entry")
    for stage, ties in zip(stages[1:], TIES):
        for i, j, side in ties:
            if stage is h and (i, j) == (1, 0):
                continue
            require(sp.cancel(stage[i, j] - side * stage[0, 0]) == 0, "prescribed algebraic tie")

    checks = {}

    def check(name, bound, sign=1):
        require(bound.low > 0 if sign == 1 else bound.high < 0, f"{name}: {bound}")
        checks[name] = bound.as_strings()

    for name, expr in denominators.items():
        check(f"denominator {name}", enclose(expr, BOX), sign=-1 if name == "y" else 1)
    for k, stage in enumerate(stages[:3], start=1):
        for i in range(stage.rows):
            for j in range(stage.cols):
                for side in (-1, 1):
                    slack = stage[0, 0] - side * stage[i, j]
                    if sp.cancel(slack) != 0:
                        check(f"stage {k} entry {i+1},{j+1} side {side:+}", mean_value(slack, BOX))

    p = h[0, 0] + h[1, 0]
    delta = h[1, 1] - h[0, 0]
    growth = h[1, 1] + h[0, 0]
    check("p5 along p=0", mean_value(growth, BOX))
    p_x = enclose(sp.diff(p, X), BOX)
    p_y = enclose(sp.diff(p, Y), BOX)
    check("p_x", p_x)
    check("p_y", p_y)
    for name, expr in (("H22-H11", delta), ("p5", growth)):
        d_x = enclose(sp.diff(expr, X), BOX)
        d_y = enclose(sp.diff(expr, Y), BOX)
        check(f"implicit x derivative of {name}", d_x + point(Q(-1)) * d_y * p_x * p_y**(-1))

    for name, x, yspan, delta_sign in (
        ("left", BOX[X].low, LEFT_Y, -1), ("right", BOX[X].high, RIGHT_Y, 1),
    ):
        require(BOX[Y].low < yspan.low < yspan.high < BOX[Y].high, "nested y bracket")
        face = dict(BOX)
        face[X], face[Y] = point(x), yspan
        check(f"{name} face H22-H11", mean_value(delta, face), sign=delta_sign)
        for label, endpoint, sign in (("low", yspan.low, -1), ("high", yspan.high, 1)):
            endpoint_box = dict(face)
            endpoint_box[Y] = point(endpoint)
            check(f"{name} face p at {label} y", mean_value(p, endpoint_box), sign)

    return {
        "box": {str(var): span.as_strings() for var, span in BOX.items()},
        "left_y_bracket": LEFT_Y.as_strings(),
        "right_y_bracket": RIGHT_Y.as_strings(),
        "strict_interval_checks": checks,
    }


if __name__ == "__main__":
    print(json.dumps(prove(), indent=2))

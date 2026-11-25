"""Exact certificates for the relaxed fifth-pivot problem."""

import argparse
import json
from math import lcm
from pathlib import Path

import sympy as sp


T = sp.Symbol("t")
HERE = Path(__file__).resolve().parent
Q = sp.Rational
BASE_INTERVAL = (Q("4.132517078632472854223346853277"), Q("4.132517078632472854223346853278"))
SIGNED_ENTRIES = {
    (0, 0): 1, (0, 1): 1, (1, 3): -1, (1, 4): 1, (2, 0): -1,
    (2, 3): -1, (2, 4): -1, (3, 2): 1, (3, 3): 1, (3, 4): 1,
    (4, 1): -1, (4, 2): 1, (4, 3): -1, (4, 4): 1,
}
TIES = (
    ((1, 0, 1), (1, 2, -1), (2, 1, 1), (3, 0, -1)),
    ((0, 2, -1), (1, 0, 1), (2, 1, -1)),
    ((0, 1, 1), (1, 0, -1)),
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def raw_matrix(x, y, z, reduce=sp.cancel):
    """Reconstruct (2.2)–(2.9) by solving the ties sequentially."""
    x, y, z = map(sp.sympify, (x, y, z))
    b = 1 + z
    q = 2 - z * (x + z)
    t = reduce(-z / y)
    v = reduce(z * (x + z) / (1 + x))
    d = -1 + x * z
    e = 1 - x * v
    f = reduce(b * (q - b - d - y * (z + v)) / (e - d))
    u = reduce((b - q) * b / f)
    return sp.Matrix([
        [1, 1, t, -z, v],
        [x, 1 + x + z, u + x * t, -1, 1],
        [-1, z, q + u - t, -1, -1],
        [y, f + y, 1, 1, 1],
        [z, -1, 1, -1, 1],
    ])


def schur_stages(matrix, reduce):
    stages = [matrix.applyfunc(reduce)]
    while stages[-1].rows > 1:
        m = stages[-1]
        stages.append((m[1:, 1:] - m[1:, :1] * m[:1, 1:] / m[0, 0]).applyfunc(reduce))
    return stages


def determinant_stage(matrix, k, reduce):
    """Use bordered determinants of the original matrix, with no Schur recursion."""
    head = list(range(k))
    den = reduce(matrix.extract(head, head).det(method="berkowitz"))
    return sp.Matrix(
        5 - k, 5 - k,
        lambda i, j: reduce(matrix.extract(head + [k + i], head + [k + j]).det(method="berkowitz") / den),
    )


def poly_interval(poly, low, high):
    left = right = Q(0)
    for coefficient in poly.all_coeffs():
        products = (left * low, left * high, right * low, right * high)
        left, right = min(products) + coefficient, max(products) + coefficient
    return left, right


class IsolatedRoot:
    def __init__(self, poly, low, high):
        self.poly = poly.monic()
        require(low < high, "nonempty root interval")
        require(self.poly.eval(low) != 0 and self.poly.eval(high) != 0, "root away from endpoints")
        require(self.poly.count_roots(low, high) == 1, "unique root in specified interval")
        require(self.poly.is_irreducible, "one irreducible polynomial for the algebraic field")
        self.low, self.high = low, high

    def reduce(self, expr):
        num, den = sp.cancel(expr).as_numer_denom()
        num, den = (sp.Poly(term, T, domain=sp.QQ) for term in (num, den))
        den = den.rem(self.poly)
        require(not den.is_zero, "defined algebraic denominator")
        return (num.rem(self.poly) * den.invert(self.poly)).rem(self.poly).as_expr()

    def interval(self, expr, width=Q(1, 10**20)):
        reduced = sp.Poly(self.reduce(expr), T, domain=sp.QQ)
        low, high = self.low, self.high
        while True:
            left, right = poly_interval(reduced, low, high)
            if right - left <= width:
                return left, right
            low, high = self.poly.refine_root(low, high, eps=(high - low) / 8)

    def sign(self, expr):
        reduced = self.reduce(expr)
        if reduced == 0:
            return 0
        width = Q(1, 10**12)
        while True:
            low, high = self.interval(reduced, width)
            if low > 0:
                return 1
            if high < 0:
                return -1
            width /= 100


def verify(matrix, reduce, sign, boundary=False):
    require(matrix.shape == (5, 5), "matrix size")
    for (i, j), value in SIGNED_ENTRIES.items():
        require(sign(matrix[i, j] - value) == 0, f"initial sign at ({i+1}, {j+1})")
    stages = schur_stages(matrix, reduce)
    active = []
    for k, m in enumerate(stages):
        direct = determinant_stage(matrix, k, reduce)
        pivot_sign = sign(m[0, 0])
        require(pivot_sign != 0, f"nonzero pivot {k+1}")
        absolute_pivot = pivot_sign * m[0, 0]
        stage_active = []
        for i in range(m.rows):
            for j in range(m.cols):
                require(sign(m[i, j] - direct[i, j]) == 0, f"determinant identity at stage {k+1}")
                for side in (-1, 1):
                    slack_sign = sign(absolute_pivot - side * m[i, j])
                    require(slack_sign >= 0, f"pivot dominance at stage {k+1}, ({i+1}, {j+1})")
                    if slack_sign == 0:
                        stage_active.append([i + 1, j + 1, side])
        active.append(stage_active)
    for k, ties in enumerate(TIES, start=1):
        for i, j, multiple in ties:
            require(sign(stages[k][i, j] - multiple * stages[k][0, 0]) == 0, f"prescribed tie at stage {k+1}")
    h = stages[3]
    ratio = reduce(h[1, 1] / h[0, 0])
    require(sign(ratio + 1) >= 0 and sign(ratio - 1) == (0 if boundary else -1), "stage-4 ratio")
    return stages, ratio, active


def descartes_variations(poly, low, high):
    """Count sign changes after mapping the positive line onto (low, high)."""
    low, high = map(Q, (low, high))
    common = lcm(int(low.q), int(high.q))
    numerator = sp.Poly(int(high * common) * T + int(low * common), T, domain=sp.ZZ)
    denominator = sp.Poly(common * (T + 1), T, domain=sp.ZZ)
    transformed = sp.Poly(0, T, domain=sp.ZZ)
    multiplier = sp.Poly(1, T, domain=sp.ZZ)
    for coefficient in poly.all_coeffs():
        transformed = numerator * transformed + coefficient * multiplier
        multiplier *= denominator
    signs = [sp.sign(c) for c in transformed.all_coeffs() if c]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def benchmark():
    coeffs = [int(line) for line in (HERE / "p5_coefficients.txt").read_text().splitlines()]
    poly = sp.Poly.from_list(coeffs[::-1], T, domain=sp.ZZ)
    low, high = BASE_INTERVAL
    require(poly.degree() == 61, "benchmark degree")
    for a, b in ((4, 5), (low, high)):
        require(descartes_variations(poly, a, b) == 1, "benchmark isolation")
        require(poly.eval(a) * poly.eval(b) < 0, "opposite signs at benchmark endpoints")
    require(descartes_variations(poly, -5, 0) == 0, "negative benchmark roots")
    require(poly.eval(-5) and poly.eval(0), "negative benchmark endpoints")
    return {"degree": 61, "roots_in_[4,5]": 1, "roots_in_[-5,0]": 0, "growth_interval": list(map(str, (low, high)))}


def controls():
    matrix = raw_matrix(-Q(1, 2), -Q(1, 2), Q(0))
    stages, _, active = verify(matrix, sp.cancel, sp.sign, boundary=True)
    pivots = [str(m[0, 0]) for m in stages]
    require(pivots == ["1", "1", "2", "2", "4"], "exact boundary fixture pivots")
    rejected = []
    fixtures = (
        (matrix, False, "stage-4 ratio"),
        (matrix.copy(), True, "initial sign"),
        (matrix.copy(), True, "pivot dominance at stage 1"),
    )
    fixtures[1][0][0, 1] = 0
    fixtures[2][0][0, 2] = 2
    for trial, boundary, expected in fixtures:
        try:
            verify(trial, sp.cancel, sp.sign, boundary)
        except ValueError as error:
            require(str(error).startswith(expected), f"wrong control failure: {error}")
            rejected.append(str(error))
        else:
            raise ValueError(f"unexpectedly accepted {expected} control")
    return {"boundary_pivots": pivots, "active_signed_entries": active, "rejected_controls": rejected}


def witness(x, z, low, high):
    raw = raw_matrix(x, T, z)
    stages = schur_stages(raw, sp.cancel)
    h = stages[3]
    num, _ = sp.cancel(h[0, 0] + h[1, 0]).as_numer_denom()
    poly = sp.Poly(num, T, domain=sp.QQ).clear_denoms()[1].primitive()[1]
    root = IsolatedRoot(poly, low, high)
    matrix = raw_matrix(x, T, z, root.reduce)
    stages, ratio, active = verify(matrix, root.reduce, root.sign)
    growth = root.sign(stages[-1][0, 0]) * stages[-1][0, 0]
    bounds = root.interval(growth)
    if bounds[0] > BASE_INTERVAL[1]:
        comparison = "above"
    elif bounds[1] < BASE_INTERVAL[0]:
        comparison = "below"
    else:
        comparison = "unresolved"
    return {
        "x": str(x), "z": str(z),
        "y_polynomial_descending": list(map(str, poly.all_coeffs())),
        "y_isolating_interval": list(map(str, (low, high))),
        "pivot_intervals": [list(map(str, root.interval(m[0, 0]))) for m in stages],
        "ratio_interval": list(map(str, root.interval(ratio))),
        "growth_interval": list(map(str, bounds)),
        "growth_vs_benchmark": comparison,
        "active_signed_entries": active,
        "matrix_intervals": [[list(map(str, root.interval(matrix[i, j]))) for j in range(5)] for i in range(5)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witness", metavar=("X", "Z", "Y_LOW", "Y_HIGH"), type=Q, nargs=4)
    args = parser.parse_args()
    certificate = {"sympy": sp.__version__, "benchmark": benchmark(), "controls": controls()}
    if args.witness:
        certificate["witness"] = witness(*args.witness)
    print(json.dumps(certificate, indent=2))


if __name__ == "__main__":
    main()

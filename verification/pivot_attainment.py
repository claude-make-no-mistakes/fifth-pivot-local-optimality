"""Certify the comparisons that force a local fifth-pivot maximum at g*."""

import sys
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pivot"))
from comparator import prove as certify_comparator
from verify import Q, require


def main():
    certificate = certify_comparator()
    g, z = sp.symbols("g z")
    coefficients = certificate["coefficients_ascending_z_g"]
    eliminant = sum(value * z**i * g**j for i, row in enumerate(coefficients)
                    for j, value in enumerate(row) if value)
    lower, upper = Q("4.13"), Q("4.14")
    leading = sp.Poly(sp.Poly(eliminant, z).LC(), g)
    factored = (g-8) * (g-6) * (g-4)**2 * (g-2) * (g**2-12*g+40)
    require(leading == sp.Poly(factored, g) and leading.count_roots(lower, upper) == 0,
            "constant degree in z on the whole comparator range")
    cut = Q("4.13251707863247283")
    actual_lower, actual_upper = map(Q, certificate["box_growth_interval"])
    require(lower < actual_lower <= actual_upper < upper, "all comparator pivots lie in (4.13,4.14)")
    left, right = map(Q, certificate["z_strip"])
    middle = Q("0.45322491")
    require(left < middle < right, "test value is interior to the z interval")

    for parameter, preferred in ((left, (lower, cut)), (right, (lower, cut)),
                                 (middle, (cut, upper))):
        fiber = sp.Poly(eliminant.subs(z, parameter), g, domain=sp.QQ)
        require(all(fiber.eval(point) != 0 for point in (lower, cut, upper)),
                f"nonzero interval endpoints at z={parameter}")
        require(fiber.count_roots(lower, upper) == 1 and fiber.count_roots(*preferred) == 1,
                f"one comparator root in {preferred} at z={parameter}")

    print("Both endpoint comparators are below the interior comparator by exact Sturm counts.")


if __name__ == "__main__":
    main()

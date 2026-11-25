"""Eliminate x and y to verify the local full-equality comparator bound."""

import json

import sympy as sp

from region import BOX, X, Y, Z, enclose, mean_value, stages_and_denominators
from verify import BASE_INTERVAL, HERE, Q, benchmark, descartes_variations, require


def prove():
    g = sp.Symbol("g")
    x, y, z = X, Y, Z
    (_, _, third, fourth), denominators = stages_and_denominators()
    a, b, q = 1 + x, 1 + z, third[0, 0]
    s, w, c = -third[0, 1], third[2, 0], third[2, 2]
    k = denominators["k"]
    h = fourth[0, 0]

    # At full equality h=g/2, w=g/2-c, and -h=H21=-q+sw/q.
    require(sp.cancel(fourth[1, 0] + q - s * w / q) == 0, "H21 identity")
    require(sp.cancel(fourth[1, 1] - c - w) == 0, "H22 identity")
    require(sp.cancel(q + s - (2 + z - z**2)) == 0, "coefficient of g/2")
    require(sp.cancel(q**2 + s * c - (4 + 2 * z * (1 - x - 2 * z))) == 0, "affine x identity")
    require(BOX[z].low > 0, "safe x recovery")
    recovered_x = 1 - 2 * z + 2 / z - g * (2 + z - z**2) / (4 * z)

    # H11 is affine in y once x has been fixed.
    beta = sp.diff(h, y)
    alpha = h.subs(y, 0)
    require(y not in beta.free_symbols and sp.cancel(h - alpha - beta * y) == 0, "affine y identity")
    beta_bounds = enclose(beta, BOX)
    require(beta_bounds.low > 0, "safe y recovery")
    recovered_y = sp.factor(((g / 2 - alpha) / beta).subs(x, recovered_x))
    require(sp.cancel(w - 1 - z**2 / y - (b - q) / k) == 0, "G31 identity")

    # Multiply w-(g/2-c) by a common denominator before substituting.
    # The first term uses the explicit G31 formula so each denominator is visible.
    product = (1 + c - g / 2) * a * k * y + z**2 * a * k + (b - q) * a * y
    require(sp.cancel(product - a * k * y * (w - (g / 2 - c))) == 0, "cleared G31 identity")
    recovered_a, recovered_c, recovered_q = (sp.factor(expr.subs(x, recovered_x)) for expr in (a, c, q))
    recovered_k = sp.factor(k.subs(x, recovered_x).subs(y, recovered_y))
    cleared = -1024 * z**4 * (z - 1)**2 * recovered_a * (
        (1 + recovered_c - g / 2) * recovered_k * recovered_y
        + z**2 * recovered_k + (b - recovered_q) * recovered_y
    )
    eliminant = sp.Poly(sp.cancel(cleared), z, g, domain=sp.ZZ)
    require(eliminant.degree_list() == (10, 7), "eliminant degrees")

    coeffs = [int(line) for line in (HERE / "p5_coefficients.txt").read_text().splitlines()]
    p5 = sp.Poly.from_list(coeffs[::-1], g, domain=sp.ZZ).as_expr()
    disc = eliminant.discriminant()
    expected = (
        -295147905179352825856 * (g - 6) * (g - 4)**32 * (g - 2)**6
        * (g**2 - 12 * g + 40) * p5
    )
    require(sp.Poly(disc.as_expr() - expected, g).is_zero, "discriminant identity")
    root_info = benchmark()

    g_low, g_high = BASE_INTERVAL[0], Q("4.14")
    z_low, z_high = BOX[z].low, BOX[z].high
    require(4 < g_low < g_high < 5, "critical g interval")
    growth_bounds = mean_value(fourth[0, 0] + fourth[1, 1], BOX)
    require(growth_bounds.high < g_high, "growth below top fiber")
    f = eliminant.as_expr()
    checks = {}

    def no_roots(name, expr, variable, low, high):
        polynomial = sp.Poly(expr, variable, domain=sp.QQ).clear_denoms()[1]
        variations = descartes_variations(polynomial, low, high)
        endpoint_signs = [int(sp.sign(polynomial.eval(endpoint))) for endpoint in (low, high)]
        require(variations == 0 and all(endpoint_signs), f"no roots on closed interval: {name}")
        checks[name] = {"descartes_variations": variations, "endpoint_signs": endpoint_signs}
        return polynomial

    no_roots("left z edge", f.subs(z, z_low), g, g_low, g_high)
    no_roots("right z edge", f.subs(z, z_high), g, g_low, g_high)
    leading = sp.Poly(f, z).LC()
    no_roots("leading coefficient in z", leading, g, g_low, g_high)
    upper = no_roots("top g fiber", f.subs(g, g_high), z, z_low, z_high)
    require(upper.count_roots(z_low, z_high) == 0, "Sturm cross-check of top fiber")
    checks["top g fiber"]["sturm_roots"] = 0

    return {
        "degrees_z_g": list(eliminant.degree_list()),
        "coefficients_ascending_z_g": [
            [int(eliminant.coeff_monomial(z**i * g**j)) for j in range(8)] for i in range(11)
        ],
        "benchmark": root_info,
        "discriminant_matches_published_factorization": True,
        "coefficient_of_y_interval": beta_bounds.as_strings(),
        "box_growth_interval": growth_bounds.as_strings(),
        "g_strip": list(map(str, (g_low, g_high))),
        "z_strip": list(map(str, (z_low, z_high))),
        "closed_interval_root_checks": checks,
    }


if __name__ == "__main__":
    print(json.dumps(prove(), indent=2))

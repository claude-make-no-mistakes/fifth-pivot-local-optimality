# Formal Verification Status

Lean 4.34.1 (pinned in [`lean-toolchain`](lean-toolchain)) verifies two standalone files with zero `sorry`s and only standard Lean axioms (`propext`, `Classical.choice`, `Quot.sound`):

- [`formalization/RatioBound.lean`](formalization/RatioBound.lean): checks the rational order-arithmetic deductions behind the global necessary product inequality (`Pivot.necessary_product`) and the strict comparator bound implication (`Pivot.local_strict`).
- [`formalization/BenchmarkSigns.lean`](formalization/BenchmarkSigns.lean): recomputes by kernel reduction (`by decide`) all exact integer Möbius-transformed sign variations and endpoint signs for:
  - the degree-61 benchmark brackets $[4, 5]$ (`bracket_variation`, `bracket_crossing`) and $[g_L, g_H]$ (`narrow_variation`, `narrow_crossing`);
  - the comparator rectangle face exclusions, leading-coefficient exclusion, and top-fiber exclusion (`comparator_faces`);
  - the fixed-$z$ attainment fiber crossings at $z = 0.45322490$, $0.45322491$, and $0.45322501$ (`attainment_fibers`).

## Running the Lean checks

From the repository root:

```sh
lean formalization/RatioBound.lean
lean formalization/BenchmarkSigns.lean
```

## Binding and scope

[`verification/check_lean_inputs.py`](verification/check_lean_inputs.py) (invoked automatically by [`replay.py`](replay.py)) checks that `benchmarkCoeffs` and `eliminantCoeffs` in `formalization/BenchmarkSigns.lean` match `pivot/p5_coefficients.txt` and `pivot/comparator.json` identically.

The symbolic derivation of the eliminant $E(g, z)$ and its discriminant, exact rational interval arithmetic over the box $B$, Descartes' rule bridging sign variations to real roots, and the real-variable calculus steps (implicit function theorem and compactness) are verified by the Python/SymPy replay and the proofs in [`paper/paper.tex`](paper/paper.tex).

# Formal verification status

Lean 4.34.1 checks [the rational order-arithmetic steps](formalization/RatioBound.lean) behind the ratio bound and the strict comparator implication.

Run this command from the repository root with the toolchain specified in the toolchain file.

~~~sh
lean formalization/RatioBound.lean
lean formalization/BenchmarkSigns.lean
~~~

Lean also recomputes the integer sign variations and endpoint signs for the benchmark brackets, the comparator faces and the attainment fibers. The Python replay independently derives the eliminant and its discriminant. Applying Descartes' rule to the signs and proving the real-variable comparison and attainment remain outside Lean.

# Formal verification status

Lean 4.34.1 checks [the rational order-arithmetic steps](formalization/RatioBound.lean) behind the ratio bound and the strict comparator implication.

Run this command from the repository root with the toolchain specified in the toolchain file.

~~~sh
lean formalization/RatioBound.lean
~~~

The paper proves the normalized elimination identities, the real-variable comparison, and the attainment argument. The exact Python replay checks the interval and polynomial certificates, including the discriminant and root counts. The real-variable theorem and the certificates remain outside Lean.

# A Sharp Local Bound for the Fifth Complete Pivot

We study a one-tie relaxation of the extremal complete-pivoting pattern of Chen, Edelman, and Urschel (IMA Journal of Numerical Analysis, 2026; [doi:10.1093/imanum/drag039](https://doi.org/10.1093/imanum/drag039), [arXiv:2602.20390v2](https://arxiv.org/abs/2602.20390v2)) for real $5 \times 5$ matrices. Retaining all fourteen prescribed initial $\pm 1$ entries and all stage-2, stage-3, and off-diagonal stage-4 ties while relaxing the final diagonal tie $H_{22} = H_{11}$ to $-1 \le r = H_{22}/H_{11} < 1$, we prove:

- **Sharp local bound on a rational box:** On $B = [-0.6177, -0.6174] \times [-0.7794, -0.7789] \times [0.45322490, 0.45322501]$ for $(x, y, z) = (A_{21}, A_{41}, A_{51})$, every relaxed matrix satisfies $0 < p_5 < g_*$, where $g_* \approx 4.13251707863247$ is the unique root in $[4, 5]$ of the degree-61 benchmark polynomial $P_5$. The supremum over the relaxed family in $B$ equals $g_*$, and $g_*$ is attained in $B$ when $r = 1$ is allowed.
- **Certified relaxed witness:** An explicit algebraic witness in $B$ at $x = -386/625$, $z = 18129/40000$, and a cubic irrational $y \in (-0.779073, -0.779072)$ has $r \approx 0.999939320464152$ and $p_5 \approx 4.1324652818749806$.
- **Global necessary condition:** Across the entire relaxed family (including outside $B$), every matrix has $p_5 > 0$, and any counterexample with $p_5 > g_*$ must satisfy $0 < z < 1$ and $(1 + r)(2 + z(1 - z)) > g_*$ (in particular $r > 4g_*/9 - 1 \approx 0.83667$ and $z(1 - z) > g_*/2 - 2 \approx 0.06626$).

The region outside $B$ meeting the necessary condition, as well as the general $5 \times 5$ growth-factor conjecture, remains open in this work.

## Documents and Data

- [Paper (PDF)](paper/paper.pdf) · [LaTeX source](paper/paper.tex)
- [Proof notes](pivot/RESULTS.md) · [Sharpness and attainment notes](verification/PIVOT-ATTAINMENT.md) · [Frozen target specification](pivot/TARGET.md)
- [Formal verification status](FORMALIZATION.md)
- Archived exact certificates: [baseline.json](pivot/baseline.json), [witness.json](pivot/witness.json), [region.json](pivot/region.json), [comparator.json](pivot/comparator.json), [p5_coefficients.txt](pivot/p5_coefficients.txt)

## Reproducing the Verification

### Exact Python/SymPy certificates

Requires Python 3.12 and the pinned dependencies in `requirements.txt` (`sympy==1.14.0`, `mpmath==1.3.0`):

```sh
python3 -m pip install -r requirements.txt
python3 -B replay.py
```

`replay.py` verifies the SHA-256 hashes of the frozen inputs (`verification/FROZEN.sha256`), regenerates and compares all four archived JSON certificates (`pivot/verify.py`, `pivot/region.py`, `pivot/comparator.py`), checks that the Lean polynomial constants match the benchmark and derived eliminant (`verification/check_lean_inputs.py`), checks the direct rational enclosure on the witness sub-box (`verification/pivot_subbox.py`), and checks the exact Sturm counts for local attainment (`verification/pivot_attainment.py`).

### Lean 4 formal checks

With Lean 4.34.1 (specified in `lean-toolchain`):

```sh
lean formalization/RatioBound.lean
lean formalization/BenchmarkSigns.lean
```

See [FORMALIZATION.md](FORMALIZATION.md) for the exact scope of the Lean checks.

# Verification record

An independent model reviewer checked the mathematical arguments and replayed the exact certificates.

The reviewer matched the 62 integer coefficients of the benchmark polynomial with equation (2.15), independently isolated the witness's cubic root by a Sturm count, and verified its pivot inequalities and quoted decimal interval. The proof of the necessary condition \(r>4g_*/9-1\) applies throughout the specified relaxation, including any points missed by the local parametrization. [Pivot results](pivot/RESULTS.md) · [Benchmark polynomial in the paper](https://arxiv.org/html/2602.20390v2#S2.SS1)

The first box proof used the published Theorem 2.2 for the full-equality comparator. The reviewer found that its printed case discussion skips a displayed factor. The revised proof derives a local eliminant from the Schur identities and certifies the comparator bound directly. The reviewer rederived the formulas for recovering \(x,y\), verified the discriminant identity by 127 exact Sylvester determinant evaluations against a degree bound of 126, and independently reproduced the face and top-fiber root exclusions. The rendered equation (2.13) has an extra \(97792g\) at \(z^0\). The revised proof uses the independently derived polynomial and checks its discriminant against the benchmark. [Elimination proof](pivot/RESULTS.md) · [Printed equations](https://arxiv.org/html/2602.20390v2#S2.SS1)

The reviewer independently reproduced the fixed-\(z\) root comparisons that force the maximum of the comparator curve into the interior of the box. The leading coefficient stays nonzero throughout the comparator range, so the discriminant identifies the maximum as \(g_*\). Approaching it from \(r<1\) gives the supremum of the relaxed family. A separate direct interval calculation, also reviewed, excludes the smaller box containing the cubic witness using the Schur formulas alone. [Sharpness proof](verification/PIVOT-ATTAINMENT.md) · [Small-box enclosure](verification/pivot_subbox.py)

The literature checks did not establish publication priority.

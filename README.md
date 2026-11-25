# A Sharp Local Bound for the Fifth Complete Pivot

We relax one equality in the extremal pattern of Chen, Edelman and Urschel. Over a specified rational box the fifth pivot has supremum \(g_*\), the degree-61 benchmark, and attains it when the original tie is allowed. Exact interval and polynomial root certificates establish the result. A necessary inequality restricts challengers throughout the relaxed family.

[Paper PDF](paper/paper.pdf) · [LaTeX source](paper/paper.tex) · [Proof notes](pivot/RESULTS.md) · [Verification data](pivot/comparator.json)

## Replay the certificates

Use Python 3.12 and SymPy 1.14.0.

```sh
python -m pip install -r requirements.txt
python -B replay.py
```

The general conjecture remains open in this work. See [formal verification status](FORMALIZATION.md).

[Verification record](REVIEW.md) · [Draft post for X](POST.md)

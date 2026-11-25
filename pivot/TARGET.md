# Frozen pivot target

Use the real, normalized, completely pivoted matrices of equation (2.1). Keep all 14 signed initial entries, the four stage-2 ties, the three stage-3 ties, and the stage-4 ties \(H_{12}=H_{11}\) and \(H_{21}=-H_{11}\). Replace \(H_{22}=H_{11}\) with \(-1\le r=H_{22}/H_{11}<1\).

A counterexample requires every pivot to be nonzero, every stage entry to satisfy pivot dominance in the given elimination order, every retained tie to hold exactly, and \(|p_5|>g_*\), where \(g_*\) is the unique root in \([4,5]\) of the paper's equation (2.15). A tie allows the displayed pivot order. A candidate with \(r<1\) and \(0<|p_5|<g_*\) certifies nonemptiness only. A bounded exclusion must state its exact domain and account for all zero denominators there. A numerical sweep can suggest candidates but cannot certify either a witness or an exclusion.

Verify candidates by comparing successive Schur updates with independent bordered determinants, checking both signed dominance slacks at all entries, and comparing exact isolating intervals with the benchmark. Allow \(r=1\) only for a boundary control. Every sign or tie must be proved by rational arithmetic, polynomial identities, and real-root isolation.

[Paper, equation (2.1) and Conjecture 2.1](https://arxiv.org/html/2602.20390v2) · [Authors' exact computation, pinned revision](https://github.com/jamiecjx/growth-5x5/tree/bd6d13af571d0c454178dac24e753f7cac7a3727)

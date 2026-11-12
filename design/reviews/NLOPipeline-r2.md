# the-enemy — NLO plan (now design/NLOPIPELINE.md), round 2
Attacked commit: 33e790e (plan file, before the record) · Date: 2026-11-12 · Class: plan
Verdict: 0 blocker, 1 major, 3 minor, 0 nit. Report returned inline (plan mode); copied here. Resolutions: NLOPIPELINE.md §5.

## Held
- Code claims of round 1, the fingerprint argument, BMWflow claims (loop_shape, fp, `_fill_p0`, `grid_log0`).
- N3/N4 rederived: O[cR] = E[c]O[R] + O[c]E[R]; the per-pattern memo is valid (s_Q does not depend on ϖ); ϖ → 0 limit
  equals `integratedOf[pr, 1]`.
- N7 structure: the vertex with the external h leg carries p·q, so E[F_{h̄h}] = O(p²) at every ϖ.

## Findings
- E1 [major] A ϖ² breaker cannot make `nloLimit` fail (it vanishes at ϖ → 0); `shiftToZero` had no breaker.
- E2 [minor] Single-family Q terms: the shift w → w − ν must act on the whole numerator polynomial, not on n.
- E3 [minor] The Adimensional LO branch is shared: the NLO changes must be gated, or LO fingerprints change.
- E4 [minor, suspected] No Mathematica-side reference for the BMWflow ϖ = 0 column comparison.

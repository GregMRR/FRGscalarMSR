# the-enemy — design MOMENTUMFUNCTIONS (plan read-claude-md-i-greedy-platypus.md), round 4
Attacked commit: 21719c9 (plan file outside the repo; code unchanged)  ·  Date: 2026-10-23
Verdict: 0 blocker, 2 major (minors/nits not sought this round, by instruction)

(Copy to design/reviews/MOMENTUMFUNCTIONS/r4.md at implementation start. Suite not rerun: code unchanged.)

## Attack surface (re-attack of r3 resolutions, then the rewrite)
- E25 → C1 (r → −r symmetrisation ≡ p → −p at fixed r): held. Any integrand built from scal (and wedge in d = 2) obeys
  L(r,p) = L(−r,−p), so the symmetrised integrand is even in the vector p; its Taylor series has (p²)^a (p·r)^b with b even, hence
  even powers of ap at fixed c as well. Frequency relabelling of the second ordering is absorbed by the complete ω integral (also for
  ∂ν^m at ν = 0). Tadpole (one line) routed through the same substitution integrates to the same value: held.
- E26 → C1 `dsR[c,d][s]` + E2 Oracle ∂sR per line: held. Note (not a finding): `Projection.functionalOrder` is tautological after
  the symmetrisation (the symmetrised form is even whatever the attribution), so the Oracle is the only guard on attribution; E2 has it.
- E27 → C4 (every ŷ = 0 object from today's p → 0 operator): held for η rows, B2 checks, regression; D5 consistency
  β_F̂(0) = η(F̂(0) − 1) holds in the integrated sense. Gap for the functional row's own ŷ = 0 value → **E33**.
- E28 → D7: held for AtDimension = 1 (pointwise TRS identity at c = ±1 re-derived from r3-probe-trs-pointwise; η balance, α̂,
  Rule direction, regulator, misuse aborts all specified). AtDimension ≠ 1 → **E34**.
- E29 → B2: held. Checked against Symmetries.wl:26-40 (omitted Momentum = All) and KPZSetup.wl:14-16: the KPZ tags use explicit
  `"Momentum" -> 0`, catch-all `Neglected[_]` does not match ≥ n₀, Galilean tag on an f^λ slot is order n₀ = 0 explicit: matches.
  Slot count (A4) uses unordered pairs (Symmetries.wl:71), so Fnu in (h,hbar) is one slot: no false `functionSlot` abort.
- E30 → D6: held. grep: no existing test/benchmark uses "Equation" without "Mode" -> "Running" (only
  ShallowWaterDE4Adimensional.wl:15, with Running), so the behaviour change breaks nothing.
- D3 sign of 2ŷ∂ŷF̂ (∂s = k∂k, ∂s ŷ = −2ŷ at fixed p): held.
- D4 K_d = S_{d−1}/(2π)^d with the project's S_d = 2π^{d/2}/Γ(d/2) (CONVENTIONS.md:46-47): S_{d−1}∫(1−c²)^{(d−3)/2}dc = S_d, held.
- C2 subtraction with un-averaged Taylor coefficients (Declared, Protected, Neglected lower orders): held.
- C3 combined residues: the ω integral is analytic at coincident pole positions (also accidental ones at s+ ≠ s− from non-monotone
  F̂), so a fully cancelled rational form has no (X+ − X−) factor: held provided "cancelled before storing" is complete.
- B3 AllMomenta Protected checked pointwise: a pointwise identity in q-routing at generic p implies it in symmetrised midpoint: held.

## Findings

### E33 [major] The functional row is 0/0 at ŷ = 0 and loses all precision as ŷ → 0; the plan fixes no stored ŷ = 0 value and no small-ŷ form
- Where: plan C2 ("divide by p^{2n₀}"), C4 (list of ŷ = 0 objects omits the functional row's own ŷ = 0 value), D2/D3 records
  (only `"Argument"`; no ŷ = 0 entry, no loop variable/measure for it).
- Evidence:
  - For n₀ ≥ 1 the stored row is (I − Σ_{n<n₀} T_n)/ap^{2n₀}. The numerator contains F̂(ar² + ŷ/4 ± ar ap c) − F̂(ar²)-type
    differences and is not symbolically divisible by ap^{2n₀}, so at ap = 0 the expression is Indeterminate. A downstream grid in ŷ
    normally starts at ŷ = 0 (and F̂(0) is the normalised point of D5).
  - Toy n₀ = 2 (B(s) = 1/(1+s), symmetrised midpoint, r = 0.7, c = 0.3, the exact T₀, T₁ subtracted), in double precision vs exact rationals:

    | ap | double | exact |
    |---|---|---|
    | 1e-2 | 0.0124480868 | 0.0124480876 |
    | 1e-3 | 0.0124544659 | 0.0124482000 |
    | 1e-4 | −0.4088 | 0.0124482011 |

    The loss goes as ε/ŷ^{n₀}. For KPZ f_ν (n₀ = 1, T₀ ≡ 0) it is still ε/ap, from the cancellation of the two O(p) halves (r3 E25 table).
- Consequence: the implementer has to decide what the record holds at ŷ = 0. The choices are the pointwise limit of the midpoint
  row (needs one more shape derivative, r3 E27) or the p → 0 operator row at order n₀ (different loop variable and measure, same
  integral). Without a decision, a naive export gives NaN at F̂(0), or garbage at small ŷ, to the downstream integrator.
- Hint: store per functional row an `"AtZero"` piece (the order-n₀ p → 0 operator row, with its own aq and LoopMeasure) and either
  a small-ŷ Taylor form or a documented ŷ_min, and say so in D2/D4 and ShowBeta.
- Resolution: _(author)_

### E34 [major] `Adimensional.related` checks only c = ±1, which is the whole measure only at d = 1. With AtDimension ≠ 1 a false relation passes silently
- Where: plan D7 ("Check `Adimensional.related`: β_FD̂ = α̂ β_Fnû pointwise at c = ±1"). The hint table (`relatedDimension`:
  "Theory d ≠ AtDimension") implies AtDimension may take other values, and "Dim → 1 in every row" contradicts that.
- Evidence:
  - Generic-p integrands are Dim-free pointwise (Dim appears only in Projection.wl:98 for the p → 0 angular average, and in the
    measure). At c = ±1 every scal reduces to the 1-d product.
  - So the c = ±1 check gives the same verdict in every d. For KPZ, TRS holds pointwise at collinear configurations
    (r3-probe-trs-pointwise: ~1e-25), so `Related[FD -> Fnu, "AtDimension" -> 2, …]` passes the check and imposes FD̂ = Fnû,
    η_D = η_ν at d = 2, where TRS is false.
  - The F̂ linear parts carry no explicit d here, so BetaL does not rescue the check.
- Consequence: a nonsensical ansatz goes through with no message. The implementer must also decide what AtDimension ≠ 1, or an
  omitted AtDimension (an FDT-type relation valid in every d), means.
- Hint: either restrict AtDimension to 1 (abort otherwise), or check at generic c (several random c ∈ (−1,1)) whenever d ≠ 1 and
  substitute Dim → AtDimension rather than → 1.
- Resolution: _(author)_

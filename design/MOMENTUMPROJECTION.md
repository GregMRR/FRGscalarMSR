# Momentum projection at any order (p^{2n}) — design

**Status: implemented** (single commit after ba4a80d). Code: `Kernel/Projection.wl` (`projCAverage`, `directionalAverage`,
`explicitProjection`, `checkExplicitVectors`, `checkNormalisation`), Core step 0, Ansatz, Symmetries, Oracle; tests
`tests/unit/ProjectionTest.wl`, `tests/acceptance/ProjectionRegression.wl`, `tests/reference/OldProjection.wl`.

## For the implementing agent (read first)
- Repo: `/home/greg/Physics/Sources/FRGscalarMSR` (Mathematica package; entry `FRGscalarMSR.wl`, modules in `Kernel/`,
  suite `math -script tests/runTests.wl`). Read its `CLAUDE.md`, `docs/CONVENTIONS.md`, `docs/MODULES.md` (Projection,
  Oracle sections) and `design/LEVICIVITA.md` (format to copy, and its F2/F4 wedge facts) before editing.
- Background reports with the evidence and rerunnable probes: `design/reviews/MomentumProjection-r1.md` and `-r2.md`
  (the-enemy rounds 1–2). Finding labels E1–E22 below refer to them.
- All design decisions are settled (user-approved in a grilling session, 2026-10-21); do not reopen them. Decided: general n,
  tested to p⁴ (+ one n = 3 operator test); one projection method for all n; no end-to-end p⁴ acceptance test; keep the
  one-`jac` limit; generic shapes only in the library; **one single final commit**.
- The top priority is the Invariant below: compare against the old operator before trusting anything new.
- Suggested order: step 0 → 1 (+ ProjectionTest, regression vs old operator) → 2 → 3 → 5 → 4, 6, 6b → docs → full suite → commit.
- Before committing: the tree had unrelated uncommitted work (Kernel/Adimensional.wl, Dimensions.wl, JuliaExport.wl, tests,
  docs/handoff/); ask the user how to handle it.

## Context
`ProjectOnMonomials` projects (iν)^m (p²)^n with the script's "Laplacian trick" (`Kernel/Projection.wl:45-57`: bare-symbol
`D[e, p, p]` + `dimFactors` + `lam`, /2d). It is exact only for n = 1 (iterating it gives 8(d²+2) instead of 8d(d+2) at
n = 2), so n ≥ 2 is blocked in five places. Goal: one projection valid for every n, with checks that keep it under control.
Design only for now; implementation later, one final commit. Reviewed by the-enemy, round 1 (1 blocker, 3 major, 6 minor,
2 nit), and round 2 (0 blocker, 4 major); reports in
`design/reviews/MomentumProjection-r1.md`, `-r2.md`. All findings (E1–E22) are addressed below.

## Invariant (top priority)
Every existing pipeline (n ≤ 1: KPZ, LPA, DE4, shallow water, Jacobian2D, Adimensional, JuliaExport round trip) must produce
exactly the same results as before: same parsed monomials, classification, integrands (after Expand), flows, equations,
normalisations and exports (incl. Julia names `Df`, `D2f`). Any difference is a bug in the generalisation, not an accepted
change. (Enemy probe P2: new operator ≡ old on all 16 KPZ monomials, Expand-identical.)

## Facts (labelled M1–M5 to avoid clashing with LEVICIVITA F1–F5)
- **M1 (operator).** For each fixed q, the (p²)^n coefficient = average over directions e of the s^{2n} Taylor coefficient of
  f(p = s e) (equals M2 pointwise because (∇²_p)^n commutes with joint rotations). Substitute scal[p,p] → s²,
  scal[p,q] → s c, wedge[p,q] → s ω; take D[·,{s,2n}]/(2n)! at s = 0; then ⟨c^{2k}⟩ = q^{2k}(2k−1)!!/[d(d+2)⋯(d+2k−2)].
  Frequency part unchanged.
- **M2 (equivalence).** Same as (∇²_p)^n at p = 0 / N_n, N_n = 4^n n! (d/2)_n (N_1 = 2d, N_2 = 8d(d+2)). Verified d = 1..5.
- **M3 (ε, d = 2).** After the LEVICIVITA F2 reduction, odd-ω terms average to 0 (reflection of e about q): the F4 drop is
  exact for every n.
- **M4 (shapes).** Δ depends on (p+q)²: p^{2n} brings Δ/shape derivatives up to order 2n (deltaCompletion, Assemble,
  FrequencyIntegral already general). Requires C^{2n} shapes (E8).
- **M5 (parser).** With ≥ 2 summed indices, `canonical` (Ansatz.wl:89-91) is not a canonical form (lap[lap[A]] → {S1,S1,S2,S2},
  div[grad[lap[A]]] → {S2,S2,S1,S1}). Cosmetic for vertices (sums) but duplicates the monomial table.

## Changes
0. **Sound exact zero test (E1, blocker).** `ExactlyZeroQ` (Core.wl:85-86) maps `GammaAtom[..]` / `DeltaCoefficient[..]` to
   numbers *inside* `Derivative[k][…]`, which then evaluates to 0: derivative terms are invisible. Fix `atomRule` to replace each
   whole atom — `f[..][x]`, `Derivative[k][f[..]][x]`, `Derivative[k][f][x]` — by its own random rational (one per distinct
   atom). Add a unit test in CoreTest: `ExactlyZeroQ[Derivative[1][GammaAtom[hb,h,0]][QQ]]` is False. Existing checks become
   stronger; results unchanged (invariant). Any existing check that now fails is reported, not masked.
1. **Projection (`Kernel/Projection.wl`).** Replace `dimFactors`/`lam` and the `n == 1` branches with the M1 operator for all n;
   no separate `/(2 dim)`. Private helper `directionalAverage`.
   - `Projection.evenInQ` → `Projection.noOddC`: assert no odd power of c at order s^{2n}; labelled an internal invariant
     (unfalsifiable by degree counting: c-degree is 2n − 2a, r2 E18) — never discard odd terms silently (E6).
   - `Projection.complete`: also free of s, c, ω.
   - New `Projection.explicitVectors` (runtime, **n ≥ 2 only**: n ≤ 1 is certified by the regression vs the old operator, and
     this keeps existing builds' timing and behaviour unchanged, r2 E13/E22): explicit components, (∇²_p)^n at p = 0 / N_n,
     compared with the sound zero test (step 0). Dimensions: the declared d when the theory fixes it (also without ε,
     r2 E13), otherwise Dim = 1, …, n+1 (n+1 independent d-structures, E4). Cost ~0.03–0.3 s per KPZ p⁴ integrand.
   - `Projection.wedgeOdd` (Projection.wl:38): explicit (∇²)^n for all n, not only n == 1 (E5).
   - `Projection.LHS` stays but no longer tests the d-normalisation (E2): add `Projection.normalisation`, a runtime probe that
     applies the full `project` to synthetic integrands — (p·q)^{2k}(p²)^{n−k}, k = 0..n — against a reference written
     independently of `directionalAverage`: ⟨c^{2k}⟩ = q^{2k} Γ(k+½)Γ(d/2)/(√π Γ(k+d/2)) (r2 E17), for each n in the truncation.
2. **Ansatz (`Kernel/Ansatz.wl`).** Lift `Ansatz.momentumOrder` (:196). Fix powers (r2 E14): `S[Power[x_, n_]]` (:68)
   reuses one dummy for every factor, so `lap[h]^2`, `(grad h.grad h)^2` abort at `Ansatz.indexPairs`; use
   `Expand[Times @@ Table[S[x], {n}]]` with fresh dummies per factor (index-free powers give the same result: invariant).
   `canonical` (E9): for each permutation of the dummy
   labels — relabel, sort each factor's index list (derivatives commute; EpsIndex kept as it is), sort the factor list —
   and keep the minimum; ε: also over the slot swap, sign −1 for a swap, zero iff the minimum over swapped relabellings
   equals the minimum over unswapped ones. Update usage strings (`FieldFactor`, `SpaceIndex`) and the `singleEpsilon` hint
   (:178, drop "p^4 is unsupported"; keep one-jac limit).
3. **Symmetries (`Kernel/Symmetries.wl`).** Remove `Symmetries.momentumOrder` cap (:60); `declOutside` scan (:94) up to the
   ansatz's max MomentumOrder/2.
4. **Regulator smoothness (E8, r2 E16).** The library only accepts generic declared shapes (`Regulator.onlyDeclaredShapes`,
   Regulator.wl:74-77), so no explicit shape can reach the projection; generic shapes are fine at any order. Explicit shapes
   are substituted only after the library, in test helpers (KPZ.wl:13-16, KPZAdimensional.wl:13: Litim as `1/y − 1` on y < 1).
   - Projection card: ASSUMPTION line when some n ≥ 2: "p^{2n} flows contain shape derivatives up to order 2n: any explicit
     shape substituted later must be C^{2n}".
   - Abort where shapes become explicit: the KPZ Litim substitution helper aborts if the expression it substitutes into contains
     any `Derivative[k][shape]` with k ≥ 3 (p² flows stop at r″), with
     "an explicit non-smooth (Litim) shape cannot be used with p^4 or higher flows". Existing KPZ tests unaffected (they stop at p²).
   - expectFail test for that abort.
5. **Oracle (`Kernel/Oracle.wl:76-90`).** Nested central stencils over direction tuples (c₁…cₙ), divide by N_n (replaces
   `(2 d)^n`). For n ≥ 2: PrecisionGoal ≈ 8 + 10n + margin, WorkingPrecision raised, and Richardson extrapolation (h, h/2)
   to restore headroom against the 10⁻⁸ tolerance (r2 E19). The memo key (Oracle.wl:76) gains the precision so n ≥ 2
   stencils never reuse n ≤ 1 points (r2 E15). n ≤ 1 records keep today's precision, stencil and values (invariant). Oracle card
   states that at n ≥ 2 the d-dependence is certified by `explicitVectors` (the Oracle runs only d = 2, 3).
6. **JuliaExport (`Kernel/JuliaExport.wl:22`).** Keep `D`, `D2` and the abort for k ≥ 3 unchanged in behaviour; only its message
   now says p⁴ flows need r‴/r⁗ and a FlowRG interface extension (E7, r2 E20).
6b. **Adimensional.wl:64** keeps its own atom pattern (harmless: `AdimAtom` has no argument, shapes are symbols); add a comment
   noting the restriction, no code change (r2 E21, invariant).
7. **Closure/Assemble.** No code change (project already returns the coefficient); fix the doc text MODULES.md:460.

## Tests
- Unit `tests/unit/ProjectionTest.wl` (new):
  - p⁴: ⟨f((p+q)²)⟩ coefficient = f″/2 + 2q² f‴/d + 2q⁴ f⁗/[d(d+2)];
  - p⁶: f‴/6 + q² f⁗/d + 2q⁴ f⁽⁵⁾/[d(d+2)] + 4q⁶ f⁽⁶⁾/[3d(d+2)(d+4)];
  - N_n vs explicit (∇²)^n (p²)^n, n = 1..3, d = 1..4; wedge-odd term averages to 0 in d = 2;
  - mutation test: a wrong ⟨c⁴⟩ (e.g. 3/(7d−6)) injected into `directionalAverage` only (not the Γ-function reference) must
    make `explicitVectors` and `normalisation` fail (guards E1/E4/r2 E17);
  - expectFail tests (house convention, r2 E18) from synthetic inputs for `Projection.explicitVectors`,
    `Projection.normalisation`, `Projection.wedgeOdd` at n = 2, and the Litim-at-p⁴ abort.
- Regression (invariant): old Laplacian `project` kept test-only; new ≡ old (Expand-identical) for every integrand of every
  acceptance setup; full existing suite green with unchanged outputs.
- CoreTest: ExactlyZeroQ sees derivative atoms (step 0).
- AnsatzTest: `lap[h]^2`, `lap[h] lap[h]`, `(grad h.grad h)^2` parse (r2 E14); `lap[lap[A]]` ≡ `div[grad[lap[A]]]`; `(grad A.grad A)^2` with swapped dummies → one monomial; ε sign/zero rule
  with 2 dummies; p⁴ monomials accepted.
- OracleTest: N_n stencil normalisation on a known polynomial; precision adequate at n = 2.
- No end-to-end p⁴ acceptance test (decided).

## Docs
This file (`design/MOMENTUMPROJECTION.md`, facts M1–M5, steps, checks) and `design/reviews/MomentumProjection-r1.md`,
`-r2.md` (already in place). Round 2: 0 blocker, 4 major, 3 minor, 3 nit — all addressed
above; r2 confirmed step 0 breaks no existing check and new ≡ old on KPZ 16/16, jac-2D 36/36, shallow-water DE4 56/56.
Update CONVENTIONS (Projections), MODULES (Core ExactlyZeroQ, Ansatz, Symmetries, Projection incl. :386 LHS text, Closure :460,
Oracle, JuliaExport), FUTURE (remove p⁴ item; add FlowRG r‴/r⁗, lift one-jac limit, shallow-water p⁴ acceptance),
FRGDESIGN E2 (correct "24d²" → 8(d²+2)).

## Delivery
Single final commit after the full suite passes. Pending uncommitted work in the tree must be handled first.

## Verification
`math -script tests/runTests.wl`: whole suite green with unchanged outputs (invariant), new tests green, mutation test
red→caught. Manually run a small KPZ + κ lap[lap[h]] ansatz at `"Momentum" -> 2` to see LHS, normalisation, explicitVectors,
Oracle d = 2, 3 pass and to gauge runtime.

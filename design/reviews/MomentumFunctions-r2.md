# the-enemy — design MOMENTUMFUNCTIONS (plan read-claude-md-i-greedy-platypus.md), round 2
Attacked commit: 21719c9 (plan file outside the repo; code unchanged)  ·  Date: 2026-10-23
Verdict: 1 blocker, 2 major, 6 minor, 1 nit

(Plan mode allows writing only this file. Copy it to design/reviews/MOMENTUMFUNCTIONS/r2.md, and save the appendix probes as
r2-probe-*.wl. The suite was not rerun: it was green in round 1 (486/486, 2 skipped) and the code has not changed.)

## Attack surface
- Round-1 resolutions:
  - E1 (D10, plan:40-41) — rederivation (r2-probe-subtraction) — held for n₀ ≤ 1, **fails for n₀ ≥ 2** → E15
  - E2 (D8, plan:34-37) — held: lower orders are classified, and the KPZ shift tag keeps its meaning
  - E3 (FrequencyIntegration, plan:77-79; D14) — held
  - E4 (D6, plan:27-29) — held for (a)–(d); the level at which "one term" is judged is unspecified → E21
  - E5 (D9, plan:38-39) — held
  - E6 (D5, plan:21-26) — held in substance; 4 of the 5 rejections have no check name → E22
  - E7 (D12, plan:45-53) — Related → E16 (blocker); FixedValue / AllP / "Dimension" / where tags go → E18
  - E8 (functionBare) — held
  - E9 (D7 aborts, plan:32-33) — held (the check names are missing → E22)
  - E10 (D11 measure) — held: K_d is defined, d = 1 is a separate case, the measure is stored per row
  - E11 (cascade, plan:66-87) — held
  - E12 (D13) — held. Side effect: the LO approximation on sibling slots is silent unless the user opts in → E20
  - E13 (D15) — held as a commitment. The arXiv IDs, equation numbers and errata are still to be written; I could not verify them offline.
  - E14 (notation) — held, but D11 uses both ap and ŷ for the same argument → E23
- D5 DSL on KPZ LO, on a mixed ansatz with a mass term, on a p⁴ function over a constant p², and on f^λ — written out below — findings E15, E17, E19, E20
- D6 vs Dimensions (probe r2-probe-flam) — E17
- D7 Normalise vs the existing amplitude idiom (`nk ank`, "Mode", "Equation", "Dimensionless"; MODULES:674-707) — E19
- D8 + "TruncationOrders" for a mixed ansatz — E20
- D12 Related (probe r2-probe-related) — E16. AllP / FixedValue / "Dimension" / which function takes the tags (Symmetries.wl:56-58, 33-34) — E18
- Abort hints for every new check, and the existing hints that the new atoms reach — E22, E17
- Outputs: ShowMonomials, ShowClassification, ShowEquation/EquationOf, ShowBeta, atoms, "State" — E23

## Findings

### E15 [major] (resolution of E1 fails for n₀ ≥ 2) After subtracting the direction-averaged lower orders, the symmetrised integrand is not O(p^{2n₀}) pointwise
- Where: plan:34-36 (D8: "F's flow = [flow(p) − Σ_{n<n₀} flow_n p^{2n}]/p^{2n₀}"), plan:40-41 (D10: "the symmetrised, subtracted integrand is O(p^{2n₀}) pointwise"), plan:74-75 (Projection: "subtract lower orders; symmetrise; divide"). Today's p^{2n} coefficient is direction-averaged (MODULES:378-381, `projCAverage`).
- Evidence (r2-probe-subtraction.wl): generic rotation-invariant integrand I = f(q², |p+q|²) + g(q²)(p·q)², with the lower orders taken from today's operator and the subtracted piece symmetrised as well (the best case).
  - n₀ = 1: the s¹ coefficient is 0. Held: subtracting I(0,q) and then symmetrising works, mass term included.
  - n₀ = 2: the s¹ coefficient is 0, but the s² coefficient is `(−2413161750121283 + 83 d (10400200000000 + 9954083736401 x²))/(9928281844741383 d)`, which is 0.0891 at d = 3, x = 1 and 0.0059 at x = 0. The s³ coefficient is −0.0143.
  - The s² remainder is nonzero even after the angular average at fixed aq: it cancels only over the whole q-integral, because the shift q → −p−q moves weight between shells.
  - Reason: the un-averaged s² coefficient depends on c, while today's operator replaces c² by 1/d.
- Consequence: D8 explicitly allows a p⁴ function over a constant p² (`-nk hbar lap[h] + tk Ft hbar lap[lap[h]]`, the functional analogue of KS). There `Projection.functionalOrder` aborts on a valid ansatz. If the check were relaxed, F̂ would carry 1/ap² and 1/ap pieces. The KPZ_LO benchmark (n₀ ∈ {0, 1}) cannot see this.
- Hint: route the loop at the midpoint, r = q + p/2. There the symmetrised integrand is even in p at fixed r, so subtract its *un-averaged* Taylor coefficients pointwise; they integrate to the same flow_n. Alternatively, restrict D8 to n₀ ≤ 1 with a named abort.
- Resolution: _(author)_

### E16 [blocker] The Related example is wrong for α ≠ 1, its "When" rules are dimensionful, and nothing requires the output to obey them
- Where: plan:49-50 (`Related[FD -> α Fnu, "Reason" -> …, "When" -> {nk -> α Dk, rn -> α rD}]`, "identity checked on the dimensionless β's (β_FD̂ = α β_Fnû) … Fnu stays in the state").
- Evidence (r2-probe-related.wl): time reversal written exactly as in the KS benchmark (Checks.wl:20, 32: ν = αD, τ = αζ, r_ν = αr_D), on Γ2 + R: Re(Γ11 + R11)/p² = T·(−(Γ02 + R02)/2) for every p.
  - With the plan's three substitutions there is no T for α ≠ 1: KS-style regulator `{al Dk − Dk T, al Dk − al Dk T}` → T: {}; KPZ-style regulator (`nk q² rn`) `{al² Dk − Dk T, …}` → T: {}. Only α = 1 works.
  - The consistent rules are FD → Fnu (no α), nk → α Dk, and rn → α rD (KS regulator) or rn → rD (KPZ regulator). Both give T = α.
  - The correct rn rule therefore depends on how the user wrote the regulator. The library knows it (Regulator stage); the user is asked to retype it.
  - In KPZ, nk and Dk are independent references: `Dimensions.referencesIndependent` passes for {nuk, Dk} (tests/acceptance/KPZAdimensional.wl:9). So α = nk/Dk is dimensionful, and `nk -> α Dk` is not homogeneous. It also has nothing to act on, since the dimensionless β's no longer contain references. Its dimensionless content is η_ν = η_D, which neither the plan nor the "When" list states. With the default Fixed references (D7, MustFix) the check cannot hold unless η_ν = η_D is imposed.
  - "FD replaced by α Fnu, Fnu stays in the state": the plan does not say that the output is evaluated under the When conditions (rn = rd, η_ν = η_D, d = 1). Dropping FD's equation is valid only on that submanifold.
- Consequence: as specified, the check aborts on correct physics (d = 1 KPZ, f_D̂ = f_ν̂) unless α = 1. If it is "fixed" by tuning α, the reduced system is exported outside the conditions that justify it.
- Hint: state the relation once at the Γ2 level and let the library derive the rest (F map, regulator relation, η equality):
  - before: `Related[FD -> α Fnu, "Reason" -> "TRS", "When" -> {nk -> α Dk, rn -> α rD}]` (plus a separate `"Dimension" -> 1`)
  - after:  `Related[nk Fnu -> α Dk FD, "Regulator" -> True, "Dimension" -> 1, "Reason" -> "time reversal (d = 1)"]`
  - The library derives FD̂ = Fnû, η_D = η_ν, and the R relation from the Regulator stage (aborting with a named check if R violates it); it imposes them on the output; and it prints them as ASSUMPTION lines.
- Resolution: _(author)_

### E17 [major] A function with no reference prefactor (f^λ in the iω slot, the natural next KPZ step) leads to a dead end: the existing hint says to do what D5 forbids
- Where: plan:28-29 (D6 "otherwise F is dimensionful with derived dimension"), plan:24 (D5 "Dimensionless" rejects functions), plan:57-58 (D15 "f^λ ≡ 1").
- Evidence (r2-probe-flam.wl, today's code, with a constant `Flam` standing in for the function):
  `hbar (Flam dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2` →
  `FAIL Dimensions.freedoms … "undetermined" -> {"time", h, hbar, Flam, lamk} … hint: declare the amplitudes in "Dimensionless" (e.g. Ak in nk Ak)`.
  Under D5 the user follows that hint and then hits the "Dimensionless rejects functions" abort. The user has no way to say that F^λ is dimensionless, although D6's rule on dimensions has to be able to express this.
- Consequence: an ansatz that the D-rules allow (one function, n₀ = 0, m = 1, with D2 and D6 both satisfied) cannot get dimensions, and the hints send the user in a loop.
- Hint: allow functions in `"Dimensionless"`, or let D6 infer "a function whose constant counterpart would be numeric is dimensionless". Teach `Dimensions.freedoms` to name functions in its hint:
  - before: `"MomentumFunctions" -> {Flam}` … `DeriveDimensions[an, "References" -> {nk -> "nu", Dk -> "D"}]` → freedoms abort → `"Dimensionless" -> {Flam}` → abort
  - after:  the same call with `"Dimensionless" -> {Flam}` accepted (or not needed). The freedoms hint ends with "…or, for a MomentumFunction, declare it in "Dimensionless"".
- Resolution: _(author)_

### E18 [minor] Tag DSL: AllP collides with All, two spellings appear, FixedValue duplicates Protected, and the plan does not say which function takes Related or FixedValue
- Where: plan:46-53. Symmetries.wl:33-34 (`Symmetries.tagOptions` allows only Frequency/Momentum/Reason), Symmetries.wl:56-58 (`Symmetries.arguments` allows only Protected/Neglected/NonRenormalised), Symmetries.wl:26-29 (`All` = every order up to the truncation).
- Evidence:
  - (a) `"Momentum" -> All` (today: each p^{2n}, n ≤ mp, at p → 0) and `"Momentum" -> AllP` (generic p) differ in meaning by one letter. AllP is a new bare symbol where the library otherwise uses strings ("Running", "Fixed"). D12 also spells the same tag positionally: `Neglected[…, AllP]` (plan:46-48).
  - (b) `FixedValue[Fnu, 0 -> 1]`. The check is "flow at ŷ = 0 is 0", so the value 1 is never used, and the plan does not say whether this is the dimensionful ∂s F(0) or β_F̂(0). These differ by the canonical term x·F̂(0) when F is dimensionful. Its physics is the Galilean Ward identity, which today is written `Protected[{hbar,h}, "Frequency" -> 1, "Momentum" -> 0]`. The user has to change spelling when λ becomes f^λ.
  - (c) `"Dimension" -> 1` reuses the name of the Theory option, which takes `Dim | FixedDimension[2]`. Nothing says what happens when the Theory declared d = 2, or whether Dim → 1 is substituted too ("checks draw c ∈ {±1} only").
  - (d) Neither Related nor FixedValue is assigned to a function. The Related check needs dimensions (Adimensional), so it cannot go in DeclareSymmetries without being carried through four stages.
- Hint (before → after):
  - `Neglected[{hbar,h}, "Frequency" -> 1, "Momentum" -> AllP, …]` → `Neglected[{hbar,h}, "Frequency" -> 1, "Momentum" -> "AllMomenta", …]` (a string)
  - `FixedValue[Flam, 0 -> 1, "Reason" -> "Galilean"]` → `Protected[{hbar,h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean"]`, which on a functional slot means "the dimensionful F(0) does not flow" (the existing tag, unchanged across ansätze)
  - `Related[…, "Dimension" -> 1]` → `"AtDimension" -> 1` (abort if the Theory fixed another d)
  - State that `Related` is an `AdimensionalEquations` option argument, next to "State".
- Resolution: _(author)_

### E19 [minor] "Normalise" introduces a third way to express "a reference times an amplitude", and the natural migration from `nk ank` aborts twice
- Where: plan:31-33 (D7), plan:24 (D5: "Dimensionless" and "Equation" reject functions); MODULES:674-675, 704-706; tests/acceptance/ShallowWaterDE4Adimensional.wl:14-15 (`nk -> {"nu", "Mode" -> "Running", "Equation" -> nk ank}`, `"Dimensionless" -> {…, ank}`).
- Evidence: a physicist upgrading the constant amplitude `nk ank` to `nk Fnu` keeps the working options. Both lines abort: `"Equation" -> nk Fnu`, and `"Dimensionless" -> {Fnu}` (which D6 already implies). The replacement also repeats itself: Normalise requires `"Mode" -> "Running"`, and Mode Fixed with Normalise aborts. So Mode carries no information.
  - KPZ LO, as the plan has it:
    `"References" -> {nk -> {"nu", "Mode" -> "Running", "Normalise" -> Fnu}, Dk -> {"D", "Mode" -> "Running", "Normalise" -> FD}}`
  - after (Mode inferred; "Equation" accepts the functional combination with the same meaning as for `nk ank`, F̂(0) = 1 ⇔ amplitude counted as 1):
    `"References" -> {nk -> {"nu", "Equation" -> nk Fnu}, Dk -> {"D", "Equation" -> Dk FD}}` (or keep "Normalise" as an alias)
  - Accept `"Dimensionless" -> {Fnu}` when it agrees with the derived dimension; abort only when it contradicts it.
- Hint: one concept, one key. Infer Mode from Normalise/Equation.
- Resolution: _(author)_

### E20 [minor] What `"TruncationOrders" -> <|"Momentum" -> mp|>` means for a mixed ansatz, and which slots stay at p → 0 silently, is left to guesswork
- Where: plan:34-37 ("Truncation orders apply to the lower orders only"), plan:54-55 (D13), Symmetries.wl:73-82, 91-95.
- Evidence: KPZ LO under the design is
  `hbar (dt[h] - nk Fnu lap[h] - lamk/2 grad[h].grad[h]) - Dk FD hbar^2`, `"MomentumFunctions" -> {Fnu, FD}`, with today's `"TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>`.
  - (hbar,hbar) p¹ lies inside FD's range (n₀ = 0). The plan never says whether it becomes class Functional or is caught by the user's catch-all `Neglected[_]`. In the second case it is reported as dropped while FD absorbs it.
  - The user cannot tell whether mp = 0 or 1 is the right choice.
  - The LO approximation itself (f^λ ≡ 1: the (hbar,h) iω¹ slot frozen at all p) is reported only if the user adds a `Neglected[…, AllP]` tag. Without it D13 runs the slot at p → 0 with no WARN, so the main approximation of the scheme is silent by default.
- Hint: print every order ≥ n₀ of a functional slot as one row of class `Functional` (and do not let tags match there). Ignore mp on functional slots (card line "mp applies to function-free slots"). For every (pair, m) sibling of a functional entry that has no AllP tag, emit ASSUMPTION `Symmetries.siblingAtZeroMomentum`:
  - before (silent): no tag → nothing on the card
  - after: `ASSUMPTION Symmetries.siblingAtZeroMomentum: {hbar,h} (I w)^1 computed at p -> 0 only (add Neglected[…, "Momentum" -> "AllMomenta"] to see its size)`
- Resolution: _(author)_

### E21 [minor] D6's "one term" does not say whether it is judged on the monomial table or on the Γ2 slot
- Where: plan:27-28 (`Ansatz.functionalTerm` is listed under Ansatz, plan:64-65).
- Evidence: `-nk Fnu hbar lap[h] - nk Fnu grad[hbar].grad[h]` gives two monomials (the canonical form does not integrate by parts, MODULES:117-120) and one Γ2 term 2 nk Fnu p². At the monomial level it aborts, although the ansatz is meaningful. At the Γ2 level it passes. A sum whose pieces partly cancel by integration by parts is judged differently by the two levels.
- Hint: check on the Γ2 slot coefficient (after Vertices), and quote the slot in the FAIL block.
- Resolution: _(author)_

### E22 [minor] Four rejections have no check name, and no hint texts are specified; two existing hints mislead on the new inputs
- Where: plan:23-25 ("every option that takes couplings rejects them with a named check (… `NonRenormalised`, `"References"`, `"Dimensionless"`, `"Equation"`)": only `Regulator.noFunctions` is named). CLAUDE.md "New check: … named `Module.check` + an expectFail".
- Evidence: the abort chain of E17. `Symmetries.contradiction`'s hint ("drop the tag or the coupling", Symmetries.wl:78-79) is what a user gets for Protected on a functional slot, where the fix is FixedValue/Related, or the tag at "Momentum" < n₀.
- Hint: name them `Symmetries.nonRenormalisedFunction`, `Dimensions.functionReference`, `Dimensions.functionDimensionless` (if E19 is not adopted), `Adimensional.functionEquation`. Put this check → hint table in the design record:

| check | trigger | hint the user should read |
|---|---|---|
| Ansatz.functionQuadraticOnly | `Fnu hbar h^2` | "MomentumFunctions multiply quadratic terms only (Γ3/Γ4 couplings stay constant, D1)" |
| Ansatz.functionNotInOperator | `lap[Fnu h]` | "write the function outside the operators: `Fnu lap[h]`" |
| Ansatz.functionBare | `Fnu[0]` | "write `Fnu` bare: it already means Fnu(p²)" |
| Ansatz/Symmetries.functionalTerm | `nk Ak Fnu`, `Fnu^2`, `1/Fnu`, F at two powers | "one term ref × Fnu × (p²)^n₀: absorb `Ak` / the power into the function" |
| Symmetries.functionalSlot | constant at ≥ n₀, or two functions | "`-nk hbar lap[h]` is F(0)'s job: drop it" / "Fnu + Gk is one function: keep one" |
| Symmetries.functionSlot (0 slots) | `Fnu div[hbar grad[h]]` | "this term vanishes in Γ2 (total derivative / odd ω on a diagonal entry): drop it" |
| Symmetries.functionSlot (2 slots) | `nk Fnu …, Dk Fnu …` | "one function per slot: use two functions and `Related[…]`" |
| Projection.functionalOrder | E15 / lower orders nonzero | quote the first non-vanishing order p^{2j}, j < n₀, and its slot |
| Regulator.noFunctions | `nk Fnu q^2 rn` | "regulators take constant couplings: use the reference `nk`" |
| Dimensions.freedoms (existing) | E17 | add "…or declare the MomentumFunction in "Dimensionless"" |
| Ansatz.pipeline (existing, adapted) | user `"Pipeline"` contradicts the functions | "drop `"Pipeline"`: it is set from `"MomentumFunctions"`" |
- Resolution: _(author)_

### E23 [minor] The output naming and addressing of functional equations are unspecified
- Where: plan:42-44 (D11), plan:82-86 (Assemble, Adimensional); MODULES:500-546 (record keys `Flowing`, `Momentum`; `EquationOf[eqs, comb]`), MODULES:519 (ScaleAtom[a,b,m,j]), 750 (ShowBeta[adim, c]), 734 ("State").
- Evidence:
  - The dimensionless symbol of F is not named. A reference prefactor makes F dimensionless, so by Dimensions.names it "keeps its name". Fnu[q²] (dimensionful argument) and Fnu[aq²] (F̂) then share one symbol.
  - ScaleAtom/AdimAtom have no momentum argument, yet they must exist at s₁ and at s₂.
  - The plan does not say whether a functional equation is addressed as `EquationOf[eqs, nk Fnu]`, `EquationOf[eqs, Fnu]`, or `ShowBeta[adim, Fnu]`, nor whether functions belong in "State".
  - D11 writes the argument as ap² while BetaL uses ŷ.
- Hint (after):
  - `EquationOf[eqs, Fnu]` and `ShowBeta[adim, Fnu]` accept the function symbol.
  - The record has `"Momentum" -> "Functional"` and `"Argument" -> ap^2`.
  - F̂ is named `aFnu` (or `"Names"`).
  - Atoms are `AdimAtom[a,b,m,j][aq^2]` vs `[ap^2 + aq^2 + 2 ap aq c]`.
  - `ShowMonomials` prints `-nk Fnu(p²)`; `ShowClassification` prints `{hbar,h} | (I w)^0 (q^2)^{≥1} | Functional | nk Fnu(p²)`.
  - Functions are listed in "State" by default, in declaration order.
- Resolution: _(author)_

### E24 [nit] Presentation slips in the examples
- plan:21: `{Fnu, FD, Gk}` — `Gk` has the coupling-style k suffix and is never used. `Hbar lap[H]` is written in capitals, but the KPZ fields are `hbar`, `h`.
- plan:49: `rD`, where the shapes in KPZSetup/KS are `rn`, `rd`.

## Usability walkthrough (the cases asked for)
- **KPZ LO.** `DefineAnsatz[th, "Density" -> hbar (dt[h] - nk Fnu lap[h] - lamk/2 grad[h].grad[h]) - Dk FD hbar^2, "Couplings" -> {nk, lamk, Dk}, "MomentumFunctions" -> {Fnu, FD}]` is obvious and minimal. It reads as today's KPZ with two symbols added; held.
  - The difficulty is downstream: the TruncationOrders/Functional rows and the silent f^λ (E20), the Normalise verbosity (E19), and the time-reversal tag (E16).
- **Mixed with a mass term.** `hbar (dt[h] + mk h - nk Fnu lap[h]) - Dk FD hbar^2`. mk is a lower order of Fnu's slot (D8), and the n₀ = 1 subtraction holds (probe). Clear.
- **A p⁴ function over a constant p²** (`-nk hbar lap[h] + tk Ft hbar lap[lap[h]]`). Easy to write, but it aborts on `functionalOrder` (E15).
- **Time reversal.** It cannot be written correctly with the D12 syntax for α ≠ 1, and it needs the regulator writing retyped (E16).
- **f^λ.** Easy to write, but it reaches the dimension dead end (E17).
- **Consistency with the existing library.** "MomentumFunctions" as a DefineAnsatz option, function symbols written bare, and the reuse of "References"/"Mode" all match the library's style. Where they do not match: AllP is a bare symbol (E18), Normalise duplicates Equation (E19), and FixedValue duplicates Protected (E18).

## Nonsensical / ambiguous ansatz cases (round 2)
| case | named check in the design? |
|---|---|
| function in a cubic or quartic monomial | yes: Ansatz.functionQuadraticOnly |
| function under an operator | yes: Ansatz.functionNotInOperator |
| undeclared function / function also listed in "Couplings" | yes: knownSymbols / fresh-symbol checks (D5) |
| applied symbol `Fnu[0]` | yes: Ansatz.functionBare |
| `Fnu^2`, `1/Fnu`, F at two powers, `nk Ak Fnu` | yes: Ansatz.functionalTerm (level unspecified, E21) |
| constant at power ≥ n₀ in a functional slot; two functions in one slot | yes: Symmetries.functionalSlot |
| one function in two slots / in zero slots | yes: Symmetries.functionSlot |
| lower p-orders of a functional slot that do not vanish and are not declared | yes: classification as today (D8) + Projection.functionalOrder |
| function in the regulator | yes: Regulator.noFunctions |
| NonRenormalised[Fnu], function in References/Dimensionless/Equation | rejected, but no check name (E22) |
| Normalise misuse (wrong reference, dimensionful F, Mode Fixed, with FixedValue) | rejected (D7), no check names (E22) |
| user "Pipeline" contradicting the functions | yes (D5, check name unstated) |
| Protected/Neglected AllP on a functional slot | yes: Symmetries.contradiction (hint misleads, E22) |
| catch-all `Neglected[_]` matching orders ≥ n₀ of a functional slot | unspecified (E20) |
| Related with an inconsistent relation / a regulator that breaks it / outside its conditions | none: the example itself is inconsistent; the regulator and η conditions are not derived or checked (E16) |
| Related/FixedValue given to the wrong function, `"Dimension" -> 1` with a Theory at d = 2 | unspecified (E18) |
| valid p⁴ function over a constant p² | **wrongly rejected** by functionalOrder (E15) |
| valid f^λ (function without a reference) | **wrongly rejected**, with a looping hint (E17) |
| `(nk + Dk) Fnu` | yes: functionalTerm (two terms) |
| function absent from the density | unspecified (WARN as for `Ansatz.unusedCouplings`?) |
| quadratic `jac` with a function | yes: Ansatz.jacQuadratic (existing) |

## Appendix: probes (`math -script <file>`)

### r2-probe-subtraction.wl
```mathematica
(* is [flow(p) - sum_{n<n0} flow_n p^2n]/p^(2 n0) O(1) pointwise after q <-> -p-q symmetrisation? generic I, d symbolic;
   lower orders from today's operator (n=0: p=0 value; n=1: direction-averaged s^2 coefficient), subtracted piece symmetrised *)
f[a_, b_] := 1/((1 + a) (2 + b)^2); g[a_] := 1/(3 + a);
qq = 7/10;
P2[s_, x_] := s^2 + qq^2 + 2 s qq x;
I1[s_, x_] := f[qq^2, P2[s, x]] + g[qq^2] (s qq x)^2;
I2[s_, x_] := f[P2[s, x], qq^2] + g[P2[s, x]] (-s^2 - s qq x)^2;   (* q -> -p-q *)
sym[s_, x_] := (I1[s, x] + I2[s, x])/2;
c0F[Q2_] := f[Q2, Q2];
c2F[Q2_] := Module[{e = f[Q2, sg^2 + Q2 + 2 sg Sqrt[Q2] xx] + g[Q2] (sg Sqrt[Q2] xx)^2},
   Expand[SeriesCoefficient[e, {sg, 0, 2}]] /. xx^2 -> 1/d /. xx -> 0];
Lfun[Q_, s_] := c0F[Q^2] + s^2 c2F[Q^2];
Lsym[s_, x_] := (Lfun[qq, s] + Lfun[Sqrt[P2[s, x]], s])/2;
rem1 = Series[sym[s, x] - (Lfun[qq, 0] + Lfun[Sqrt[P2[s, x]], 0])/2, {s, 0, 2}];
Print["n0=1  coefficient of s^1 (must be 0): ", Simplify[SeriesCoefficient[rem1, 1]]];
rem2 = Series[sym[s, x] - Lsym[s, x], {s, 0, 3}];
Print["n0=2  s^1: ", Simplify[SeriesCoefficient[rem2, 1]]];
Print["n0=2  s^2 (must be 0 for O(p^4)): ", Simplify[SeriesCoefficient[rem2, 2]]];
Print["n0=2  s^2 at d=3, x=1 and x=0: ", N[SeriesCoefficient[rem2, 2] /. d -> 3 /. x -> 1], "  ", N[SeriesCoefficient[rem2, 2] /. d -> 3 /. x -> 0]];
Print["n0=2  s^3: ", N[SeriesCoefficient[rem2, 3] /. d -> 3 /. x -> 1/2]];
Print["n0=2  angular average at fixed aq of the s^2 remainder: ",
  Simplify[Expand[Expand[Normal[SeriesCoefficient[rem2, 2]]] /. x^4 -> 3/(d (d + 2)) /. x^3 -> 0 /. x^2 -> 1/d /. x -> 0]]];
```
Output: n0=1 s¹ 0; n0=2 s¹ 0; s² = (−2413161750121283 + 83 d (10400200000000 + 9954083736401 x²))/(9928281844741383 d); 0.0891, 0.0059; s³ −0.0143; angular average 200000000(−22736 + 12367 d)/(28447798982067 d) ≠ 0.

### r2-probe-related.wl
```mathematica
cond[r11_, r02_] := (nk Fnu[y] + r11) - T (Dk FD[y] - r02/2);   (* Re(G11+R11)/p^2 - T (-(G02+R02)/2) *)
plan = {FD -> (al Fnu[#] &), nk -> al Dk, rn -> (al rD[#] &)};
ks = cond[Dk rn[y], -2 Dk rD[y]] /. plan;        (* KS regulator, benchmarks/KuramotoSivashinsky/Pipeline.wl:21 *)
kpz = cond[nk rn[y], -2 Dk rD[y]] /. plan;       (* KPZ regulator, tests/reference/KPZSetup.wl *)
solveT[e_] := Solve[Thread[DeleteCases[Flatten[CoefficientList[Expand[e /. {Fnu[y] -> u, rD[y] -> v}], {u, v}]], 0] == 0], {T}];
Print[solveT[ks], solveT[kpz], solveT[ks /. al -> 1], solveT[kpz /. al -> 1]];
Print[solveT[cond[Dk rn[y], -2 Dk rD[y]] /. {FD -> (Fnu[#] &), nk -> al Dk, rn -> (al rD[#] &)}],
      solveT[cond[nk rn[y], -2 Dk rD[y]] /. {FD -> (Fnu[#] &), nk -> al Dk, rn -> (rD[#] &)}]];
```
Output: `{} {} {{T->1}} {{T->1}}`, then `{{T->al}} {{T->al}}`.

### r2-probe-flam.wl
```mathematica
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
th = Block[{Print}, DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
an = Block[{Print}, DefineAnsatz[th, "Density" -> hbar (Flam dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2,
     "Couplings" -> {Flam, nuk, lamk, Dk}]];
r = CheckAbort[DeriveDimensions[an, "References" -> {nuk -> "nu", Dk -> "D"}], "ABORTED"];
Print["result: ", If[r === "ABORTED", r, "accepted"]];
```
Output: `FAIL check: Dimensions.freedoms … "undetermined" -> {"time", h, hbar, Flam, lamk} … hint: … declare the amplitudes in "Dimensionless" (e.g. Ak in nk Ak) …`, result ABORTED.

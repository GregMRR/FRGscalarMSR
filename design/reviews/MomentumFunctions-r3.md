# the-enemy — design MOMENTUMFUNCTIONS (plan read-claude-md-i-greedy-platypus.md), round 3
Attacked commit: 21719c9 (plan file outside the repo; code unchanged)  ·  Date: 2026-10-23
Verdict: 1 blocker, 3 major, 4 minor, 0 nit

(Plan mode allows writing only this file. Copy it to design/reviews/MOMENTUMFUNCTIONS/r3.md, and the probes below to
r3-probe-*.wl; they currently sit in the session scratchpad
/tmp/claude-1000/-home-greg-Physics-Sources-FRGscalarMSR/c908efa9-dd0e-454d-bfd9-7d02dc48fa7d/scratchpad/. The suite was not
rerun, as instructed: it was green in round 1 and the code is unchanged.)

## Attack surface
- R1 midpoint routing (plan:64-67):
  - evenness in p at fixed r, against the library's own LoopIntegrand (Loops.wl:28-34) — numerical probe r3-probe-midpoint — **finding E25**
  - regulator: dsR is an abstract atom with no argument until Assemble (Loops.wl:31, FrequencyIntegration.wl:65, Assemble.wl:96) — **E26**
  - p → 0 constant slots (today's routing; Taylor-expanded F′, F″ atoms) — held
  - lower orders subtracted pointwise with un-averaged Taylor coefficients, n₀ ≥ 2 — held *once the symmetrisation is restored*: the even integrand minus its Taylor polynomial is O(p^{2n₀}), and these coefficients integrate to flow_n
  - p → 0 rows built from functional rows (R4 η row, Protected at order n₀, D15 regression) — probe r3-probe-shape-order and a hand toy — **E27**
  - two-family frequency integration: families s± = ar² + ap²/4 ± ar ap c; the coincident surface is c = 0 (and the whole line ar = 0); it is absent at d = 1 — held, provided the families are combined as in E3
  - Projection.LHS — held, but it is blind to routing: Γ2 has no loop variable, so the only check that certifies the subtraction is `Projection.functionalOrder`, and that check can be fooled through E26
  - KPZ-LO comparison "reference shifted to midpoint" — E25 (it also needs the r → −r symmetrisation)
  - D14 Oracle "(q, p, c)" — E25/E26 (it must use the same routing and symmetrisation, with ∂sR per line)
- R2 Related (plan:68-71):
  - pointwise d = 1 identity at generic p, against the library's loop — probe r3-probe-trs-pointwise — **held**: the identity holds at every q, raw and symmetrised, and the mutant η_D ≠ η_ν is detected
  - well-definedness in general: what α is, which side is eliminated, a side with no reference, the regulator, AtDimension, misuse — **E28**
- R3 "Dimensionless" accepts consistent functions — held (closes E17)
- R4 "AllMomenta", Protected at n₀, "Equation" -> nk Fnu — the Equation semantics hold (β_F̂(0) = (F̂(0) − 1)η, the analogue of `nk ank`); the tag semantics → **E29**; Mode inference vs today's code → **E30**
- R5 Functional rows, sibling ASSUMPTION — held, except the conflict with R4 → E29
- R6 functionalTerm on the Γ2 slot — held
- R7 hint table — incomplete for the round-2 checks → E32 (updated table below)
- R8 names/addressing — held. One addition (E25): rows must record their routing
- R9 — held
- Plan text vs R1 (D10, D11, D14, D15, cascade) — **E31**

## Findings

### E25 [blocker] (R1 resolution of E15 fails) Midpoint routing alone does not make the integrand even in p. The library's loop puts ∂sR on one line only, so the q ↔ −p−q symmetrisation that R1 drops is still required
- Where: plan:64-66 (R1: "the integrand is even in p at fixed r … no separate q ↔ −p−q symmetrisation"); Kernel/Loops.wl:28-34. `bubble = Tr[dsr . g0 . Γ3 . gm . Γ3 . g0]` has ∂sR at (w, q) only, and the p+q line carries none (Loops.wl:3-4, "Formula").
- Evidence (r3-probe-midpoint.wl):
  - Setup: the KPZ (hbar,h) bubble from `LoopIntegrand`, a momentum-dependent Γ2 (the same atoms as r1-probe-pointwise), d = 1, c = +1, ν = 0, ar = 0.7. ∂sR is attached at q (the correct attribution), then q = r − p/2 is substituted and ω is integrated numerically.

    | p | raw I/p² | (I(r,p)+I(r,−p))/2/p² | odd part (I(p)−I(−p))/2/p |
    |---|---|---|---|
    | 1e-1 | 0.1370 | −0.24076 | 0.03778 |
    | 1e-2 | 3.133 | −0.23782 | 0.03371 |
    | 1e-3 | 33.43 | −0.23779 | 0.03367 |
    | 1e-4 | 336.4 | −0.23779 | 0.03367 |

  - The raw midpoint ratio grows ∝ 1/p, exactly as the q-routed raw ratio did in r1.
  - Reason: by transposing the trace, the integrand with ∂sR on the other line is L(−p−q, p); in midpoint variables that is L(−r − p/2, p). Point inversion L(q, p) = L(−q, −p) holds for anything built from scal, and in d = 2 also from wedge. With it, the line-symmetric integrand is L_mid(r, p) + L_mid(r, −p), and this sum, not L_mid alone, is even in p.
- Consequence: as written, `Projection.functionalOrder` aborts on KPZ_LO (n₀ = 1). If the check were relaxed, F̂ would carry 1/ap pieces and the ŷ = 0 limit would be undefined. At d = 1 the values at c = +1 and c = −1 would differ.
  - The "KPZ-LO reference shifted to midpoint routing" comparison is pointwise only after the same r → −r symmetrisation on the reference.
  - The D14 Oracle at "(q, p, c)" compares against a differently routed integrand.
- Hint: keep the symmetrisation. In midpoint variables it is r → −r, which equals p → −p at fixed r. Then subtract the un-averaged Taylor coefficients. Apply r → −r to the transcribed reference too.
- Resolution: _(author)_

### E26 [major] dsR carries no line momentum until Assemble. After the (needed) symmetrisation its two copies sit at s− and s+, and `functionalOrder` cannot detect a wrong attribution
- Where: Kernel/Loops.wl:31 (abstract `dsR[a,b]`, implicitly at q²); FrequencyIntegration.wl:65 (memoised as an argument-free atom); Assemble.wl:96 (`DsRRules[reg, qU^2]`, one momentum for the whole row). plan:98 ("Loops — unchanged"), plan:107 ("Assemble — ScaleAtoms at s₂ too": dsR is not mentioned).
- Evidence: in L_mid(r, p) + L_mid(r, −p) the first term needs ∂sR(s−) and the second ∂sR(s+). With an argument-free atom both copies are the same symbol. Assemble then substitutes a single argument (today qU²), which is wrong for at least one copy.
  - The evenness/order check sees dsR as a p-independent atom, so it passes with the wrong attribution: the formal expression is even either way.
  - Today this cannot happen: there is one loop momentum, and dsR always sits on it.
- Consequence: a wrong functional integrand passes every named check proposed in the plan.
- Hint: give dsR its line argument (`dsR[a,b][s]`, or substitute `DsRRules[reg, s]` per line) in Loops/Projection before the symmetrisation. Have the Oracle evaluate ∂sR numerically on each line at (r, p, c).
- Resolution: _(author)_

### E27 [major] Every ŷ = 0 quantity taken from a midpoint row differs pointwise from today's p → 0 integrand and needs one more shape derivative: the D15 regression, the R4 η row and the R4 Protected-at-n₀ check break
- Where: plan:59-60 (D15 "p → 0 regression against today's KPZ equations with the c-integral done"); plan:76 (R4 "Equation" -> nk Fnu ⇒ η from Fnu's equation at ŷ = 0, a series limit, D7 [E1]); plan:75 (R4 Protected at order n₀); Kernel/Projection.wl:122-123 (`Projection.shapeSmoothness`: "C^{2n}"); tests/reference/KPZSetup.wl:5-7 (Litim: r''' and r'''' are lost in the termwise integrals).
- Evidence:
  - (a) Pointwise mismatch (hand toy, re-derivable in three lines). For L = B(|q+p|²) (∂sR line constant):
    - today's direction-averaged p² coefficient is B′ + 2q²B″/d;
    - the symmetrised midpoint one is B′/4 + r²B″/(2d).
    - With B = 1/(1+s), d = 1, r = q = 0.7 the two are 0.142 and 0.0355.
    - Both integrate to 0 over r ∈ ℝ: ∫(B′/4 + r²B″/2) = −π/8 + π/8, and ∫(B′ + 2r²B″) = −π/2 + π/2. They differ by a total r-derivative.
  - (b) Shape order (r3-probe-shape-order.wl): the KPZ (hbar,h) p² coefficient at fixed loop momentum, d = 1, before the ω integral, contains `rn/rd` derivatives of order {1} in today's routing and {1, 2} in symmetrised midpoint routing. Shifting the ∂sR line Taylor-expands ∂sR itself.
- Consequence:
  - The D15 regression holds pointwise only for n₀ = 0 (f_D). For f_ν it fails unless the integrals are done numerically.
  - The η row of `"Equation" -> nk Fnu` acquires r″ terms. With the Litim shape used by today's exact KPZ test (KPZAdimensional), these are δ terms that the termwise integration drops.
  - A true Ward identity tagged Protected at order n₀ ≥ 1 can abort. Today's checks are pointwise in q-routing, where the midpoint coefficient is that integrand plus a nonzero total derivative.
  - `Projection.shapeSmoothness` understates the requirement for midpoint lower-order subtractions.
- Hint: compute every ŷ = 0 object (the η row, the Protected-at-n₀ check, the regression, ShowBeta's "F̂(0)" line) with today's p → 0 operator. It is the same integral, and it uses the routing and measure of the constant rows. Keep midpoint only for the generic-p functional row, and extend the shapeSmoothness ASSUMPTION for subtracted orders.
- Resolution: _(author)_

### E28 [major] (R2 resolution of E16 is KPZ-specific) Related's α is undefined, the side being eliminated is ambiguous, and the "derived" relations hold only for ref×F ↔ ref×F with α̂ = 1
- Where: plan:68-71 (`Related[nk Fnu -> α Dk FD, "AtDimension" -> 1, …]`; "derives FD̂ = Fnû, η_D = η_ν and the regulator relation … abort if R violates it"); D9 hint (plan:38-39).
- Evidence:
  - Held part (r3-probe-trs-pointwise.wl): constant-coupling KPZ, d = 1, rn = rd, η_ν = η_D = 0.3, ν/D = 2, p = 0.6. I11/p² − (ν/D)(−I02/2) after the ω integral is ~1e-25 at q = 0.7 and 1.3, both raw and symmetrised. The mutant η_D = 0.5 gives −7.9e-3 and −3.0e-3. So the pointwise check at c = ±1 is sound.
  - (a) What α is. TRS gives nk(Fnu + rn) = T·Dk(FD + rd) with T fixed. Write α̂ := α Dk/nk; this must be a k-independent number, so the general content is:
    - Fnû = α̂ FD̂;
    - rn = α̂ rd (for the regulator written as in KPZSetup);
    - η_ν = η_D.

    "FD̂ = Fnû" is the case α̂ = 1. α is dimensionful in KPZ (nk and Dk are independent references), and its dimension, declaration and value are not specified. A numeric α in KPZ is inhomogeneous, and there is no named abort for it.
  - (b) Direction. In Mathematica `a -> b` replaces a: `nk Fnu -> α Dk FD` eliminates Fnu. D12's `FD -> α Fnu` eliminated FD ("Fnu stays in the state"). R2 flips the sides, yet still states FD̂ = Fnû without saying which equation is dropped.
  - (c) Other legal ansätze. D6 allows a function with no reference, e.g. an FDT-type `Flam hbar dt[h]` ↔ `Dk FD hbar²`. Then k-independence of α forces η_D = 0, not "η_D = η_ν". The rule is that the η's of the references on the two sides balance, and the plan does not state it.
  - (d) The regulator. With symbolic shapes (rn, rd) the relation cannot be checked, only imposed. "Abort if R violates it" is decidable only for the structural prefactors: the KS-style `Dk q² rn` needs rn = (nk/Dk)·rd, which is dimensionally impossible.
  - (e) Unspecified:
    - AtDimension with a symbolic Dim: is Dim → 1 substituted in every row, including the constant and η rows and the d = 1 measure?
    - AtDimension with FixedDimension[2].
    - With D7's default (F̂(0) free), the initial condition FD̂ = Fnû/α̂ that the downstream integrator must respect.
  - (f) Misuse with no named check:
    - a function related to a constant coupling;
    - both sides in one slot;
    - chains or cycles (A→B, B→C);
    - an eliminated function listed in a user "State".
- Consequence: an ansatz option that can be read two ways (which function survives, what α means). Outside KPZ-with-α̂ = 1 the derivation either does not exist or produces the wrong η relation.
- Hint: relate the *functions* and let the library derive the rest:
  - before: `"Related" -> {Related[nk Fnu -> α Dk FD, "AtDimension" -> 1, "Reason" -> "time reversal (d = 1)"]}`
  - after: `"Related" -> {Related[FD -> Fnu, "AtDimension" -> 1, "Reason" -> "time reversal (d = 1)"]}`. FD is eliminated (Rule semantics). An optional number `FD -> 2 Fnu` sets α̂; it is dimensionless, since both functions carry their references or are declared dimensionless (`Adimensional.relatedDimension` otherwise).
  - The library derives, prints as ASSUMPTION lines, and imposes:
    - the η balance from the two slots' references (η_D = η_ν here; η_D = 0 if one side has none);
    - the shape relation from the Regulator stage (imposed; abort `Adimensional.relatedRegulator` only on a prefactor mismatch);
    - Dim → 1 in every row (abort `Adimensional.relatedDimension` if the Theory fixed d ≠ 1).
  - It checks β_FD̂ = α̂ β_Fnû pointwise at c = ±1.
- Resolution: _(author)_

### E29 [minor] (R4 vs R5) Whether a tag matches order n₀ of a functional slot is contradictory, and what "F(0) does not flow" means is ambiguous for a referenced function
- Where: plan:74-75 (R4: `"Momentum" -> n<… ≤ n₀>` at order n₀ means "F(0) does not flow"; the text is garbled); plan:77 (R5: "tags do not match there"); Symmetries.wl:36 (an omitted "Momentum" is `All`, Symmetries.wl:25).
- Evidence:
  - R5 says no tag matches orders ≥ n₀. R4 says a tag at n₀ does.
  - `Protected[{hbar,h}, "Frequency" -> 0, "Reason" -> "shift"]` with no Momentum is `All`, which contains n₀ = 1 of Fnu's slot. One reading silently imposes "Fnu(0) does not flow", and the other ignores the tag.
  - For `nk Fnu`, "F(0) does not flow" can mean ∂s(nk F(0)) = 0 (the loop part vanishes, β_F̂(0) = η F̂(0)) or β_F̂(0) = 0. With a Fixed nk these differ by η F̂(0), which is r2 E18(b) again.
  - Neglected at n₀, or Protected at n > n₀, has no defined meaning and no check name.
- Hint: only an explicit integer `"Momentum" -> n₀` matches order n₀ (All and patterns never do). Its meaning is "the loop part of the slot's flow at p = 0 vanishes", checked with today's operator (E27). Everything else at ≥ n₀ aborts `Symmetries.functionalTag`.
- Resolution: _(author)_

### E30 [minor] R4's "Mode inferred from Equation" contradicts today's code, which silently ignores "Equation" without "Mode" -> "Running" (earlier-step defect)
- Where: Kernel/Adimensional.wl:152-157 (eqSpec is read only in the `mode === "Running"` branch; otherwise the reference is Fixed); Kernel/Dimensions.wl:110 (accepts "Equation" without "Mode"); plan:76.
- Evidence: today `nk -> {"nu", "Equation" -> nk ank}` makes nk Fixed with no message. Under R4, `nk -> {"nu", "Equation" -> nk Fnu}` makes it Running. One key would mean two things depending on whether the amplitude is constant.
- Hint: infer Running from "Equation" in both cases (or abort `Dimensions.equationNeedsRunning`). Abort `Dimensions.equationModeFixed` on "Mode" -> "Fixed" with "Equation".
- Resolution: _(author)_

### E31 [minor] The plan body still describes the superseded q-routing; an implementer following the cascade gets something else
- Where: D10 (plan:40-41), D11 (plan:42-44: variables (aq, ap, c), F̂(ap² + aq² + 2 ap aq c)), D14 (plan:56, "(q, p, c)"), D15 (plan:58-60), Theory (plan:88: reserves ap, c, not ar), Projection (plan:99-100), FrequencyIntegration (plan:102: s₁ = q², s₂ = (p+q)²), Assemble (plan:107: "(q, p, c)"), Docs (plan:117: CONVENTIONS lacks ar and the routing).
- Hint: rewrite these in place in design/MOMENTUMFUNCTIONS.md. "R supersedes D" is not acceptable in a record that will become authoritative. Record the routing per row (`"Routing" -> "Midpoint"`, `"LoopVariable" -> ar`) so that ShowBeta can print "lines at |r ± p/2|".
- Resolution: _(author)_

### E32 [minor] (R7) The check → hint table does not cover the round-2 additions
- Where: plan:80-81. Missing: `Adimensional.related`, `Adimensional.relatedRegulator`, `Symmetries.siblingAtZeroMomentum`, the functional tags (E29), the Equation/Mode misuse (E30), and the Related misuses (E28 f).
- Hint: use the updated table below.
- Resolution: _(author)_

## Usability (objective 2): before → after
- **Related** (E28)
  - before: `AdimensionalEquations[eqs, dims, "Related" -> {Related[nk Fnu -> α Dk FD, "AtDimension" -> 1, "Reason" -> …]}]`
  - after: `"Related" -> {Related[FD -> Fnu, "AtDimension" -> 1, "Reason" -> "time reversal (d = 1)"]}`
  - Card:
    `ASSUMPTION Adimensional.related: FD eliminated: aFD(y) = aFnu(y); etaD = etanu; rd = rn; Dim -> 1 (time reversal (d = 1)); initial data must satisfy aFD = aFnu`
- **Tags on a functional slot** (E29)
  - before: `Protected[{hbar,h}, "Frequency" -> 0, "Reason" -> "shift"]` (ambiguous at n₀)
  - after: `Protected[{hbar,h}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "shift"]` (as in KPZSetup today)
  - On All, the FAIL text reads: "Momentum -> All on Fnu's slot: write the orders below n₀ = 1 explicitly; order n₀ is Fnu(0) (Protected only with "Momentum" -> 1)"
- **Equation** (E30)
  - `nk -> {"nu", "Equation" -> nk Fnu}` and `nk -> {"nu", "Equation" -> nk ank}` both mean Running, with the amplitude counted as 1. "Mode" is never needed with "Equation".
- **Outputs** (E25, E31)
  - `ShowBeta[adim, Fnu]` prints `∫ dar dc K_d ar^(d−1)(1−c²)^((d−3)/2) · BetaNL(ar, ap, c), lines at s± = ar² + ap²/4 ± ar ap c (symmetrised r → −r); aFnu(0) row from the p → 0 operator (aq)`
  - The user can then see that the function row and the η row use different loop variables, and why.
- **Held**: the DefineAnsatz line for KPZ LO, `"MomentumFunctions"`, the R3 "Dimensionless" acceptance, `"AllMomenta"` as a string, `EquationOf[eqs, Fnu]` / `ShowBeta[adim, Fnu]`, `aFnu` naming, ShowMonomials `Fnu(p²)`.

## Nonsensical / ambiguous ansatz cases (round 3)
| case | named check | hint |
|---|---|---|
| function in a Γ3/Γ4 monomial | Ansatz.functionQuadraticOnly | "functions multiply quadratic terms only (D1)" |
| function under an operator | Ansatz.functionNotInOperator | "write `Fnu lap[h]`" |
| undeclared / also in "Couplings" / is a field | fresh-symbol, knownSymbols, couplingsNotFields | existing |
| `Fnu[0]` | Ansatz.functionBare | "write Fnu bare" |
| `Fnu^2`, `1/Fnu`, F at two powers, `nk Ak Fnu`, `(nk+Dk) Fnu` | Vertices.functionalTerm (Γ2 level, R6) | "one term ref × Fnu × (p²)^n₀"; quote the slot |
| constant at ≥ n₀, or a second function in the slot | Symmetries.functionalSlot | as in r2 |
| function in 0 or 2 slots | Symmetries.functionSlot | 2 slots → "use two functions and Related[FD -> Fnu]" (E28 syntax) |
| lower orders not declared, not tagged | Symmetries.allTagged (existing) | existing |
| non-O(p^{2n₀}) after subtraction | Projection.functionalOrder | **false abort on KPZ_LO unless the symmetrisation is restored (E25)**; blind to dsR attribution (E26) |
| function in the regulator | Regulator.noFunctions | as in r2 |
| NonRenormalised[Fnu] | Symmetries.nonRenormalisedFunction | "use Protected at "Momentum" -> n₀" |
| function as a reference | Dimensions.functionReference | as in r2 |
| "Dimensionless" contradicting the derived dimension | Dimensions.functionDimensionless (R3) | held |
| "Equation" naming a function not multiplying that reference / a dimensionful function | Adimensional.functionEquation | "Equation -> nk Fnu needs Fnu in nk's slot term" |
| "Equation" with "Mode" -> "Fixed"; "Equation" without Mode on a constant amplitude | **none** (E30) | proposed Dimensions.equationModeFixed / equationNeedsRunning |
| user "Pipeline" contradicting the functions | Ansatz.pipeline (adapted) | as in r2 |
| Protected/Neglected "AllMomenta" on a functional slot | Symmetries.contradiction (hint updated, R7) | held |
| tag with `All`/pattern reaching order n₀; Neglected at n₀; Protected at n > n₀ | **contradictory R4/R5** (E29) | proposed Symmetries.functionalTag |
| catch-all Neglected[_] at ≥ n₀ | does not match (R5) | held |
| sibling slot left at p → 0 | ASSUMPTION Symmetries.siblingAtZeroMomentum (R5) | held |
| Related: numeric α between independent references; ambiguous eliminated side | **none** (E28 a, b) | proposed Adimensional.relatedDimension; Rule semantics |
| Related on a side with no reference (FDT type) | **wrong derived η relation** (E28 c) | η balance |
| Related with a structurally incompatible regulator | Adimensional.relatedRegulator | only on a prefactor mismatch; shapes imposed (E28 d) |
| Related with AtDimension 1 and a Theory at d = 2; symbolic Dim left in the output | **unspecified** (E28 e) | Adimensional.relatedDimension; Dim → 1 |
| Related function ↔ constant, same slot, chain/cycle, eliminated function in "State" | **none** (E28 f) | Adimensional.relatedForm / relatedState |
| valid p⁴ function over a constant p² | passes after E25's fix (held in principle; the shape must be C³ for the subtraction, E27) | — |
| f^λ (no reference) | passes (R3) | held |
| quadratic `jac` with a function | Ansatz.jacQuadratic | existing |
| function absent from the density | unspecified (r2) | WARN like Ansatz.unusedCouplings |

## Appendix: probes (`math -script <file>`)

### r3-probe-midpoint.wl
```mathematica
(* R1: is the library's KPZ loop integrand even in p at fixed midpoint r = q + p/2 (d = 1, c = +1, nu = 0)? *)
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
quiet[x_] := Block[{Print}, x];
Get["/home/greg/Physics/Sources/FRGscalarMSR/tests/reference/KPZSetup.wl"] // quiet;
{wv, qv} = LoopVariables[kpzLoops];
li = LoopIntegrand[kpzLoops, {hbar, h}, 0, pp];
li = li /. DsRRules[kpzRegulator, scal[qv, qv]];          (* d_s R attached to its line (q) before the shift *)
li = li //. DeltaRules[kpzPropagator];
atoms = {GammaAtom[h, hbar, 0] -> (# /(1 + #) + # Exp[-#] &), GammaAtom[hbar, h, 0] -> (# /(1 + #) + # Exp[-#] &),
   GammaAtom[hbar, hbar, 0] -> (-2/(1 + 2 #) - 2 Exp[-#] &)};
li = li /. atoms /. AtomRules[kpzPropagator];
num[r_, p_] := li /. scal[a_, b_] :> ((a /. {qv -> r - p/2, pp -> p}) (b /. {qv -> r - p/2, pp -> p})) /.
   {rn -> (Exp[-#] &), rd -> (Exp[-#] &), nuk -> 1, Dk -> 1, lamk -> 1, k -> 1, dsCoupling[_] -> 0};
wint[r_, p_] := With[{w = wv, e = num[r, p]}, NIntegrate[e, {w, -Infinity, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 12]];
Do[With[{a = wint[7/10, p], b = wint[7/10, -p]},
   Print["p = ", N[p], "  midpoint raw I/p^2 = ", N[a/p^2, 6], "   (I(p)+I(-p))/2/p^2 = ", N[(a + b)/2/p^2, 6], "   odd part (I(p)-I(-p))/2/p = ", N[(a - b)/2/p, 6]]],
 {p, {1/10, 1/100, 1/1000, 1/10000}}];
```
Output: raw 0.1370, 3.133, 33.43, 336.4; symmetrised −0.24076, −0.23782, −0.23779, −0.23779; odd/p 0.03778, 0.03371, 0.03367, 0.03367.

### r3-probe-shape-order.wl
```mathematica
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
quiet[x_] := Block[{Print}, x];
Get["/home/greg/Physics/Sources/FRGscalarMSR/tests/reference/KPZSetup.wl"] // quiet;
{wv, qv} = LoopVariables[kpzLoops];
li = LoopIntegrand[kpzLoops, {hbar, h}, 0, pp] /. DsRRules[kpzRegulator, scal[qv, qv]] //. DeltaRules[kpzPropagator];
li = li /. AtomRules[kpzPropagator] /. k -> 1;
at[q_, p_] := li /. scal[a_, b_] :> ((a /. {qv -> q, pp -> p}) (b /. {qv -> q, pp -> p}));
shapeDerivs[e_] := Union[Cases[e, Derivative[n_][rn | rd] :> n, {0, Infinity}, Heads -> True]];
today = D[at[x, p], {p, 2}] /. p -> 0;
mid = D[(at[x - p/2, p] + at[-x - p/2, p])/2, {p, 2}] /. p -> 0;
Print["today's routing: ", shapeDerivs[today], "   symmetrised midpoint: ", shapeDerivs[mid]];
```
Output: `{1}` vs `{1, 2}`.

### r3-probe-trs-pointwise.wl
```mathematica
(* R2: d = 1 time reversal of the Gamma2 flows, pointwise in q at generic p (constant-coupling KPZ, rn = rd, eta_nu = eta_D) *)
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
quiet[x_] := Block[{Print}, x];
Get["/home/greg/Physics/Sources/FRGscalarMSR/tests/reference/KPZSetup.wl"] // quiet;
{wv, qv} = LoopVariables[kpzLoops];
prep[pair_] := Module[{li = LoopIntegrand[kpzLoops, pair, 0, pp]},
   li = li /. DsRRules[kpzRegulator, scal[qv, qv]] //. DeltaRules[kpzPropagator]; li /. AtomRules[kpzPropagator]];
nuv = 2; Dv = 1; eta = 3/10;
pars = {rn -> (Exp[-#] &), rd -> (Exp[-#] &), nuk -> nuv, Dk -> Dv, lamk -> 1, k -> 1,
   dsCoupling[nuk] -> -eta nuv, dsCoupling[Dk] -> -eta Dv, dsCoupling[lamk] -> 0};   (* mutant: dsCoupling[Dk] -> -(eta + 1/5) Dv *)
l11 = prep[{hbar, h}] /. pars; l02 = prep[{hbar, hbar}] /. pars;
num[e_, q_, p_] := e /. scal[a_, b_] :> ((a /. {qv -> q, pp -> p}) (b /. {qv -> q, pp -> p}));
wint[e_] := With[{w = wv}, NIntegrate[e, {w, -Infinity, Infinity}, WorkingPrecision -> 25, PrecisionGoal -> 12]];
f[q_?NumericQ, p_?NumericQ] := wint[num[l11, q, p]]/p^2 - (nuv/Dv) (-wint[num[l02, q, p]]/2);
fs[q_?NumericQ, p_?NumericQ] := (f[q, p] + f[-p - q, p])/2;
p0 = 3/5;
Print["raw: ", N[{f[7/10, p0], f[13/10, p0]}, 6], "  sym: ", N[{fs[7/10, p0], fs[13/10, p0]}, 6]];
```
Output: all ~1e-25 (0 to working precision); the mutant gives raw {−0.00794, −0.00296} and sym {−0.00095, −0.00111}.

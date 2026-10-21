# the-enemy — plans/look-at-sync-cnrs-projects-shallow-water-twinkly-candle.md, round 1
Attacked commit: 6208fe0 (FRGscalarMSR; README.md modified, uncommitted)  ·  Date: 2026-10-21
Verdict: 2 blocker, 4 major, 5 minor, 1 nit

Plan mode was active, so this report could not go to `reviews/<target>/r1.md`. The evidence is given as inline
commands, rerunnable read-only. Baseline: `math -script tests/runTests.wl` → 343/343 passed, exit 0.

Paths: P = ~/.claude/plans/look-at-sync-cnrs-projects-shallow-water-twinkly-candle.md, SW = ~/sync_cnrs/Projects/shallow_water.

## Attack surface
- P:52 monomial dimension equation: independent solve of the LPA-plus system with both signs (E-a) — finding E1
- P:17 D1, "N+1 freedoms": counted by hand for LPA-plus (Ak redundant, aDphik needed) and DE4plus (aDhk needed, ank comes out 0) — held
- P:18 D2: DE4plus needs {Ak, aDhk, aDphik} declared; the abort-and-name rule covers this — held
- P:19 D3 Leibniz, NonRenormalised ⇒ ∂s c = 0: hand derivation for Hk·lk and nk·Ak — held
- P:20 D4 against the user's real LPA regulator — finding E2
- P:21 D5/D6 — held
- P:23 D7.1–2, 3, 5 — held; D7.4 HomClean ⇒ unit-free — held; D7.6 — finding E6; the extra random rescaling — finding E8
- P:24 D8 against DE4plus — finding E5
- P:25 D9 (old dsMh ∝ η_φ): checked against the old notebook rule `dsMh -> -2*Dhk*q^2*(-etaphi*armh…)`, and ∂s Dhk = −η_h Dhk — held; η_h = η_φ neutralises it in betal and betanl — held
- P:26 D10 naming — finding E10
- P:27 D11 signs: running η = −RHS/ref gives Ml = I, Mnl = L, bnl = −L0, the same as the old `(Id+∫M)η = −∫b`; fixed β_a = aη + RHS/ref, the same as the old `Ak*etanu` — held. Classification applied to DE4plus — finding E3
- P:28 D12 against FlowRG `build_model` (fixed_etas accepts any value, rows masked) — held
- P:29 D13, P:32 D16 against equations-interface.md — held; P:73 "omit all-zero parts" — finding E12
- P:30 D14 — finding E11
- P:65 measure k^d S_d aq^{d−1}/(2π)^d: rederived, correct for isotropic integrands — held (S_d undefined — E10)
- P:63–68 the l/nl split, the linear part −x + Σyη: checked against the old LinEQ and `rhsHEQl = aH(−1+etah/2−etanu/2)` — held
- P:99–110 Levi-Civita section: EpsIndex canonical form (Ansatz.wl:89), theory["Dimension"] = 2, the F4 isotropy (Projection.wl:26–42) — held
- P:113 hand values (old Dimensionality): g, H, l reproduced with the correct sign — held
- P:117 KPZ acceptance reference — finding E4
- P:117–119 DE4plus acceptance — findings E3, E5, E7
- P:119, P:129 LPA end-to-end / FlowRG load — finding E2
- Monomial coefficients that are sums — finding E9

## Findings

### E1 [blocker] The sign of the measure in the per-monomial dimension equation is wrong
- Where: P:52 "[coef] + Σ[fields] − nt·[t] + n_x − [measure] = 0 (measure = [t] − Dim)". The plan says this goes into CONVENTIONS.md.
- Evidence: Γ = ∫dt dᵈx 𝓛 is dimensionless (P:9), so [𝓛] = −[measure]. The equation therefore needs **+[measure]**. The LPA-plus system
  (SW/ShallowWaterLPA.nb ansatz; refs nk, Dhk, Dphik; Ak, aDphik dimensionless) was solved both ways, with vectors over (k, nk, Dhk, Dphik):
  - sign +: g {1,1,−1/2,1/2}, H {Dim/2,−1/2,1/2,0}, l {1−Dim/2,3/2,0,−1/2}, [H̄]+[H] = Dim. This equals the old notebook's Dimg, DimH, Diml.
  - sign − (as written in P:52): g is unchanged, but H {(−4−Dim)/2,−3/2,1/2,0}, l {(6+Dim)/2,5/2,0,−1/2}, [H̄]+[H] = −4−Dim−2·nk. These are wrong.
  - Rerun: the `math -noprompt` heredoc `solve[sign_]` (12 monomial equations of the LPA-plus density; see the reply).
- Consequence: every field and coupling dimension except a few (g, mk, Dmk) is wrong. The wrong convention would be written into
  CONVENTIONS.md. Note that a CheckDimensions on gk alone, as in the workflow example P:44, passes under the wrong sign.
- Hint: write "+ [measure]", or define the measure weight as d − [t].
- Resolution: _(author)_

### E2 [blocker] D4 forbids the regulator of the plan's own main target (LPA-plus, swLPAplusequations.wdx)
- Where: P:20 (D4) vs P:7 (context: AdimensionalLPAsw.nb "is waiting for this"), P:119 (LPA end-to-end from swLPAplusequations.wdx), P:129 (load in FlowRG).
- Evidence: SW/ShallowWaterLPA.nb:836 has `{Hbar, H} -> nk (1 + Ak) q^2 arnh[q^2/k^2]`. Importing SW/swLPAplusequations.wdx gives
  DsRRules containing `nk*q^2*arnh[q^2/k^2]*dsCoupling[Ak]`, and the DsCouplings of the equations are {dsCoupling[Ak], dsCoupling[Dphik], dsCoupling[nk]}.
  Ak is not a reference, so D4 aborts at the first step of AdimensionalEquations.
- Consequence: the acceptance in step 4, the verification at P:129 and the user workflow the plan was written for cannot run. The case that
  is needed is exactly the FUTURE.md item "non-reference regulator couplings via X = (η…, ∂sc…)" that the plan defers.
- Hint: either implement the joint linear solve for (η, ∂s c_reg) now, or change the acceptance and context to a regulator with reference-only prefactors and say so.
- Resolution: _(author)_

### E3 [major] D11 classification on DE4plus (and the old LPA) cannot reproduce the old η systems, and the acceptance only compares β
- Where: P:27 (D11), P:117–118.
- Evidence: in the library DE4plus ansatz (tests/reference/ShallowWaterSetup.wl), nk occurs only as nk·Ak and nk·ank, Dhk only as Dhk·aDhk,
  and Dphik only as Dphik·aDphik and Dphik·D2phik. None of them occurs alone in a flowing combination, so all three would be classified Fixed.
  - Read strictly, "two candidate equations ⇒ abort" applies to nk and Dphik, and DE4plus aborts.
  - Otherwise MUST_FIX = (:nu,:phi,:h) and the η system is three placeholder rows.
  - The old SW/ManualCodeOld/DE4plus_largscale_forcing_LLM/eqs/Equations2D.jl solves a dynamic η system from the nan, DphiaDphi and DhaDh rows
    (`etaMnl`, `etabnl`). In effect it keeps ank, aDphik and aDhk at 1: β_ank = (ank − 1)η_ν.
  - The acceptance compares only betal and betanl, so this η-system regression goes unchecked. Neither reading is stated.
  - The same holds for the library LPA test ansatz (ank, aDhk) against LPA_KPZ_thermal_forcing (etanuexp etc.).
- Consequence: the DE4plus example and acceptance either abort or produce a file whose η's must all be hand-fixed. That is a different
  physics prescription from the reference, and nothing flags it.
- Hint: state the DE4plus classification and MUST_FIX explicitly, and either add "η fixed by amplitude ≡ 1" or compare the η system and record the expected difference.
- Resolution: _(author)_

### E4 [major] The KPZ acceptance reference is not a KPZ equation set
- Where: P:117 "KPZ vs ManualCodeOld/LPA_KPZ_thermal_forcing/eqs (full)".
- Evidence: `eqs/Equations2Dgeneric.jl` holds the 2-field shallow-water LPA, with couplings A, ag, aH, al and functions rhsgEQnl, rhsAEQl,
  etanuexp, …. It uses the old export format (alphabetical arguments, no COUPLINGS/ETAS). The folder's CLAUDE.md §1–2 says so. The 1-scalar
  KPZ references are Canet Eq. (20) (tests/acceptance/KPZ.wl) and `lpa/kpz.jl`.
- Consequence: a "full" comparison is undefined, or would be engineered on a subspace (ag = aH = 0) without stating it.
- Hint: compare KPZ with Canet / kpz.jl. Compare the library LPA ansatz with Equations*Dgeneric.jl through the documented old-format bridge.
- Resolution: _(author)_

### E5 [major] D8 makes the planned DE4plus export, example and round trip impossible
- Where: P:24 (D8), P:80–81 (export section in ShallowWaterDE4plus example), P:117–118, P:30.
- Evidence: the old DE4plus Equations2D.jl contains D2armh ×12, D2armphi ×12, D2arnh ×16 and D2arnphi ×16 (from the g2 and D2phi flows).
  docs/FUTURE.md says "g2 and D2phi need r″". (KPZ and LPA-plus integrands contain only first derivatives, checked.)
- Consequence: ExportJulia aborts on DE4plus. The example notebook section fails, and the D14 round trip cannot cover DE4plus. The step 4
  comparison against a .jl has no defined route; the old file itself does not load in FlowRG.
- Hint: drop DE4plus export from the examples and tests, and state that the DE4plus acceptance is done in Mathematica on the Adimensional stage after parsing the old .jl.
- Resolution: _(author)_

### E6 [major] Check D7.6 states an invariance that is false
- Where: P:23 (6) "numeric rescaling of x, t, fields leaves the density invariant"; P:55, P:103.
- Evidence: with the correct weights (E1), the KPZ/LPA monomial H̄ ∂tH scales as λ^{[H̄]+[H]−[t]} = λ^{Dim+2} (k-part) ≠ 1.
  Only dt dᵈx 𝓛 is invariant, and only when the couplings are rescaled as well. With the couplings fixed, the monomials do not even scale alike.
- Consequence: as written, the "independent" check fails, or gets adjusted until it passes. It also cannot catch E1 unless it includes the Jacobian of the measure.
- Hint: rescale x, t, fields and couplings by their weights; assert the density changes by exactly λ^{d−[t]} (i.e. the action is invariant).
- Resolution: _(author)_

### E7 [minor] The (2π)^{−d} acceptance outcome is left open but can already be decided
- Where: P:118–119 "this also settles whether the old file had the (2π)^{-d}".
- Evidence:
  - The old notebook has `IntegrationPrefactor[d_] := 2*Pi^(d/2)/Gamma[d/2]*aq^(d-1)`, with no (2π)^{−d}.
  - The 2D file has `juliaPi*(…)` prefactors.
  - DE4plus_largscale_forcing_LLM/flowsFAST/*.wdx have the same md5 as tests/reference/flowsFAST, which the library Integrand (per dᵈq/(2π)ᵈ) matches exactly.
- Consequence: an acceptance criterion with an open outcome invites tuning a factor afterwards.
- Hint: state up front that betanl_new = (2π)^{−d} betanl_old exactly; the comparison must hit that ratio.
- Resolution: _(author)_

### E8 [minor] The "random rescaling of k and refs" check is tautological if applied to the final β
- Where: P:23 "Plus a random rescaling of k and refs leaving every β unchanged".
- Evidence: check 7 / HomClean already asserts that β is free of k and refs. Substituting k → λk then compares an expression with itself.
- Hint: apply it to the substituted expression before HomClean extracts the unit monomial (ratio to [c] times λ-power).
- Resolution: _(author)_

### E9 [minor] "[coef]" is undefined when a monomial coefficient is a sum
- Where: P:52.
- Evidence: `DefineAnsatz[… hbar (dt[h] - nuk (1 + Ak) lap[h] …)]` gives the coefficients {−nuk − Ak nuk, 1, −Dk, −lamk/2}. The parser keeps sums.
- Hint: one homogeneity equation per coefficient term (this forces [Ak] = 0), or abort with a hint.
- Resolution: _(author)_

### E10 [minor] Loose declarations
- Where: P:26, P:65.
- Evidence:
  - ank comes out with dimension 0 in DE4plus but is not declared. "a<>name for dimensionful" makes it unclear whether it becomes `aank` or `ank`; the old file uses `ank`.
  - S_d is never defined (should be 2π^{d/2}/Γ(d/2)).
  - Nothing says where DeriveDimensions gets the Ansatz when only a .wdx is loaded. It is embedded at equations["Closure"]["Integrated"]["Projection"]["Loops"]["Vertices"]["Ansatz"], but the plan names no accessor.
- Resolution: _(author)_

### E11 [minor] The D14 round trip can be skipped silently while the suite exits 0, and the FlowRG location is undeclared or stale
- Where: P:30, P:125; design/FRGDESIGN.md:9 (`~/sources/FlowRG`, which does not exist; the real path is ~/Physics/Sources/FlowRG.jl);
  examples/makeNotebooks.wl:15, FRGscalarMSR.wl:5 and tests/runTests.wl:4 (`~/sources/FRGscalarMSR`, which does not exist).
- Consequence: if FlowRG is looked up on the stale path, the round trip is never run and is still counted green.
- Hint: declare the FlowRG path once, and count a skip as a visible non-pass in the summary.
- Resolution: _(author)_

### E12 [nit] "omit all-zero parts" (P:73) should exempt betal, which FlowRG requires (equations-interface.md rule 3).
- Resolution: _(author)_

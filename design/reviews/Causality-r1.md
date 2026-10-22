# the-enemy — /home/greg/.claude/plans/i-want-to-add-streamed-fairy.md, round 1
Attacked commit: 1b9fae1  ·  Date: 2026-10-22
Verdict: 1 blocker, 3 major, 5 minor, 0 nit

(Plan mode was active: this report could not be written to `reviews/i-want-to-add-streamed-fairy/r1.md`. The evidence is
given as an inline `math -noprompt -run` probe, E1-probe below, so that it can be rerun.)

## Attack surface
- Confirm green: the target is a plan, so it has no machine checks of its own. Nothing to run.
- Context: "degree 2 → c1/c2>0, c0/c2>0; degree 1 → c0/c1>0" (Kernel/Propagator.wl:86-91). I re-derived it by hand with u = iω:
  Im ω > 0 ⇔ Re u < 0, so a real quadratic needs all coefficients of the same sign. The atoms are real (`Propagator.atomsReal`). — held
- KPZ claim `nuk q²(1+rn) > 0 → 1 + rn(aq²)`. Probe output: `{nuk*q^2*(1 + rn[q^2/k^2]) > 0}`, c = {GammaAtom, 1, 0}. — held as the Explicit form, but see E7
- SW claim, condition 1: in DE4 it reduces to `Ak + ank + arnh + arnphi` (probe). — held for DE4; wrong file, see E2
- SW claim, condition 2. — finding E1
- Test placement in `ShallowWaterLPAAdimensional.wl`. — finding E2
- Decision 6 / step 5: `expectFail[…, "Adimensional.validityMinimal"]` "on a corrupted condition". — finding E3
- Decision 5 / handoff: "> 0 for every aq in the domain", "evaluated on its aq quadrature nodes", compared with the ASSUMPTION "at every q". — finding E4
- Algorithm 9.1 and 9.3 (FactorTermsList, numeric/aq content). — finding E5
- Decision 8: positive constants dropped; omitted function versus the round-trip `z(:causality, k, …)`. — finding E6
- Decision 3, references > 0. Adimensional already uses fractional powers of the references (e.g. Dphik^(1/2)), and ∂s ref = −η ref keeps their sign. — held (implicit already; making it explicit is right)
- Decision 5: FlowRG accepts an extra function. equations-interface.md rule 6 says "may define anything else", and load.jl only requires REQUIRED_CONSTS and betal. — held
- Signature `causality(couplings…, aq)` without etas: the conditions come from Δ, which has no ∂sR. — held
- Step 3: `noUnits`/`knownSymbols`/`systemHeads` cover the new part through `all`. The validity atoms are ScaleAtom[…,0,0], so they contain no shape derivative. — held
- Decision 7, name "Causality". — finding E8
- Handoff "exposes the minimum margin". — finding E9

## Findings
### E1 [blocker] The quoted SW condition 2 drops a relative aq² (and uses the wrong coupling names)
- Where: plan "Current conditions", SW line; step 5, the ShallowWater test "gives the two conditions above".
- Evidence: E1-probe output, DE4 setup (`tests/reference/ShallowWaterSetup.wl` with the `de4Dims` of ShallowWaterDE4Adimensional.wl:14). Units are substituted and the result is Factor'ed:
  `aq^2*k^4*nk^2*(agk*aHk*alk + ag2k*aHk*alk*aq^2 + Ak*ank*aq^2 + ank*aq^2*arnh[aq^2] + Ak*aq^2*arnphi[aq^2] + aq^2*arnh[aq^2]*arnphi[aq^2])`
  So the minimal form is `agk aHk alk + aq² (ag2k aHk alk + (Ak+arnh)(ank+arnphi)) > 0`. The plan writes
  `ag al aH (…) + (A+arnh)(an+arnphi)`, which has no aq² in front of the second term. A by-hand dimension count gives the same result:
  [gk Hk lk q²] = k⁴nk² aq², while [nk² q⁴] = k⁴nk² aq⁴. The names are also wrong: the probe gives agk, aHk, alk, ank, ag2k, not ag, al, aH, an.
- Consequence: an "exact" acceptance test written from the plan fails. If it is instead tuned until it passes, the reference stops being independent.
- Hint: take the expected values from a hand dimension count (above). Do not copy them from old notation.
- Resolution: _(author)_

### E2 [major] The SW test is placed in the LPA file, but the conditions quoted are those of the DE4 ansatz
- Where: plan step 5, `tests/acceptance/ShallowWaterLPAAdimensional.wl`.
- Evidence: E1-probe, LPA ansatz of that file (it has `mk dt[H]` and no ank or g2k):
  `{aq^2 k^2 nk (1 + Ak + aHk alk amk + arnh + arnphi), aq^2 k^4 nk^2 (agk aHk alk + aq^2 (Ak + arnh)(1 + arnphi))}`.
  The first condition has an `amk aHk alk` term, and `ank` is replaced by 1. Neither matches the plan.
- Consequence: the acceptance criterion is false for the file it names.
- Hint: put the DE4 expectations in ShallowWaterDE4Adimensional.wl and the LPA ones (above) in the LPA file.
- Resolution: _(author)_

### E3 [major] The planned expectFail cannot trigger `Adimensional.validityMinimal`
- Where: plan step 2 (the check) and step 5 ("expectFail … on a corrupted condition").
- Evidence: a corrupted Equations stage is rejected first by `Stage.fingerprint` (Kernel/Core.wl:110). If the stage is re-stamped
  (as VerticesTest.wl:31 does), the stripping still yields stripped·minimal = original for whatever input it gets. The identity
  is an internal consistency check, so no input can make it fail. The mutation pattern in the repo replaces a private helper
  instead (ProjectionTest.wl:60, `Block[{FRGscalarMSR`Private`projCAverage = mutAvg}, …]`). The second half of the check,
  "stripped a positive monomial", has no planned mutation at all.
- Consequence: the named check ships without any test that shows it can fail. This violates the CLAUDE.md "New check" rule in substance.
- Hint: give the stripping helper a distinctive name and mutate it with Block twice: once with a wrong factor, once with a negative stripped content.
- Resolution: _(author)_

### E4 [major] The domain of the exported condition is a silent weakening of the declared assumption
- Where: plan Decision 5 ("for every aq in the domain") and the handoff ("on its aq quadrature nodes"). Compare
  Kernel/Propagator.wl:92 and Kernel/Assemble.wl:126 ("at every q"), docs/MODULES.md:258,534, and README.md:427.
- Evidence: FlowRG's domain is the integration support of the widest regulator (FlowRG.jl src/load.jl:154), which is (0,1) for
  Litim. Gauss nodes never include the aq → 0 end, and that end is where condition 2 reduces to `agk aHk alk > 0` (E1). A sign
  change between two nodes also goes unseen.
- Consequence: three statements of one assumption disagree ("every q", "domain", "nodes"), and no recorded decision says which
  holds or why. A violation outside the support, or between nodes, is reported as causal.
- Hint: record the reason in design/CAUSALITY.md. The residue formula is only needed where ∂sR ≠ 0. Then state that the check samples and is not a proof, and align the ASSUMPTION texts.
- Resolution: _(author)_

### E5 [minor] Algorithm 9.3 as written can strip a negative number and misses the aq content
- Where: plan Decision 9.1 and 9.3.
- Evidence: E1-probe prints `FactorTermsList[-2x-2y] = {-2, x+y}` and `FactorTermsList[-2 x aq^2 - 4 y aq^4] = {-2, aq^2 x + 2 aq^4 y}`.
  The numeric content can be negative, and FactorTermsList does not extract the aq power. Rule 9.1 does not cover a c_j that is
  a negative number. The "stripped positive" check would catch the sign error, if that check is tested (E3).
- Hint: strip |numeric content| and aq^(min exponent), with exponents counted outside the shape arguments rn[aq²].
- Resolution: _(author)_

### E6 [minor] Constant conditions that are zero or negative are not specified, and the round-trip count depends on where constants are dropped
- Where: Decision 8, step 3, step 4.
- Evidence: only "positive-constant conditions dropped" is specified. A constant ≤ 0 means Δ is acausal at every state, and it should abort, not be exported.
  `z(f, k, …)` returns `zeros(k)` when the function is absent (JuliaRoundTrip.wl). The count matches only if `JuliaParts["causality"]`
  already has the constants removed; dropping them only in ExportJulia makes the counts differ.
- Hint: drop the constants in Adimensional (the "Compact"/"Explicit" keys), and abort with a named check on a nonpositive constant.
- Resolution: _(author)_

### E7 [minor] The KPZ expectation `{1 + rn[aq^2]}` does not name which form it applies to
- Where: plan step 5, AdimensionalTest.
- Evidence: the Compact form of the KPZ condition is the atom `AdimAtom[h, hbar, 0, 0]`, whose value is `1 + rn[aq^2]` (probe output).
  Only the Explicit key equals `1 + rn[aq^2]`.
- Hint: test both keys: Compact → {AdimAtom[h,hbar,0,0]}, Explicit → {1 + rn[aq^2]}.
- Resolution: _(author)_

### E8 [minor] The name "Causality" already means something else in this codebase
- Where: Decision 7 and 8.
- Evidence: Kernel/Symmetries.wl:77 defines a class "Causality" (the Γ_φφ flows). Kernel/Closure.wl:38 has `Closure.causality`,
  Kernel/Vertices.wl:132 has `Vertices.causality`, and Kernel/Oracle.wl:65 selects `#Class == "Causality"`. The new stage key
  "Causality" and the INFO `JuliaExport.causality` mean causal poles.
- Hint: name it "CausalPoles"/`causal_poles`, or define the two meanings in CONVENTIONS.
- Resolution: _(author)_

### E9 [minor] The "minimum margin" has no meaning after stripping
- Where: Decision 5 ("a margin, not a Bool"), the handoff ("exposes the minimum margin").
- Evidence: the minimal form is defined only up to positive factors (aq², k, references, numbers), and those factors are dropped
  arbitrarily. A minimum taken over entries or over aq therefore depends on that choice.
- Hint: the contract should promise the sign only. If a margin is wanted, normalise it explicitly.
- Resolution: _(author)_

## E1-probe (rerun from the repo root, about 4 min)
```
CODE='testRoot="/home/greg/Physics/Sources/FRGscalarMSR/tests/"; Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
SetAttributes[quiet,HoldAll]; quiet[e_]:=Block[{Print},e];
adimCond[eqs_,dims_]:=Module[{n=dims["Names"],v=dims["Vectors"],b=dims["Basis"],cps=AnsatzOf[eqs]["Couplings"],r=dims["References"],tu},
 tu=Join[(#->n[#] FRGscalarMSR`Private`basisForm[v[#],b])&/@Complement[cps,r],{eqs["LoopMomentum"]->dims["Scale"] aq}];
 Factor[Together[First/@(eqs["ValidityConditions"]["Explicit"]/.tu)]]];
(* LPA: ansatz/regulator/symmetries copied from tests/acceptance/ShallowWaterLPAAdimensional.wl:9-28, dims from :34 *)
Get[testRoot<>"reference/ShallowWaterSetup.wl"];
de4Dims=quiet[DeriveDimensions[swaAnsatz,"References"->{Dphik->{"phi","Mode"->"Running"},Dhk->{"h","Mode"->"Running"},nk->{"nu","Mode"->"Running","Equation"->nk ank}},"Dimensionless"->{Ak,ank,aDhk,aDphik}]];
Print[InputForm[adimCond[swaEquations,de4Dims]]]; Print[FactorTermsList[-2x-2y]]; Exit[]'; math -noprompt -run "$CODE"
```

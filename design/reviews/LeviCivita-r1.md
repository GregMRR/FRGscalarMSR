# the-enemy — /home/greg/.claude/plans/i-want-to-use-greedy-clover.md, round 1
Attacked commit: 1fbb12e  ·  Date: 2026-10-21
Verdict: 0 blocker, 3 major, 3 minor, 2 nit

Note: plan mode is active, so this report could not be written to the default review location
(`reviews/i-want-to-use-greedy-clover/r1.md`). The evidence probes are inline `math -noprompt -run '…'` commands, quoted
below, and each one runs standalone. The target is a design plan with no machine checks of its own. The regression suite
`math -script tests/runTests.wl` was not run, because it writes a temp `.wdx` file (StagesTest.wl:7) and plan mode allows
no writes. The plan references "line" numbers below are section names, because the plan file has no stable line numbering.

## Attack surface
- F1 (wedge bilinear/antisymmetric) — hand + probe — held
- F2 (εε = δδ − δδ, d = 2) — explicit 2-vectors, probe P1 — held
- F3 (∂→−iq, jac → −wedge; Γ2 untouched) — plane waves with System D, probe P1; wedge[q,−q]: OrderedQ[{q,−q}] = False so the swap rule fires, then wedge[q,q] = 0 — held
- F4 (e = 1 terms drop out of n = 0 and n = 1 projections, all m) — explicit components with generic f(p², p·q, q², ν), probe P1 — held (pointwise in q, so the Oracle at fixed qVec agrees)
- F4 premise "Γ3 linear in wedge" — depends on "at most one ε per monomial" — finding E1
- Core `wedge` spec "linear … like scal" — compared with Core.wl:136-140 (scal gets both-slot linearity from Orderless) — finding E6
- `reduceWedge2D`, `ExactlyZeroQ` "no change" (Core.wl:68-82 randomises wedge atoms; only F2-reduced input is fed in) — held
- `displayString` a^b — finding E7
- Ansatz: rank, toFactors/derive (Ansatz.wl:52-56 appends any non-"t" index), jac[A,A] = 0, L/L′ sign rule for a single ε — held
- Ansatz canonical with ≥ 2 ε (collapse EpsIndex[_,s] → EpsIndex[s]) vs "caught by momentumOrder" — probe P2 — finding E1
- Ansatz couplingsOutsideOperators (Ansatz.wl:137, `First[#]` only) — finding E4
- Ansatz `"Epsilon"` key read downstream — finding E6
- Vertices factorSymbol/(−i)^n, planeWave in 2D (D-based, independent of the parser; catches the swap sign because the toy puts Psi at leg 3 and H at leg 1), Γ2 checks — held
- Symmetries / Loops / FrequencyIntegration / Closure "no change" — held (Γ2 wedge-free; linearity collapses wedge[p+q,q])
- Projection: d = 2, LHS with the same dim, Laplacian/(2d) = angular average for p_i p_j M_ij — held; detection — E6
- Assemble "keeps Dimension -> Dim" — finding E3
- Oracle: d = 2 only, explicit wedge — E8 (nit); it shares Gamma3Matrix and the wedge rules with the symbolic chain (independent of F2/F4/Projection only, as the plan says) — held
- Unit-test hand value −c wedge[q_H, q_Psi] for c Hbar jac[H,Psi] — VertexValue leg convention: −c wedge[q2, −q1−q2] = −c wedge[q1,q2], consistent with F3 — held
- Acceptance toy model — findings E2, E5
- Regression constraint (byte-identical ε-free path) — sign = 1 path, GroupBy keys, card text — held, except E6 (Missing key)
- Docs list — finding E3 (FRGDESIGN.md decision not updated)

## Findings
### E1 [major] The canonical form collapses ε labels, so a nonzero jac×jac (p⁴) monomial is declared zero and dropped before momentumOrder can catch it
- Where: plan §Ansatz "canonical" (map EpsIndex[_, s] → EpsIndex[s], then L′ === L ⇒ zero) and §Ansatz "Checks" ("jac × jac is p⁴ and is caught by momentumOrder anyway"); Ansatz.wl:150-153 (canonical/drop) run before Ansatz.wl:162 (momentumOrder) and :164 (indexPairs).
- Evidence (probe P2):
  `math -noprompt -run 'fac={FF[Hbar,0,{}],FF[A,0,{EpsIndex[j1,1]}],FF[B,0,{EpsIndex[j1,2]}],FF[B,0,{EpsIndex[j2,1]}],FF[A,0,{EpsIndex[j2,2]}]}; c=fac/.EpsIndex[_,s_]:>EpsIndex[s]; Print[Sort[c]===Sort[c/.{EpsIndex[1]->EpsIndex[2],EpsIndex[2]->EpsIndex[1]}]]; w[a_,b_]:=a[[1]]b[[2]]-a[[2]]b[[1]]; Print[Factor[w[{a1,a2},{b1,b2}] w[{b1,b2},{a1,a2}]]]; Exit[]'`
  → `True` (declared zero), while the true value is `-(a2*b1 - a1*b2)^2` ≠ 0. So `c Hbar jac[A,B] jac[B,A]` gets only an "Ansatz.cancelled" INFO line, and neither momentumOrder nor indexPairs sees it.
- Consequence: an unsupported p⁴ term disappears without an abort. The pipeline then returns equations for an ansatz that silently lacks that term.
- Hint: keep the ε identity (j) through the zero/sign test, or run momentumOrder/indexPairs on the raw factor lists before `canonical`.
- Resolution: _(author: "fixed in <commit>: …" or "rebutted: …")_

### E2 [major] The acceptance toy never generates an e = 1 term, so F4 and Projection.wedgeOdd are never tested end to end
- Where: plan §Tests "Acceptance (tests/acceptance/Jacobian2D.wl)"; §Projection "Projection.wedgeOdd … evidence … must be exactly 0".
- Evidence: in `Hbar(dt[H] + λ jac[Psi,H] − ν lap[H]) − D Hbar² + Psibar(dt[Psi] − μ lap[Psi]) − E Psibar²`, the only cubic monomial is the jac term and there is no quartic, so every Γ3 entry is pure wedge and every bubble term has e = 2 exactly. The dropped part is therefore identically 0, and the wedgeOdd "evidence" compares 0 with 0. Mutation: replacing the F4 drop with "keep e = 1 terms" or "drop e = 0 terms with a wedge-free coefficient" would not change the toy's output. (A non-ε cubic like KPZ's grad[H].grad[H] would not help either: without H–Psi mixing in G it never meets the jac vertex in a bubble.) Each Γ3 is also ∝ wedge[p,q], so every n = 0 flow vanishes identically. Only the p² flows (ν, and a generated p² noise term) are nonzero, which makes the Oracle comparison very narrow.
- Consequence: the one non-trivial manipulation that has no independent check elsewhere (dropping e = 1 terms) ships untested.
- Hint: add a non-ε cubic term with the same field content as the jac term, e.g. `κ Hbar grad[Psi].grad[H]`, so that Γ3_{Hbar Psi H} mixes scal and wedge and the bubble has e = 1 cross terms.
- Resolution: _

### E3 [major] The ε path silently fixes d = 2, but the output still advertises a symbolic dimension (convention and decision drift)
- Where: plan §Theory ("Assemble keeps "Dimension" -> Dim, but the equations become free of Dim"); FRGDESIGN.md:73 ("d is a symbol (Dim) throughout"), :184 (Assemble contract "Dim symbolic"), :186 (E12: Oracle in d = 2 and d = 3), Q16 (Equations is the input of Adimensional/JuliaExport); CONVENTIONS.md "The space dimension d is a user symbol … throughout"; Assemble.wl:128 card "Dim symbolic", :133.
- Evidence: Projection substitutes dim = 2, so the projected integrands are valid only at d = 2. The Equations stage still carries `"Dimension" -> Dim`, and the q-measure ∫dᵈq/(2π)ᵈ is left to downstream numerics (Loops.wl:6), which read that key. The claim "equations become free of Dim" is also false whenever a coefficient contains Dim, which Ansatz.wl:129-132 allows. The Docs list updates CONVENTIONS/MODULES but not the FRGDESIGN decisions it overrides.
- Consequence: any consumer that integrates with the stage's Dim (the Q16 FlowRG export, or a user doing d = 3) gets wrong flows without any check firing.
- Hint: make the Equations (and Symmetries) stage carry the fixed dimension 2 (or a "FixedDimension" key that consumers assert), change the Assemble card text, and record the exception in FRGDESIGN.md and CONVENTIONS.md.
- Resolution: _

### E4 [minor] Only jac's first argument goes through the coupling-outside-operator check
- Where: Ansatz.wl:137 (`FreeQ[First[#], …]`); plan §Ansatz updates operatorArity and operatorOnField but not couplingsOutsideOperators.
- Evidence: `jac[H, nk Psi]` passes the check, against CONVENTIONS.md "Couplings … stand outside the operators" and MODULES.md:77 "never sit under an operator". The result is numerically harmless (Ansatz.wl:56 gives derive[c] = 0), but it breaks the stated DSL contract.
- Hint: test every argument (`Level[#, {1}]`) for the jac head.
- Resolution: _

### E5 [minor] The acceptance toy, as written, cannot run, and its criterion is underspecified
- Where: plan §Tests "Acceptance".
- Evidence: `math -noprompt -run 'Print[Context[D]," ",Context[E]," ",MemberQ[{"k","q","p","y","t","x","\[Nu]","\[Omega]"},SymbolName[\[Nu]]]];Exit[]'` → `System` `System` `True`. `D` is also in systemOps (Ansatz.wl:26), so Ansatz.noSystemDerivatives, assertFreshSymbols .system, and .reserved for ν all abort. The plan also gives no DeclareSymmetries tags or TruncationOrders and no regulator. The generated p² term in Γ2_{Hbar Hbar} is outside the ansatz and needs a Neglected/Protected tag, or Symmetries.allTagged aborts.
- Hint: rename the couplings (e.g. nuk, Dk, Ek, muk, lamk) and specify the tags, truncation orders and regulator, as KPZSetup.wl:7-12 does.
- Resolution: _

### E6 [minor] Several spec items would break the ε path or the regression constraint if implemented as literally written
- Where: plan §Core `wedge` ("linear over sums and explicit numbers, like scal"); §Projection "Detect ε via loops["Vertices"]["Ansatz"]["Epsilon"]"; §Oracle ("Dimensions" defaults to {2}, an explicit request with 3 aborts).
- Evidence: (a) scal gets linearity in both slots only from Orderless (Core.wl:136-140). `math -noprompt -run 'Print[{OrderedQ[{p,p+qL}],OrderedQ[{p,-p-qL}]}];Exit[]'` → `{True, True}`, so wedge[p, −p−qL], which Gamma3Matrix produces for the last leg, is "ordered": the swap rule does not fire, and with first-slot rules alone it never reduces. wedgeOdd then aborts. `wedge[0, b] = 0` is also unspecified (scal has it). (b) `If[<|"Stage"->"Ansatz"|>["Epsilon"], 1, 0]` stays unevaluated as `If[Missing["KeyAbsent","Epsilon"], 1, 0]` (probe output), which breaks the ε-free path. (c) `OptionValue` cannot distinguish the default {2,3} from an explicit {2,3}.
- Hint: specify sum and number rules in both slots plus `wedge[0,_] = 0`, use `TrueQ[Lookup[…, "Epsilon", False]]`, and detect an explicit "Dimensions" with `FilterRules[{opts}, "Dimensions"]`.
- Resolution: _

### E7 [nit] Displaying wedge[a,b] as `a^b` is read as a power
- Where: plan §Core displayString; Core.wl:47 already prints scal[a,a] as `a^2`.
- Evidence: `wedge[q1, q2]` would print `q1^q2`, which looks like an exponent, in the vertex table the user checks by hand (Q2).
- Hint: print `a×b` or `wedge(a,b)`.

### E8 [nit] Oracle `dot` with a₁b₂ − a₂b₁ on the symbolic pVecO emits Part::partd
- Where: plan §Oracle; Oracle.wl:63-72 (pVecO is substituted only after `dot`).
- Evidence: `cross[x_,y_]:=x[[1]]y[[2]]-x[[2]]y[[1]]; cross[pV,{3/10,-1/2}] /. pV->{1,2}` gives the right value (−11/10), but prints Part::partd twice.
- Hint: keep the wedge inert (`crossVec[a_List, b_List]`) until pVecO is a list.

## Probe P1 (F2, F3, F4 — all held)
`math -noprompt -run 'w[a_,b_]:=a[[1]]b[[2]]-a[[2]]b[[1]]; {a,b,c,d}=Table[{x[i,1],x[i,2]},{i,4}]; Print[Expand[w[a,b]w[c,d]-((a.c)(b.d)-(a.d)(b.c))]]; X={x1,x2};kA={ka1,ka2};kB={kb1,kb2};A=Exp[-I kA.X];B=Exp[-I kB.X]; Print[Simplify[(D[A,x1]D[B,x2]-D[A,x2]D[B,x1])/(A B)+w[kA,kB]]]; P={p1,p2};Q={q1,q2}; g=w[P,Q] f[P.P,P.Q,Q.Q,nu]; Print[D[g,{nu,2}]/.{p1->0,p2->0}]; Print[Simplify[(D[g,{p1,2}]+D[g,{p2,2}])/.{p1->0,p2->0}]]; Print[Simplify[D[D[g,{p1,2}]+D[g,{p2,2}],nu]/.{p1->0,p2->0}]]; Exit[]'`
→ 0, 0, 0, 0, 0.

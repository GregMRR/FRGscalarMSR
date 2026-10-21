# the-enemy — plans/i-want-to-add-kind-wreath.md (momentum projection at any order), round 1
Attacked commit: 6208fe0 (dirty tree: uncommitted Kernel/Adimensional, Dimensions, JuliaExport, …)  ·  Date: 2026-10-21
Verdict: 1 blocker, 3 major, 6 minor, 2 nit

Intended location: design/reviews/MomentumProjection-r1.md (project convention, cf. design/reviews/LeviCivita-r1.md).
Plan mode was active, so this file is the only one written. The probes are inline below. Run each with
`math -noprompt -run '<code>'` from the repo root.
Step 2 (confirm green): the target is a design with no checks of its own. The full suite was not run because it writes temp
files, which plan mode forbids. The KPZ pipeline up to Loops builds cleanly in every probe.

## Attack surface
- Plan F1 ⟨c^{2k}⟩ table — explicit (∇²)^k e₁^{2k}/N_k, k = 1..3, d = 1..5 (P1) — held
- Plan F2 N_n = 4^n n! (d/2)_n — explicit (∇²)^n (p²)^n, n = 1..3, d = 1..4 (P1) — held
- Plan unit-test closed form p⁴ — explicit (∇²)²F((p+q)²)/N₂ in d = 1..4, which pins all 3 d-structures (P1) — held
- "Analogous p⁶ case" — derived f‴/6 + q²f⁗/d + 2q⁴f⁽⁵⁾/[d(d+2)] + 4q⁶f⁽⁶⁾/[3d(d+2)(d+4)], checked d = 1..4 (P1) — held (formula not in plan: E12)
- F3 (odd-ω drop exact for all n) — reflection of e about q keeps c and flips ω, and the circle measure is reflection invariant — held
- F1 applied pointwise in q to a non-isotropic integrand — equals F2 pointwise, since (∇²_p)^n commutes with joint rotations — held (wording: E12)
- F4 / downstream generality — deltaCompletion general in j (Projection.wl:23-24), Assemble jmax general (Assemble.wl:84), FrequencyIntegral any pole order (FrequencyIntegration.wl:24-46) — held
- F5 — code trace: Ansatz.wl:90-91 relabels by Union (creation-counter order). lap[lap[A]] → {S1,S1,S2,S2} and div[grad[lap[A]]] → {S2,S2,S1,S1} — held. Consequence is cosmetic: vertices are sums (Vertices.wl:30-39)
- Invariant (n ≤ 1 identical after Expand) — old `project` vs my F1 implementation on all 16 KPZ monomials (P2) — held (Expand-identical and ExactlyZeroQ-equal)
- Closure normalisation, item 6 — project returns the coefficient before and after the change — held. Doc text: E12
- Symmetries declOutside, item 3 — Symmetries.wl:94 fixed n ≤ 2, and MomentumOrder/2 bounds the Γ2 p-power — held
- Cited line numbers (Projection.wl:45-57, Ansatz.wl:89/196, Symmetries.wl:60/93, Oracle.wl:86-90, JuliaExport.wl:22, Closure.wl:52, MODULES.md:460) — held
- Background claim "iterating gives 24d²" (plan:5, FRGDESIGN.md:45) — iterated dimFactors (P5) — finding E10
- Projection.explicitVectors (exact equality at a random rational point) — mutation of ⟨c⁴⟩ (P3) — finding E1, E4
- Projection.LHS "now any n" — what it exercises — finding E2
- Oracle at p^{2n} — stencil error vs NIntegrate PrecisionGoal (P4) — finding E3
- Projection.wedgeOdd for n ≥ 2 — code read, Projection.wl:38 — finding E5
- evenInQ replacement — degree counting — finding E6
- JuliaExport D{k}f — invariant — finding E7
- Regulator smoothness at order 2n — finding E8
- canonical "sort then min over permutations" — finding E9
- Notation and docs list — findings E11, E12
- explicitVectors timing — KPZ p⁴: 0.03 / 0.1 / 0.3 s per integrand at d = 1 / 2 / 3 (P3) — held (cheap)

## Findings
### E1 [blocker] The exact random-rational zero test is blind to derivatives of GammaAtom / DeltaCoefficient atoms, so explicitVectors as planned cannot fail on a wrong p⁴ projection
- Where: plan item 1 `Projection.explicitVectors` ("random rational point … must equal the symbolic result exactly; generalises
  … `Projection.wedgeOdd`"). Earlier-step defect in Core.wl:85-86: `atomRule` matches `(h_Symbol)[___]` *inside*
  `Derivative[k][GammaAtom[a,b,m]]`. It replaces `GammaAtom[a,b,m]` by a number, and `Derivative[k][number]` evaluates to 0.
  The same happens to `Derivative[j][DeltaCoefficient[m]]`. Core.wl:28-30 (usage) claims that r′[y] atoms get random values.
  The same blindness affects the existing evidence at Projection.wl:40.
- Evidence (P3):
  `ExactlyZeroQ[Derivative[1][GammaAtom[hb,h,0]][QQ]]` → True. The undifferentiated `ExactlyZeroQ[GammaAtom[hb,h,0][QQ]]` instead
  aborts with `(9/86)[2/3]` non-numeric. I mutated ⟨c⁴⟩ to 3q⁴/(7d−6) and compared the KPZ p⁴ integrand of Γ_{h h̄} with the
  explicit (∇²)²/N₂. ExactlyZeroQ returns True for the mutation at d = 1, 2 and 3. A sound test, with every atom replaced by a
  fresh symbol, rejects the mutation at d = 1 (output below).
- Consequence: the new runtime check of the angular average compares only the derivative-free terms. Those terms never carry
  c^{2k} with k ≥ 1 from the propagator. So a wrong table, or a wrong s-order, passes.
- Hint: before any exact comparison, map every `f[..][x]`, `Derivative[..][f[..]][x]` and `Δ[..]` to its own fresh symbol, or
  fix atomRule to match whole atoms first.
- Resolution: _(author)_

### E2 [major] Projection.LHS no longer tests any d-dependent momentum normalisation
- Where: plan "Consistency checks", `Projection.LHS … now any n`. MODULES.md:386-387 says it "checks the normalisation of both
  projections, including 2d".
- Evidence: Γ2 entries depend on p only through scal[p,p], and Vertices.wl:135 asserts this. F1 maps (p²)^n → s^{2n}, and
  D^{2n}/(2n)! gives 1 with no c present. ⟨c^{2k}⟩ is never evaluated, so any table, including the E1 mutation, passes LHS.
  The old operator did test 2d, because ∇²(p²) = 2d.
- Consequence: the runtime guard for the d-normalisation moves entirely to explicitVectors (blind, E1) and the Oracle (E3, E4).
- Hint: state this in the card and docs, and rely on a sound explicitVectors (E1, E4).
- Resolution: _(author)_

### E3 [major] At p⁴ the Oracle cannot meet its 10⁻⁸ tolerance by raising the working precision: PrecisionGoal 22 is the floor
- Where: plan item 4 ("working precision raised with n so 2n-th differences at h = 10⁻⁵ keep enough digits"), Oracle.wl:84
  (`PrecisionGoal -> 22`, fixed).
- Evidence (P4): F = 1/(1+x), q = (3/10, −1/2), d = 2, nested 2nd-difference stencils, h = 10⁻⁵:
  - f relative error 10⁻²² gives a p⁴-coefficient relative error of 1.3·10⁻²;
  - f relative error 10⁻³² gives 1.4·10⁻⁹;
  - truncation alone gives 1.4·10⁻⁹, which is only 7× below the tolerance.
- Consequence: the Oracle fails for every p⁴ equation, and the tempting fix is to loosen the tolerance.
- Hint: scale PrecisionGoal (and WorkingPrecision) to roughly 8 + 10n + margin, i.e. ≥ 32 for n = 2 and ≥ 42 for n = 3.
- Resolution: _(author)_

### E4 [minor] Checking only d = 2 and 3 cannot pin the d-dependence for n ≥ 2
- Where: plan explicitVectors "Dim → 2 and 3", Oracle d = 2, 3, FUTURE.md:16 and FRGDESIGN.md:196 ("validate with the oracle in
  d = 2 and d = 3").
- Evidence: the p^{2n} result lies in span{1, 1/d, …, 1/[d(d+2)⋯(d+2n−2)]}, which has dimension n+1. The mutation 3/(7d−6)
  equals 3/(d(d+2)) at d = 2 and 3. Even with the sound atom test of P3 it passes at d = 2 and 3 and fails only at d = 1. The
  planned symbolic-d unit closed form does cover the table for the f((p+q)²) case.
- Consequence: the runtime checks certify general d only at n ≤ 1.
- Hint: add d = 1 (0.03 s/integrand) and, for n = 3, d = 4; in general use n+1 dimensions.
- Resolution: _(author)_

### E5 [major] Projection.wedgeOdd evidence is identically 0 = 0 for n ≥ 2, yet the plan lists it as the runtime check of F3
- Where: Projection.wl:38 `If[n == 1, ex = Laplacian]`. For n ≥ 2 no p-derivative is taken, so ex = (pc₁qc₂ − pc₂qc₁)·odd at pc = 0,
  which is 0 for any `odd`. The plan's item 1 only says dropWedge "keeps" the reduction and the drop.
- Consequence: F3 at n ≥ 2 has a runtime check only if explicitVectors starts from the pre-drop integrand (which the plan does
  not say) and is sound (E1). explicitVectors at d = 3 is also impossible with wedge, which the plan does not declare.
- Hint: generalise line 38 to (∇²)^n/N_n, or state that explicitVectors (d = 2, raw integrand) replaces it.
- Resolution: _(author)_

### E6 [minor] The proposed evenInQ replacement is tautological, and ⟨c^odd⟩ = 0 would mask a bug
- Where: plan item 1, "a check that only even powers of c survive before averaging, after odd ones are discarded".
- Evidence: after the substitution, s enters only as s² and s·c (ω is dropped first). So the s^{2n} coefficient has c-degree
  2n − 2a, which is always even. An odd power can only come from a substitution bug, and discarding it hides that bug.
- Hint: assert that no odd c power is present, and do not discard anything.
- Resolution: _(author)_

### E7 [minor] The "D{k}f" naming would change existing exports, against the Invariant
- Where: plan item 5 and JuliaExport.wl:24, where Derivative[1] → `D`, Derivative[2] → `D2`. Literal `D1f` breaks the round-trip
  and export invariant. Because k ≥ 3 aborts anyway, the generalisation adds nothing.
- Hint: keep `D`/`D2` and only replace the assert message.
- Resolution: _(author)_

### E8 [minor] Undeclared regularity: a p^{2n} projection needs shapes in C^{2n}
- Where: plan F1/F4. Litim is used by the Canet comparison (MODULES.md:715). FUTURE.md:37 already notes the δ in r″. At p⁴,
  r‴ and r⁗ bring δ′ and δ″, which the pointwise pipeline cannot represent.
- Hint: declare the assumption in MOMENTUMPROJECTION.md and CONVENTIONS.
- Resolution: _(author)_

### E9 [minor] The canonical-form recipe is ambiguous
- Where: plan item 2, "sort each factor's index list, then min over all permutations of dummy labels".
- Evidence: the sort must be redone after each relabelling, because the order depends on the labels. The ε zero rule (the list
  equals its slot swap modulo relabelling) must also be restated for ≥ 2 dummies.
- Resolution: _(author)_

### E10 [minor] The quoted "24d²" for the iterated Laplacian trick is not reproduced
- Where: plan:5, FRGDESIGN.md:45 (E2).
- Evidence (P5): iterating `dimFactors` on scal[p,p]² gives 8(d²+2), not 24d². Both differ from 8d(d+2).
- Hint: fix the number, or state the iteration it refers to.
- Resolution: _(author)_

### E11 [nit] The F-labels collide with LEVICIVITA's F1–F5
- Where: plan item 1. "dropWedge keeps the F2 reduction" means LEVICIVITA F2, while the plan's own F2 is the equivalence. The
  new MOMENTUMPROJECTION.md will reuse F1–F5.
- Resolution: _(author)_

### E12 [nit] The docs and change list miss stale texts and the p⁶ formula
- MODULES.md:460-461 (Closure: "(2d)ⁿ … in Value") should read N_n.
- Ansatz.wl:178 hint "(p^4 is unsupported)" is stale.
- F1 "of an isotropic f" should say "pointwise in q, by rotation invariance of (∇²_p)^n".
- The plan should give the p⁶ closed form (above).
- F5 says "encounter order"; the code uses creation-counter order.
- Resolution: _(author)_

## Probes (inline)
P1 (F1/F2/closed forms):
```
lapn[e_,ps_,n_]:=Nest[Function[g,Sum[D[g,{x,2}],{x,ps}]],e,n]; Nn[n_,d_]:=4^n n! Pochhammer[d/2,n];
closed[2,d_,Q_]:=F2/2+2Q F3/d+2Q^2 F4/(d(d+2)); closed[3,d_,Q_]:=F3/6+Q F4/d+2Q^2 F5/(d(d+2))+4Q^3 F6/(3d(d+2)(d+4));
Table[Module[{ps=Array[pp,d],qs=Array[qq,d],ex},ex=lapn[F[(ps+qs).(ps+qs)],ps,n]/.pp[_]->0/.Derivative[k_][F][_]:>ToExpression["F"<>ToString[k]];
 {n,d,Expand[ex/Nn[n,d]-(closed[n,d,Q]/.Q->qs.qs)]===0,Nn[n,d]===lapn[(ps.ps)^n,ps,n]}],{n,2,3},{d,1,4}]
```
→ all True; the ⟨c^{2k}⟩ table minus lapn[e₁^{2k}]/N_k is {0} for k = 1..3, d = 1..5.

P2 (n ≤ 1 invariant). Build KPZ (tests/reference/KPZSetup.wl) up to Loops, then compare `FRGscalarMSR`Private`project` with
`newProject`:
```
newProject[e0_,m_,n_,{nu_,p_},{w_,q_},dim_,deg_,tab_:avg]:=Module[{e=e0,cl},If[m>0,e=D[e,{nu,m}]];
 e=e/.{scal[p,p]->ss^2,scal[p,q]->ss cc}; e=(D[e,{ss,2n}]/.ss->0/.nu->0)/(2n)!; cl=CoefficientList[Expand[e],cc];
 e=Sum[cl[[k+1]] tab[k,dim,scal[q,q]],{k,0,Length[cl]-1}];
 e=e/.scal[q,q]->LoopMomentumSquared/.FRGscalarMSR`Private`deltaCompletion[deg]; e/(I^m m!)];
avg[k_,dim_,Q_]:=If[OddQ[k],0,Q^(k/2)(k-1)!!/Product[dim+2j,{j,0,k/2-1}]];
```
→ 16/16 monomials Expand-identical.

P3 (E1/E4). Same setup. Use `mut[k_,dim_,Q_]:=If[k==4,3Q^2/(7dim-6),avg[k,dim,Q]]`, and `explicit` = components substituted
into scal, (∇²)² in the p components, p → 0, deltaCompletion, divided by N₂, at q = (3/7, −2/5, 5/11)[[;;d]]. Define
`fresh[e]` to replace every `Derivative[__][_[__]][__] | (_[__])[__] | Δ[__]` by a Unique symbol. Output:
```
d=1: ExactlyZeroQ correct True, mutated True | fresh: correct True, mutated False
d=2: ExactlyZeroQ correct True, mutated True | fresh: correct True, mutated True
d=3: ExactlyZeroQ correct True, mutated True | fresh: correct True, mutated True
ExactlyZeroQ[Derivative[1][GammaAtom[hb,h,0]][QQ]] -> True
```

P4 (E3):
```
F[x_]:=1/(1+x);q={3/10,-1/2};d=2;Q=q.q;exact=F''[Q]/2+2Q F'''[Q]/d+2Q^2 F''''[Q]/(d(d+2));
st2[h_]:={{h,1/h^2},{0,-2/h^2},{-h,1/h^2}};u=IdentityMatrix[d];
run[noise_,h_]:=Module[{f},SeedRandom[3];f[p_]:=N[F[(p+q).(p+q)],80](1+noise RandomReal[{-1,1},WorkingPrecision->80]);
 N[Abs[(Total[Flatten[Table[a[[2]]b[[2]]f[a[[1]]u[[c1]]+b[[1]]u[[c2]]],{c1,d},{c2,d},{a,st2[h]},{b,st2[h]}]]]/(8d(d+2))-exact)/exact],5]];
{run[10^-22,10^-5],run[10^-32,10^-5],run[0,10^-5]}
```
→ {1.3·10⁻², 1.4·10⁻⁹, 1.4·10⁻⁹}

P5 (E10):
```
e1=D[scal[pp,pp]^2,pp,pp]/.FRGscalarMSR`Private`dimFactors[dd,qq,lam]; Factor[D[e1,pp,pp]/.FRGscalarMSR`Private`dimFactors[dd,qq,lam]/.pp->0]
```
→ 8(2 + dd²)

# the-enemy — plans/i-want-to-add-kind-wreath.md (momentum projection at any order), round 2
Attacked commit: 6208fe0 (dirty tree)  ·  Date: 2026-10-21
Verdict: 0 blocker, 4 major, 3 minor, 3 nit

Where this report lives: it belongs at `~/.claude/plans/MomentumProjection-r2.md`, beside r1, and later at
`design/reviews/MomentumProjection-r2.md`. Plan mode allowed only this file to be written. The probes are inline below. Each one
runs standalone with `math -noprompt <<'EOF' … EOF` from the repo root and writes no files.
Step 2 (confirm green): the target is a design and has no checks of its own. The probes instead build the KPZ, jac-2D and
shallow-water pipelines. With the step-0 zero test patched in, all three build without any abort (P1).

## Re-attack of round-1 resolutions
- E1 (blind zero test) → step 0. **Held.** The patched rule gives False on `GammaAtom′`, `DeltaCoefficient″` and `f′ − f`, and
  True on the chain-rule identity. The KPZ (+ Oracle), jac-2D and shallow-water pipelines build with no failing check (P1). On
  the KPZ Γ_{h̄h} loop at p⁴, the mutated ⟨c⁴⟩ = 3q⁴/(7d−6) is now rejected at d = 1 (P4). Results are unchanged: the only
  result-changing use, Assemble.wl:106, sees `ScaleAtom[..]`/`dsR[..]`, which were never blind.
- E2 → `Projection.normalisation`. Held in principle. Its reference is underspecified (E17).
- E3 → PrecisionGoal 8 + 10n + margin. Held for a single n. The fix breaks when records of different n share the memo (E15).
- E4 → explicitVectors at Dim = 1…n+1. Held for symbolic Dim (P4). It breaks for a declared d without ε (E13).
- E5 → wedgeOdd with explicit (∇²)ⁿ. Held. P2 shows jac 36/36 Expand-identical and explicitVectors d = 2 passing on the
  pre-drop integrand.
- E6 → noOddC that asserts and does not discard. Held. The check cannot be falsified, and the plan has no test for it (E18).
- E7 → keep `D`/`D2`. Held. The wording contradicts itself (E20).
- E8 → `Projection.smoothShape`. **Failed.** The check can never fire (E16).
- E9 → the canonical recipe. Held: the min over relabellings of (sort inside factors, sort factor list) is a class
  invariant, and the ε rule is right. A parser defect upstream of it is not covered (E14).
- E10, E11, E12 → fixed in the text. Held.

## Attack surface
- Invariant n ≤ 1, new M1 operator vs the old Laplacian (plan:11-15) — Expand-identity on every integrand: KPZ 16/16, jac-2D
  (ε, d = 2) 36/36, shallow-water DE4 56/56 (P2, P3) — held
- Step 0, sound ExactlyZeroQ (plan:31-35) — unit cases plus three full pipelines (P1) — held. Duplicate rule elsewhere: E21
- M1/M2/p⁴/p⁶ closed forms (plan:19-22, 68-70) — r1 P1 re-used; the Γ-function form of ⟨c^{2k}⟩ checked by hand — held
- M3 (odd-ω drop for every n) — reflection argument; jac explicitVectors d = 2 on the raw integrand (P2) — held
- explicitVectors at Dim = 1…n+1 (plan:41-43) — mutation (P4) — held. Declared d = 2 without ε (P5) — finding E13. Cost (P3) — E22
- Ansatz lift and canonical (plan:48-52, 75-76) — code trace of toFactors (P6) — finding E14
- Oracle at n ≥ 2 (plan:60-61) — memo across n (P7) — finding E15. Truncation at h = 10⁻⁵ (P8) — finding E19
- smoothShape (plan:55-59) — trace of Regulator.wl:74-77 and KPZ.wl:13-16 — finding E16
- Projection.normalisation (plan:45-47) — independence of its reference — finding E17
- New checks vs the convention "new check + expectFail" (CLAUDE.md "Where to change X") — finding E18
- JuliaExport (plan:62-63) — finding E20
- Downstream generality (Vertices.monomialAt, gammaCoefficient, Closure, Assemble jmax, Dimensions MomentumOrder) — code read — held

## Findings
### E13 [major] explicitVectors at Dim = 1…n+1 aborts every declared-d pipeline without ε (ShallowWaterDeclared2D)
- Where: plan:41-43 ("at Dim = 1, …, n+1 … with ε only d = 2"); Projection.wl:65 (`dim = symmetries["Dimension"]` is the
  number 2 under FixedDimension[2]); tests/acceptance/ShallowWaterDeclared2D.wl.
- Evidence (P5): KPZ with `FixedDimension[2]` and no jac. Every p² integrand fails the explicit comparison at d = 1, and all
  of them pass at d = 2:
  `declared d=2 KPZ, explicitVectors at d=1: failing {{{h,h},0},{{h,h},1},…,{{hbar,hbar},2}}` / `at d=2: failing {}`.
- Consequence: as written, the Invariant breaks: the existing declared-d = 2 acceptance test aborts at n = 1.
- Hint: when dim is numeric (declared), test at that d only; use 1…n+1 only for a symbolic Dim.
- Resolution: _(author)_

### E14 [major] Lifting the p² cap exposes a parser defect: a power of a contracted scalar reuses its dummy (lap[h]², (∇h·∇h)²)
- Where: Ansatz.wl:68 `S[Power[x_, n_]] := Expand[S[x]^n]` (an earlier step). The plan's item 2 (plan:48-52) and its tests
  "(grad A.grad A)^2 … → one monomial" and "p⁴ monomials accepted" (plan:75-76) do not change it.
- Evidence (P6): `(grad h.grad h)^2 → FieldFactor[h,0,{S[{dummy,1}]}]^4` puts one dummy on four slots.
  `lap[h] lap[h]` is auto-squared by Mathematica and gives the same result. `indexPairsQ → False`, so the case aborts at the
  internal `Ansatz.indexPairs`.
- Consequence: the natural p⁴ terms (lap[h̄]², (∇h)⁴, (∇²h)²) cannot be parsed. The planned AnsatzTest fails, and
  canonical's min over relabellings never gets to see two dummies.
- Hint: give each power factor its own dummies: `S[Power[x_, n_]] := Expand[Times @@ Table[S[x], {n}]]`. Add the change to
  item 2 and add a test for it.
- Resolution: _(author)_

### E15 [major] The Oracle memo mixes precisions across n: an n = 2 stencil reuses the p = 0 and ±h e_c values computed for n ≤ 1
- Where: plan:60-61 ("PrecisionGoal ≈ 8 + 10n + margin … n ≤ 1 values unchanged") together with Oracle.wl:76, whose memo key
  is `{pair, nu, pVec}` and carries no precision. The nested n = 2 stencil contains every n ≤ 1 stencil point.
- Evidence (P7): F = 1/(1+x), d = 2, h = 10⁻⁵. With every point at relative error 10⁻³⁴ the error is 1.4·10⁻⁹. With the
  n ≤ 1 points at 10⁻²² (memoised at PG 22) and the others at 10⁻³⁴ the error is **1.7·10⁻²**.
- Consequence: every p⁴ oracle comparison fails, which is the trap r1 E3 warned about (loosening the tolerance).
- Hint: put the precision (or n) in the memo key, or run every record at the precision of the largest n and drop the promise
  that n ≤ 1 values stay unchanged.
- Resolution: _(author)_

### E16 [major] (resolution of r1 E8 fails) `Projection.smoothShape` can never fire
- Where: plan:55-59. Regulator.wl:74-77 (`Regulator.onlyDeclaredShapes`, `Regulator.shapeArgument`) forces every shape to be a
  declared symbol r[q²/k²], so no UnitStep or Piecewise can reach the projection. Explicit shapes appear only after the
  library has finished: KPZ.wl:13-16 and KPZAdimensional.wl:13 substitute Litim as `1/y − 1` with no θ and integrate over y < 1
  by hand. That form would pass the planned pattern test anyway.
- Consequence: the claimed C^{2n} guard is vacuous. At p⁴, r‴ and r⁗ (δ and δ′ at y = 1 for Litim) are lost silently,
  downstream, wherever an explicit shape is substituted. The plan also gives no expectFail that shows the check can fire.
- Hint: make it a card ASSUMPTION line plus a doc entry (MOMENTUMPROJECTION/CONVENTIONS), and put any guard where shapes become
  explicit (test helpers, FUTURE FlowRG item). Otherwise drop the abort.
- Resolution: _(author)_

### E17 [minor] The reference of `Projection.normalisation` is not independent of `project` for general n
- Where: plan:45-47 ("closed forms (M1/M2)"). For (p·q)^{2k}(p²)^{n−k} the "closed form" is the M1 table itself. The
  F((p+q)²) closed form is given only for n = 2, 3 (plan:68-69). For other n it can only be built from the same table.
- Evidence: a check whose expected value comes from the same `directionalAverage` table cannot fail on a wrong table. The
  plan's mutation test (plan:71) exposes this only if the mutation goes into project's copy alone.
- Hint: write the reference independently, e.g. ⟨c^{2k}⟩ = q^{2k} Γ(k+½)Γ(d/2)/(√π Γ(k+d/2)) (equal to (2k−1)!!/[d(d+2)⋯],
  checked by hand), and say where the mutation is injected.
- Resolution: _(author)_

### E18 [minor] New checks without expectFail tests; noOddC cannot be falsified
- Where: plan:38-47, 55-58, 67-78. The convention is "New check: … + an `expectFail[…, "Module.check"]` test" (CLAUDE.md).
  noOddC, smoothShape and wedgeOdd at n ≥ 2 have none. normalisation and explicitVectors are covered only by the mutation.
- Evidence: after the substitution, s enters only as s² and s·c. The s^{2n} coefficient therefore has c-degree 2n − 2a, which
  is always even. The loop integrand holds only the vectors p and q, so no leftover scal[p, X ≠ q] can exist either.
- Hint: add expectFail tests built from synthetic inputs, and label noOddC an internal invariant.
- Resolution: _(author)_

### E19 [minor] The Oracle step h = 10⁻⁵ leaves only a factor ~2.6 of headroom at n = 2, 3, and the Oracle stays at d ∈ {2, 3}
- Where: plan:60-61 (step size not addressed); Oracle.wl:44 (h); Oracle.wl:47 (`SubsetQ[{2, 3}, dims]`); plan:91 (verification
  "Oracle d = 2, 3").
- Evidence (P8): exact functions on the oracle shapes (k = 10/17) give truncation errors of 0.6–3.8·10⁻⁹ against the 10⁻⁸
  tolerance, for n = 2, 3 and d = 2, 3. A coefficient with internal cancellations goes over the tolerance. Two dimensions
  cannot pin the 3 (or more) d-structures at n ≥ 2. explicitVectors covers this case, but the Oracle docs still call the
  Oracle "the" independent check.
- Hint: use a higher-order stencil or Richardson extrapolation, with PG raised to match. Say in the Oracle card that the
  d-dependence at n ≥ 2 is certified by explicitVectors.
- Resolution: _(author)_

### E20 [nit] JuliaExport item contradicts itself
- Where: plan:62-63: "add `D3`, `D4`, … for k ≥ 3, and abort". Names that are always aborted on are never emitted.
- Hint: "keep the abort for k ≥ 3 and change its message".
- Resolution: _(author)_

### E21 [nit] Adimensional.wl:64 keeps its own copy of the old blind atom pattern
- Where: Adimensional.wl:64 `atom = x : (Derivative[__][_Symbol][___] | (h_Symbol …)[___])`. It is harmless today:
  `AdimAtom[..]` takes no argument and shapes are symbols. Step 0 fixes only Core.wl:85.
- Hint: reuse Core's rule, or note the restriction.
- Resolution: _(author)_

### E22 [nit] explicitVectors runtime on the existing n = 1 pipelines is not stated
- Evidence (P3): shallow-water DE4, 28 p² integrands: 4.9 s (d = 1) + 10.2 s (d = 2), against a 74 s build. That is +20% on
  every rebuild in the suite. The plan quotes only the KPZ p⁴ timings.
- Resolution: _(author)_

## Probes (inline)
Common preamble (PRE) for P1–P5: load the package, then patch the step-0 rule into the private context:
```
$HistoryLength=0; Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
FRGscalarMSR`$FRGExitOnFail=False; SetAttributes[quiet,HoldAll]; quiet[e_]:=Block[{Print},e];
testRoot="/home/greg/Physics/Sources/FRGscalarMSR/tests";
Begin["FRGscalarMSR`Private`"];
valueAt[e_, seed_Integer] := Module[{pt, memo, atomRule},
   pt = RandomRationalPoint[userSymbols[e], seed]; SeedRandom[seed + 7919];
   memo[x_] := memo[x] = RandomInteger[{2, 60}]/RandomInteger[{61, 97}];
   atomRule = x : (Derivative[__][_][___] | (_[___])[___] | (h_Symbol /; ! systemSymbolQ[h])[___]) :> memo[x];
   Quiet[Check[Together[(e /. pt) /. atomRule], $Failed], {Power::infy, Infinity::indet}]];
avg[k_,dim_,Q_]:=If[OddQ[k],0,Q^(k/2)(k-1)!!/Product[dim+2j,{j,0,k/2-1}]];
mut[k_,dim_,Q_]:=If[k==4,3Q^2/(7dim-6),avg[k,dim,Q]];
newProject[e0_,m_,n_,{nu_,p_},{w_,q_},dim_,deg_,eps_:False,tab_:avg]:=Module[{e=e0,cl},
 If[eps, e = dropWedge[e, m, n, {nu, p}, {w, q}]]; If[m>0,e=D[e,{nu,m}]];
 e=e/.{scal[p,p]->ss^2,scal[p,q]->ss cc}; e=(D[e,{ss,2n}]/.ss->0/.nu->0)/(2n)!; cl=CoefficientList[Expand[e],cc];
 e=Sum[cl[[k+1]] tab[k,dim,scal[q,q]],{k,0,Length[cl]-1}];
 e=e/.scal[q,q]->LoopMomentumSquared/.deltaCompletion[deg]; e/(I^m m!)];
lapn[e_,ps_,n_]:=Nest[Function[g,Sum[D[g,{x,2}],{x,ps}]],e,n]; Nn[n_,d_]:=4^n n! Pochhammer[d/2,n];
explicit[e0_,m_,n_,{nu_,p_},{w_,q_},d_,deg_]:=Module[{pc=Array[pcc,d],qc=Array[qcc,d],vec,e},
 vec={p->pc,q->qc};
 e=e0/.{scal[a_,b_]:>Expand[(a/.vec).(b/.vec)],wedge[a_,b_]:>((a/.vec)[[1]](b/.vec)[[2]]-(a/.vec)[[2]](b/.vec)[[1]])};
 If[m>0,e=D[e,{nu,m}]]; e=lapn[e,pc,n]/.pcc[_]->0/.nu->0; e=e/.deltaCompletion[deg]; {e/(I^m m!)/Nn[n,d], qc.qc}];
End[];
```
P1 (step 0). PRE; then `ExactlyZeroQ[Derivative[1][GammaAtom[hb,h,0]][QQ]]`, `…[Derivative[2][DeltaCoefficient[1]][QQ]]`,
`…[rr'[y]-rr[y]]` → False, False, False; `ExactlyZeroQ[D[GammaAtom[a,b,0][x^2],x]-2x Derivative[1][GammaAtom[a,b,0]][x^2]]` → True.
`Get[testRoot<>"/reference/KPZSetup.wl"]; RunOracle[kpzEquations]`, the tests/acceptance/Jacobian2D.wl pipeline up to
AssembleEquations, and `Get[testRoot<>"/reference/ShallowWaterSetup.wl"]` → all "ok", no abort (shallow water 74 s).

P2/P3 (invariant, explicitVectors). PRE; for every r in `ProjectOnMonomials[loops,sym]["Integrands"]`:
`Expand[newProject[LoopIntegrand[loops,r["Pair"],nu,p],r["Frequency"],r["Momentum"],{nu,p},LoopVariables[loops],dim,deg,eps]-r["Integrand"]]===0`.
For Momentum ≥ 1: `x=explicit[…,d,deg]; ExactlyZeroQ[x[[1]]-(r["Integrand"]/.LoopMomentumSquared->x[[2]]/.dim->d)]` →
KPZ 16/16, d = 1, 2 pass. jac 36/36, d = 2 pass (d = 1 is impossible with ε, as declared). Shallow water 56/56; d = 1 passes
(28 integrands, 4.9 s), d = 2 passes (10.2 s).

P4 (r1 E1/E4 re-attack). PRE + KPZSetup; `loop=LoopIntegrand[kpzLoops,{hbar,h},nu,p]`; n = 2, m = 0. Compare `explicit` at d =
1, 2, 3 with `newProject[…,Dim,deg,False,avg]` and with `…,mut]` →
`d=1: correct True, mutated False | d=2: True, True | d=3: True, True`.

P5 (E13). PRE; KPZ density as in KPZSetup, but with `DefineTheory["Scalars"->{h},"Dimension"->FixedDimension[2]]`,
`Neglected[_,…]`, orders {2, 1}; explicitVectors as in P3 at d = 1, 2 → d = 1 fails on all 8 p² integrands, d = 2 fails on none.

P6 (E14).
```
Get["…/FRGscalarMSR.wl"]; Begin["FRGscalarMSR`Private`"]; cnt=0; nw[]:=SpaceIndex[{"dummy",++cnt}];
Expand[toFactors[hb (grad[h].grad[h])^2, h|hb, nw]]   (* FieldFactor[h,0,{S[{dummy,1}]}]^4 hb *)
st=splitTerm[Expand[toFactors[hb lap[h] lap[h], h|hb, nw]]]; indexPairsQ[<|"Factors"->st[[2]]|>]   (* False *)
```
P7 (E15).
```
F[x_]:=1/(1+x);q={3/10,-1/2};d=2;Q=q.q;exact=F''[Q]/2+2Q F'''[Q]/d+2Q^2 F''''[Q]/(d(d+2));
st2[h_]:={{h,1/h^2},{0,-2/h^2},{-h,1/h^2}};u=IdentityMatrix[d];h=10^-5;
low[p_]:=Count[p,0]>=d-1 && Max[Abs[p]]<=h;
run[eL_,eH_]:=Module[{f},SeedRandom[3];f[p_]:=N[F[(p+q).(p+q)],80](1+If[low[p],eL,eH] RandomReal[{-1,1},WorkingPrecision->80]);
 N[Abs[(Total[Flatten[Table[a[[2]]b[[2]]f[a[[1]]u[[c1]]+b[[1]]u[[c2]]],{c1,d},{c2,d},{a,st2[h]},{b,st2[h]}]]]/(8d(d+2))-exact)/exact],3]];
{run[10^-34,10^-34], run[10^-22,10^-34]}   (* {1.4*10^-9, 1.7*10^-2} *)
```
P8 (E19). Nested stencil vs exact `lapn[F[(p+q).(p+q)],p,n]` at h = 10⁻⁵, k = 10/17, F ∈ {1/(1+(x/k²)²), e^{−x/k²},
2/(2+x/k²)}, n = 2, 3, d = 2, 3 → relative errors 6.0·10⁻¹⁰ … 3.8·10⁻⁹.

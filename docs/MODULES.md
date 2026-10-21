# Modules

Each module is a grey box. This file gives, for each one, its objective, its input and output stage, its mathematical
assumptions, the manipulations it performs, its checks, and its `Show…` commands.

## Behaviour common to every module

- **Stages.** Every stage is an immutable Association:
  `<|"Stage" -> name, "Pipeline" -> "DerivativeExpansion", "Fingerprint" -> hash, "Upstream" -> {upstream fingerprints}, data…,
  "Checks" -> {…}, "Warnings" -> {…}, "Timing" -> seconds|>`.
- **Input validation.** Every constructor checks the `"Stage"` key of its input and recomputes its fingerprint. A stage edited by
  hand, or one whose fields or couplings got a value after it was built, aborts with `Stage.fingerprint`.
- **Failed checks.** A failed check prints a red FAIL block (check name, module, short offending expression, hint) and then
  calls `Abort[]`. In a `math -script` run it calls `Exit[2]` instead, so headless runs exit non-zero. The name of the failed
  check is kept in `$FRGLastFailure`.
- **Summary cards.** Every constructor prints a card: one PASS line per check, warnings, LeafCount, timing and the first 8
  characters of the fingerprint.

---

## Core (`Kernel/Core.wl`)

**Objective.** Shared machinery: assertions, exact zero tests, stage bookkeeping and cards, scalar products.

**Contents.**
- `FRGAssert[name, cond, details, hint]` returns `name` when `cond` is True; otherwise it fails as described above. `details` and
  `hint` are evaluated only on failure.
- `ExactlyZeroQ[expr]` replaces `Simplify[…] == 0` (primer rule 7). It evaluates expr at 3 random rational points and requires
  exact 0 at each:
  - every non-System symbol (couplings, `q`, `k`, `Dim`, …) gets an independent random rational (E12: `Dim` is randomised too);
  - every remaining function atom, such as `r[y]`, `r'[y]`, `r''[y]`, `scal[a,b]`, gets its own independent random rational, so
    the identity is tested for arbitrary shape functions;
  - if the result is not an exact number, it aborts with `ExactlyZeroQ.numeric`.
- `RandomRationalPoint[syms, seed]`, `PrintCard[stage]`.
- `scal[a, b]` is the scalar product: Orderless, linear over sums and over explicit numbers only (E14). The rules come from
  `DE4plus_fast.wl` §1.
- `wedge[a, b]` = ε_ij a_i b_j in d = 2 (ε₁₂ = +1, design/LEVICIVITA.md). It is antisymmetric, `wedge[a, a]` = 0, and linear
  over sums and explicit numbers. It is kept in canonical order (`wedge[b, a]` → −`wedge[a, b]`), not Orderless.
- Internal `reduceWedge2D[e]` applies the d = 2 identity (F2)
  wedge[a,b] wedge[c,d] = scal[a,c] scal[b,d] − scal[a,d] scal[b,c] until each term has at most one wedge. With a symbolic d
  the p² projection of a wedge² term would be off by d − 1, which is why ε needs the declared d = 2.

**Assumption.** An identity that holds at 3 random rational points is treated as exact. A rational function of degree D in
its variables that vanishes at a random point drawn from M values does so by accident with probability ≤ D/M (the
Schwartz–Zippel bound), so three independent points make an accidental pass negligible.

---

## Theory (`Kernel/Theory.wl`) — `DefineTheory`

**Objective.** Declare the physical scalars, their response fields, the field order and the space dimension.

**Input / output.** Options `"Scalars" -> {φ1, …}` and `"Dimension" -> Dim` (a fresh symbol) or `"Dimension" -> FixedDimension[2]`
(needed by `jac`; the stage then stores `"Dimension" -> 2` and the card says "dimension 2 (declared)"). Output: the `"Theory"` stage, with keys
`Scalars`, `ResponseFields`, `Fields` (in field order), `ResponseOf`, `Dimension`.

**Assumptions.** N = 1 or 2. With 2 scalars, the propagator denominator Δ(ω) has degree 2 in ω, which is the largest degree
v1 integrates in closed form.

**Manipulations.** Creates the symbol `φbar`, in the same context as φ, for each scalar.

**Checks.**
- `Theory.N`: N ∈ {1, 2}.
- `Theory.scalars.*`, `Theory.responseFields.*`, `Theory.dimension.*`: the symbols are plain, not System symbols, not reserved
  (`k q p y t x ω ν`), do not collide with library names, have no definitions, and are not duplicated.
- `Theory.dimension`: a symbol or `FixedDimension[2]`. A literal number, `FixedDimension[3]`, or a `Dim` that got a value aborts.
- `Theory.dimensionDistinct`: `Dim` is not a field (symbolic d only).

**Show.** The card prints the field order.

---

## Ansatz (`Kernel/Ansatz.wl`) — `DefineAnsatz`

**Objective.** Turn the action density of Γ_k into a table of monomials, coefficient × Π(derivatives of single fields), that the
Vertices module Fourier-transforms (milestone 2).

**Input / output.** The Theory stage, plus `"Density" -> expr`, `"Couplings" -> {…}` and `"Pipeline" -> "DerivativeExpansion"`.
Output: the `"Ansatz"` stage. Its key `"Monomials"` is a list of
`<|"Coefficient", "Factors" -> {FieldFactor[φ, n_t, {indices}], …}, "Fields", "FieldCount", "ResponseCount", "TimeOrder",
"MomentumOrder"|>`. When a monomial carries ε (from `jac`), the stage also has `"Epsilon" -> True`; read it with
`Lookup[ansatz, "Epsilon", False]`.

**Assumptions.**
- Couplings depend on k only. They have no t or x dependence and never sit under an operator.
- The density is a polynomial in the fields, at most quartic.
- The momentum order of every monomial is at most 2. The p-Laplacian projection used downstream is exact only for p² (E2).
- MSR structure: every monomial contains a response field. There are no linear terms.

**Manipulations.**
1. Each vector sub-expression is given a rank, 0 for scalars and 1 for vectors. Ill-formed expressions abort, for example a
   product of two gradients without a dot, or `div` of a scalar.
2. The operators are expanded by the Leibniz rule until each acts on a single field. For example
   `div[H grad[Psi]] → grad[H]·grad[Psi] + H lap[Psi]`, and `lap[H Psi] → Psi lap[H] + 2 grad[H]·grad[Psi] + H lap[Psi]`.
   Gradients carry explicit spatial indices, and contractions are repeated indices.
3. The density is expanded into monomials. Since the momentum order is ≤ 2, there is at most one summed index, which is
   renamed `SpaceIndex[1]`; equal monomials then combine, and monomials that cancel are dropped (reported as INFO).
4. `jac[a, b]` puts `EpsIndex[1]` on the factors of a and `EpsIndex[2]` on those of b. ε is antisymmetric, so a monomial and
   its copy with the two slots swapped are the same monomial with opposite sign: the smaller of the two factor lists is kept,
   with sign −1 when the swapped one is chosen, and a list equal to its swap is zero (`jac[A, A]`, `jac[A,B] + jac[B,A]` cancel).
   The momentum order of a `jac` monomial is 2, so it has no other spatial index.

**Checks** (named as in the FAIL block).
- `Ansatz.noSystemDerivatives` (E7): no `Dt`, `D`, `Grad`, `Div`, `Laplacian` or `Derivative`.
  - Caveat: System operators applied to explicit coordinates (`Laplacian[H, {x, y}]`) evaluate to 0 before the library sees
    them. Only the package operators are safe.
- `Ansatz.knownSymbols`: every symbol is a field, a declared coupling or `Dim`. An explicit `t` or `x` gets its own hint.
- `Ansatz.couplingsOutsideOperators`: no `dt[nk H]` (every argument of `jac` is checked).
- `Ansatz.operatorArity` (`jac` takes 2 arguments, the others 1), `Ansatz.operatorOnField` (every argument contains a field).
- `Ansatz.epsilonDimension`: `jac` needs `"Dimension" -> FixedDimension[2]`.
- `Ansatz.singleEpsilon`: at most one `jac` per monomial, counted in slots (so `jac[…]^2` aborts). Two would be p⁴.
- `Ansatz.jacQuadratic`: a quadratic `jac` term is a total derivative, ∂_i(ε_ij A ∂_j B), whose vertex vanishes (F5): it aborts
  instead of leaving a coupling that silently does nothing.
- `Ansatz.indexPairs` (internal): every `SpaceIndex` appears twice; `EpsIndex[1]` and `EpsIndex[2]` once each or not at all.
- `Ansatz.vectorRank`, `Ansatz.vectorSum`, `Ansatz.vectorProduct`, `Ansatz.scalarDensity`: the vector algebra is well formed.
- `Ansatz.polynomial`, `Ansatz.unknownStructure`: the density is polynomial in the fields.
- `Ansatz.MSR`, `Ansatz.noLinear`, `Ansatz.quartic`, `Ansatz.momentumOrder` (E2).
- `Ansatz.pipeline`: only `"DerivativeExpansion"` exists in v1.
- Warning `Ansatz.unusedCouplings`: a declared coupling is absent from the density.

**Show.** `ShowMonomials[ansatz]` prints the table coefficient | fields | operators | order. It is printed automatically by
`DefineAnsatz`. `MonomialString[m]` writes a monomial back in the DSL (the ε pair as `jac[a, b]`).

---

## Vertices (`Kernel/Vertices.wl`) — `BuildVertices`

**Objective.** Fourier-transform the parsed density and build every vertex Γ⁽ⁿ⁾ (n = 2, 3, 4): the Γ2 matrix and the Γ3 (and Γ4)
functions in the leg convention of CONVENTIONS.md. This replaces the hand-written `Γ3[i]` and `AnsatzBMW` of `DE4plus_fast.wl`.

**Input / output.** The Ansatz stage. Output: the `"Vertices"` stage. Its key `"VertexTable"` maps each field multiset, sorted
in field order, to Γ⁽ⁿ⁾ as a function of all n leg momenta, before the δ is stripped. Accessors:
- `VertexValue[v, {a1,…,an}, {{ω1,q1},…,{ω_{n−1},q_{n−1}}}]`: leg j at (ωj, qj), last leg at −(sum of the others).
- `Gamma2Matrix[v, ω, q]`: entry (a, b) = Γ2_ab(ω,q), leg a at (ω,q), leg b at (−ω,−q). No regulator.
- `Gamma3Matrix[v, i, {ν,p}, {ω,q}]`: entry (a, b) = Γ⁽³⁾_{iab}, leg i at (ν,p), leg a at (ω,q), leg b at −(ν+ω, p+q). This is the
  script's `Γ3[i][ν,p,ω,q]`.

**Assumptions.** Those of the Ansatz. Couplings are constants, so each monomial is a product of independent plane-wave factors.

**Manipulations.**
1. A factor ∂t^{n_t} ∂_{i…} φ at leg l becomes (iω_l)^{n_t} Π(−i q_l,i). A repeated index between legs l and m becomes
   `scal[q_l, q_m]`, and within one factor `scal[q_l, q_l]`. The ε pair, `EpsIndex[1]` at leg l and `EpsIndex[2]` at leg m,
   becomes `wedge[q_l, q_m]`; with the (−i)² prefactor, `jac[A, B]` gives −`wedge[q_A, q_B]` (F3). Γ2 never carries a wedge,
   because `wedge[q, −q]` = 0.
2. The n-th functional derivative of c ∫dt dᵈx Π φ is the sum, over every assignment of the monomial's factors to legs that carry
   the same field, of c Π(symbols). This sum is the full symmetrisation: identical fields give the multiplicity factors (2 for
   φ̄², 3! for φ̄φ³).
3. The prefactors (2π)^{d+1} δ(Σω) δ(Σq) from ∫dt dᵈx and from the derivative normalisation cancel against those in the
   definition of Γ⁽ⁿ⁾. Nothing else remains: there is no 2π factor at vertex level.

**Checks.**
- `Vertices.planeWave` (independent). For every field multiset with n = 2..4, including those with no vertex, set
  φ_a(t,x) = Σ_l ε_l e^{i(ω_l t − q_l·x)}, with explicit 3-component vectors (2 for a declared d = 2, where `jac` becomes
  ∂₁a ∂₂b − ∂₂a ∂₁b and `wedge[a, b]` becomes a₁b₂ − a₂b₁), in the **raw density** using System `D`, `Grad`,
  `Div` and `Laplacian`. The coefficient of ε_1⋯ε_n must equal the derived vertex exactly. This tests the parser and the Fourier
  transform together, independently of their code.
- `Vertices.gamma2Symmetry`, `Vertices.causality` (Γ_φφ = 0), `Vertices.gamma2Polynomial`. These are implementation checks (E11):
  near-tautological by construction and labelled as such on the card.
- `Vertices.retardedPole[φ]` (E1) is the convention check.
  - It takes the free response G_{φ̄φ}(ω) = 1/Γ2_{φφ̄}(ω), with couplings at random positive values and q² = 1/2, and requires
    its pole in the upper half plane.
  - The time-mirrored convention or a damping term of the wrong sign aborts.
  - With no damping term (pole on the real axis) or no dt term, the check is inconclusive and reported as a WARN.
- WARN `Vertices.gamma4`: Γ⁽⁴⁾ ≠ 0. The tadpole term of the flow is hand-derived and untested. The bubble is tested against
  KPZ; see Loops.

**Show.** The vertex table, with the convention header, is printed by `BuildVertices`. The commands are `ShowVertexTable[v]`,
`ShowGamma2[v]` and `ShowVertex[v, {Hbar, H, Psi}]`.

**Acceptance (milestone 2).** In `tests/acceptance/ShallowWaterVertices.wl`, the derived Γ2 (regulator off) and Γ3[1..4] equal the
script's `Γ2 /. AnsatzBMW /. form*` and `Γ3[i] /. AnsatzBMW` exactly. The reference definitions are copied mechanically, with
renamed heads, into `tests/reference/DE4plusReference.wl`.

---

## Regulator (`Kernel/Regulator.wl`) — `DefineRegulator`

**Objective.** Declare the regulator R_k per field pair and derive ∂sR = k ∂k R. The shape functions stay abstract all the way
to the export.

**Input / output.** The Ansatz stage, plus rules `{a, b} -> R_ab` and `"Shapes" -> {r1, …}`. Output: the `"Regulator"` stage.
- `RegulatorMatrix[reg, s]` and `DsRegulatorMatrix[reg, s]` give the matrices at q² = s.
- `DsRAbstractMatrix[reg]` gives the matrix of abstract symbols `dsR[a, b]`.
- `DsRRules[reg, s]` gives the rules `dsR[a, b] -> derived entry`.

**Conventions (E10).**
- R_ab is the **Γ2-level entry**: what adds to Γ2_ab, including the factor 2 of a diagonal pair. For example,
  `{Hbar, Hbar} -> -2 Dhk q^2 armh[q^2/k^2]` gives (Γ2+R)_{H̄H̄} = −2Dhk q²(aDhk + armh).
- R is ω-independent and even in q, so R_ab = R_ba, and each unordered pair is given once.
- You write the entries with your own `q` (loop momentum) and `k` (RG scale); they are recognised by name. `k` remains the
  scale symbol in every later stage. Your `q` is stored as `"LoopMomentum"`, and the assembled equations are written in
  it.
- Map to `DE4plus_fast.wl`: dsNh = `dsR[H,Hbar]`, dsMh = `dsR[Hbar,Hbar]`, dsNphi = `dsR[Psi,Psibar]`,
  dsMphi = `dsR[Psibar,Psibar]`.

**Assumptions.**
- Couplings depend on k only.
- Shapes enter only as r(y), with y = q²/k².
- The overall amplitude is absorbed into r (Q5b): there is no separate amplitude symbol.

**Manipulations.** ∂sR = Σ_c (∂R/∂c) `dsCoupling[c]` + k ∂R/∂k at fixed q.
- `dsCoupling[c]` = k dc/dk. Downstream it becomes −η c for the couplings that carry anomalous dimensions.
- The chain rule gives k ∂k r(q²/k²) = −2 y r′(y).

**Checks.**
- `Regulator.pairFields` and `Regulator.responsePair`: only φ̄φ and φ̄φ̄ pairs, because a φφ entry would break causality.
- `Regulator.pairOnce` and `Regulator.evenInQ`.
- `Regulator.knownSymbols`: only couplings, q, k and shapes. In particular there is no frequency.
- `Regulator.onlyDeclaredShapes` and `Regulator.shapeArgument`: every shape appears only as r[q²/k²].
- `Regulator.dsRChainRule`: the same ∂sR is computed by a second route (couplings as functions c(k), total derivative) and the
  two results are compared exactly.

**Show.** `ShowRegulator[reg]` prints R_ab and ∂sR_ab. `DefineRegulator` calls it automatically.

---

## Propagator (`Kernel/Propagator.wl`) — `BuildPropagator`

**Objective.** Compute G = (Γ2 + R)⁻¹ symbolically and generically. This replaces the hand-written `G2Implicit`/`proprules`.

**Input / output.** The Vertices and Regulator stages, which must come from the same Ansatz. Output: the `"Propagator"` stage.
- `PropagatorMatrix[prop, ω, s]` gives G in atoms, with `Δ[±ω, s]`.
- `FullMatrix[prop, ω, s]` gives Γ2 + R in atoms.
- `AtomRules[prop]`, `DeltaRules[prop]` and `ValidityConditions[prop]` give the atom values, the Δ rules and the validity
  conditions.

**Abstract atoms (primer rule 1).**
- Each entry of Γ2 + R is a polynomial in iω. For a before b in field order, `GammaAtom[a, b, m][s]` is the coefficient of
  (iω)^m at q² = s. The (b, a) entry is the same polynomial evaluated at −ω.
- Numeric coefficients, such as the 1 in front of iω, stay numbers.
- The atoms get their explicit values (couplings, shapes) only at the very end.

**Mathematics.** In the order (φ1..φN | φ̄1..φ̄N), Γ2 + R = [[0, B(ω)], [B(−ω)ᵀ, C]], where:
- B is the (φ, φ̄) block, with rows φ at frequency ω;
- C is the (φ̄, φ̄) noise block;
- the (φ, φ) block is 0 by causality.

Then
- G_{φ̄φ} = B(ω)⁻¹ = adj B(ω)/Δ(ω). This is the retarded response, with poles at the roots of Δ(ω).
- G_{φφ̄} = (B(−ω)ᵀ)⁻¹, with denominator Δ(−ω).
- G_{φφ} = −(B(−ω)ᵀ)⁻¹ C B(ω)⁻¹, with denominator Δ(ω)Δ(−ω).
- G_{φ̄φ̄} = 0.

Here Δ(ω) = det B(ω) = c0 + c1 (iω) + c2 (iω)², with coefficients `DeltaCoefficient[m]`. For the shallow-water ansatz, c2 = 1,
c1 = 2A and c0 = CC of the script.

**Validity conditions.** All roots of Δ(ω) must lie in the upper half ω-plane (causal poles); every frequency integral assumes
this.
- Write u = iω. A root of Δ in the upper half ω-plane is a root u with negative real part.
- For degree 2, c2 u² + c1 u + c0 has both roots with Re u < 0 exactly when its three coefficients have the same sign: c1/c2 > 0
  and c0/c2 > 0. For the shallow-water ansatz this is A > 0 and CC > 0, as in the primer.
- For degree 1, the condition is c0/c1 > 0.
- The conditions are emitted as an ASSUMPTION, to be checked at run time at every q.

**Checks.**
- `Propagator.sameAnsatz`, `Propagator.polynomial`, `Propagator.atomsReal`, `Propagator.atomForm`, `Propagator.causalBlock`.
- `Propagator.deltaDegree`: 1 ≤ deg_ω Δ ≤ 2. Degree 3 or more is a future extension; degree 0 means the frequency integrals
  diverge.
- `Propagator.inverse` and `Propagator.inverseRight`: G·(Γ2+R) = (Γ2+R)·G = 1, exactly at random rationals for every atom, ω
  and q².
- `Propagator.responseBlock`.

**Show.** `ShowPropagator[prop]` prints the atoms, Δ, its coefficients and the validity conditions. `ShowPropagator[prop, {a, b}]`
prints G_ab.

**Acceptance (milestone 3).** In `tests/acceptance/ShallowWaterPropagator.wl`:
- the atoms (Γ_k + R) equal the script's `form*`, with the regulator included;
- Δ equals the script's Δ;
- {c0, c1, c2} = {CC, 2A, 1};
- G equals `G2Explicit /. AnsatzBMW` entry by entry, with Δ abstract on both sides;
- the dsR pattern equals `dsRmat` under the E10 map.

---

## Symmetries (`Kernel/Symmetries.wl`) — `DeclareSymmetries`

This module was added to the design's file list. `ProjectOnMonomials` needs the classification, and the classification needs
no flow, so it can fail early, before any loop is computed.

**Objective.** Classify every monomial (iω)^m (q²)^n, with m ≤ `"Frequency"` and n ≤ `"Momentum"`, of every Γ2 entry. The pair
is (a, b) with a before b in field order; the (b, a) entry is the same polynomial at −ω.

**Input / output.** The Ansatz stage, plus tags `Protected[pairPattern, "Frequency" -> m | {m..} | All, "Momentum" -> …,
"Reason" -> "…"]`, `Neglected[…]` with the same options, and `NonRenormalised[coupling, "Reason" -> "…"]`, plus
`"TruncationOrders" -> <|"Frequency" -> mf, "Momentum" -> mp|>`. Output: the `"Symmetries"` stage, whose `"Classification"`
lists one record per monomial: pair, m, n, class, Γ_k coefficient and reason.

**Classes**, in order of precedence:
- `Causality`: Γ_φφ between two physical fields. Its flow is computed and must vanish.
- `Structure`: odd m in a diagonal entry. Γ2_aa(ω) = Γ2_aa(−ω), so the monomial vanishes identically and is not computed (E5).
- `Protected`: the flow must vanish exactly, otherwise `CheckClosure` aborts with the reason.
  - A monomial that carries a coupling **and** is Protected is a contradiction and aborts (`Symmetries.contradiction`).
  - A numeric coefficient that is also Protected is consistent.
- `DeclaredNumeric`: a numeric coefficient in the ansatz, such as the 1 in front of iω. Its flow must vanish.
- `Declared`: a coefficient carrying couplings. It gives the equation ∂s(coefficient as written) = RHS (Q10).
- `Neglected`: the flow is computed, reported as a warning with its size, and dropped.
- Anything else aborts (`Symmetries.allTagged`).

**Pattern matching.** A pair pattern matches the unordered pair: `{Hbar, _}` matches every entry with an H̄ leg, including
`{Hbar, Hbar}`.

**Checks.**
- `Symmetries.momentumOrder` (mp ≤ 1, E2), `Symmetries.allTagged` and `Symmetries.contradiction`.
- `Symmetries.truncationCoversAnsatz`: no ansatz monomial lies beyond the truncation.
- `Symmetries.reason`: every tag has a reason.
- `Symmetries.nonRenormalised`.
- WARN for an unused tag.
- ASSUMPTION lines for `NonRenormalised` and for the Ward identities, which are declared, not derived.

**Show.** `ShowClassification[sym]` is printed automatically.

---

## Loops (`Kernel/Loops.wl`) — `BuildLoopIntegrands`

**Objective.** Build the one-loop integrand of ∂sΓ2_ab(ν,p) at a fixed loop point (ω, q).

**Formula.** With every ∂sR entry kept abstract (`dsR[c,d]`):

  ∂sΓ2_ab(ν,p) = ∫dω/(2π) ∫dᵈq/(2π)ᵈ { Tr[∂sR G(ω,q) Γ3_a(ν,p;ω,q) G(ω+ν,p+q) Γ3_b(−ν,−p;ω+ν,p+q) G(ω,q)]
                                         − ½ Tr[∂sR G(ω,q) Γ4_ab(ν,p;ω,q) G(ω,q)] }

- Γ3_a(ν,p;ω,q) is `Gamma3Matrix`: entry (c,d) = Γ⁽³⁾_{acd}, with leg a at (ν,p), leg c at (ω,q) and leg d at −(ν+ω,p+q).
- Γ4_ab(ν,p;ω,q) has entry (c,d) = Γ⁽⁴⁾_{abcd}, with legs at (ν,p), (−ν,−p), (ω,q) and (−ω,−q).
- The 1/(2π) of the ω-measure is included in the integrand. ∫dᵈq/(2π)ᵈ is left to the numerics.

**Hand derivation (E3/E4).** The bubble part is checked against KPZ (see below); the tadpole is not.
1. **Wetterich equation.** ∂sΓ = ½ Tr[∂sR (Γ⁽²⁾+R)⁻¹]. The fields are bosonic, so this is a plain trace. In Fourier space,
   Tr A = Σ_a ∫dω dᵈq/(2π)^{d+1} A_aa(ω,q).
2. **Inverting in Fourier space.** With the derivative normalisation of CONVENTIONS.md, the second derivative at a uniform
   background is (2π)^{d+1} δ(ω+ω′) δ(q+q′) Γ2_ab(ω,q). The operator inverse then has the same δ structure, with the matrix
   G(ω,q) = (Γ2(ω,q) + R(q))⁻¹.
3. **First field derivative.** δG/δφ_a = −G (δΓ⁽²⁾/δφ_a) G = −G Γ3_a G.
4. **Second field derivative.**
   δ²/δφ_a δφ_b of ½ Tr[∂sR G] = ½ Tr[∂sR (G Γ3_a G Γ3_b G + G Γ3_b G Γ3_a G − G Γ4_ab G)].
5. **The two cubic orderings are equal.** Shift the loop variables (ω,q) → −(ω+ν, q+p) and transpose inside the trace
   (Tr Mᵀ = Tr M), using:
   - ∂sR is symmetric;
   - G_cd(ω,q) = G_dc(−ω,−q), which follows from Γ2_cd(ω,q) = Γ2_dc(−ω,−q);
   - Γ3 is symmetric in its legs.

   The second ordering then becomes the first. Together, the two give the prefactor 1, and the tadpole keeps −½.
6. **Stripping the external δ.** Removing (2π)^{d+1} δ(ν+ν′) δ(p+p′) from both sides of the result of step 4 leaves the measure
   ∫dω dᵈq/(2π)^{d+1} = ∫dω/(2π) ∫dᵈq/(2π)ᵈ.

**What this is and is not tested against.**
- `DE4plus_fast.wl` and its oracle use exactly this formula (`d3Γ`, 1/(2π) included). Agreement with `flowsFAST` therefore
  checks the implementation, **not** the normalisation.
- **KPZ cross-check (`tests/acceptance/KPZ.wl`, milestone 6).** This test is independent of the formula. L. Canet,
  arXiv:2509.01472, Eq. (20), derives the flow of g = v_d k^{d−2} λ² D/ν³ for the same ansatz with the Litim cutoff
  r(y) = (1/y − 1) θ(1 − y). The same flow is used as the reference in `LPA_KPZ_thermal_forcing/lpa/kpz.jl`.
  - The test integrates our dimensionful ν and D equations over q with that cutoff, solves the linear system for η_ν and η_D,
    and compares d − 2 − η_D + 3η_ν with Canet's ratio. The agreement is exact, for symbolic d.
  - This fixes the bubble prefactor 1, its sign and the 2π. It also checks the residues, the frequency and momentum
    projections, and the assembly.
  - It does **not** test the Γ⁽⁴⁾ tadpole (−½), because KPZ has no Γ⁽⁴⁾. It tests only the combination
    d − 2 − η_D + 3η_ν, not each η separately; a compensating error in the ν and D equations is conceivable but would have
    to hold for every d and g.
- The card prints "tadpole … hand-derived, untested" whenever Γ⁽⁴⁾ ≠ 0.

**Checks.** `Loops.sameUpstream` and `Loops.dsRSymmetric`. The normalisation is printed on the card as an ASSUMPTION.

---

## Projection (`Kernel/Projection.wl`) — `ProjectOnMonomials`

**Objective.** For each classified monomial except `Structure`, turn the loop integrand into an ω-integrand at ν = p = 0. This
is pipeline 1: the external variables are set to zero. As a result the only ω-poles left are those of Δ(±ω).

**Operator.** The coefficient of (iν)^m (p²)^n is [(∂/∂(iν))^m / m!] applied to the integrand, followed, when n = 1, by
[∇²_p / (2d)], all evaluated at ν = 0, p = 0.
- The p-Laplacian uses the script's trick (§3): vectors are bare symbols, `D[D[f, p], p]` is taken, and the derivatives of
  `scal` are replaced by `Dimensionalfactors` (∂∂ scal[p,p] → d, and so on).
- A first derivative of scal[p,q] leaves a bare vector q. It is marked λ, so that the product of two of them is λ² = q·q, and
  the result is read off as (λ⁰ part) + ½ ∂²_λ (λ² part) × q². The λ¹ part must vanish; this is asserted.
- This is valid only for p² (E2).
- Derivatives of Δ(ω, s) are completed with Δ = Σ_m c_m(s) (iω)^m, so they come out in `DeltaCoefficient[m]` atoms and their
  s-derivatives.
- **With ε** (`jac`, declared d = 2 only), before anything else:
  - `reduceWedge2D` (F2) leaves at most one wedge per term, so every term is wedge[p,q]^e f(p², p·q, q², ν) with e ∈ {0, 1}
    (wedge[q, q] = wedge[p, p] = 0 and linearity reduce any other wedge to this one);
  - the e = 1 terms are dropped (F4): they vanish at p = 0, and so does their p-Laplacian, 2 ε_ij q_j q_i ∂_{p·q} f = 0. This
    is exact for every projection with m, n ≤ 1 and any m.

**Checks.**
- `Projection.LHS`: the same operator, applied to Γ_k (regulator off), reproduces every classified Γ_k coefficient. This checks
  the normalisation of both projections, including 2d.
- `Projection.evenInQ` and `Projection.complete`: no external variable, vector, wedge, or derivative of `scal` or Δ is left.
- With ε: `Projection.epsilonDimension` (d = 2), and `Projection.wedgeOdd`. After F2 every remaining wedge must be
  ±wedge[p, q] to a power ≤ 1, and the projection of the dropped part, evaluated with explicit 2-vectors at a random rational
  point, must be exactly 0. The card counts the integrands with a nonzero dropped part.

**Show.** `ShowIntegrand[projected, nk*Ak]` or `ShowIntegrand[projected, {a, b}, {m, n}]`.

---

## FrequencyIntegration (`Kernel/FrequencyIntegration.wl`) — `IntegrateFrequencies`

**Objective.** Integrate every projected integrand over the real ω axis, exactly. The result contains no ω, no Δ and no square
roots. This replaces the notebook's lookup table and generalises the script's `Iint`/`ωIntegrate`.

**Mathematics.**
- Δ(ω) = c0 + c1 (iω) + c2 (iω)². Its roots z1, z2 lie in the upper half plane (the propagator's validity conditions). Those of
  Δ(−ω) are −z1, −z2, in the lower half plane.
- `FrequencyIntegral[n, a, b, deg]` = ∫ ω^n Δ(ω)^{−a} Δ(−ω)^{−b} dω = 2πi × (sum of the residues at z1, z2), obtained by closing
  the contour in the upper half plane.
- Write z1,2 = iα ± β, with α = c1/(2c2) and β² = (4 c0 c2 − c1²)/(4 c2²). The residue sum is symmetric under z1 ↔ z2, so it is
  even in β, and replacing β² by its value makes it rational in c0, c1, c2.
- This covers overdamped motion (β² < 0) and underdamped motion (β² > 0) alike.
- For degree 1, the single root is z = i c0/c1.
- For a = 0 there is no pole in the upper half plane, so the integral is 0.
- Convergence requires n − deg·(a+b) ≤ −2; otherwise it aborts.
- Dispatch is by structure: Δ(ω) → 1/X and Δ(−ω) → 1/Y, then `CoefficientRules` over {ω, X, Y}. Every monomial ω^n X^a Y^b is
  integrated by the memoised `FrequencyIntegral`. Nothing is skipped silently, and nothing depends on rule order.

**Modes.**
- `Declared` monomials are integrated symbolically, in `GammaAtom`, `DeltaCoefficient` and `dsR` atoms.
- `Protected`, `DeclaredNumeric` and `Causality` monomials are integrated **exactly at 3 random rational points**. The Δ
  coefficients are computed from the atoms at each point, so the relation Δ = det B is kept; ω and Δ stay symbolic.
- `Neglected` monomials are integrated at one point that is **valid**: q² > 0, and the propagator's validity conditions hold
  (every root of Δ(ω) in the upper half plane). Points are redrawn until one is valid, at most 200 times, after which the
  module aborts with `FrequencyIntegration.validPoint`. A value reported as a size must come from a physical point. Zero
  tests need no such point, because an exact zero is a zero everywhere.
- `CheckClosure` asserts the zeros and reports the sizes.

**Checks.** `FrequencyIntegration.convergence`, `.evenInBeta`, `.rational`, `.ready`, `.complete` and `.validPoint`. The unit tests compare
the results with the script's I(0,1,1) = π/(2AC) and with `NIntegrate`, in the overdamped, underdamped and degree-1 cases.

**Acceptance (milestone 4).** In `tests/acceptance/ShallowWaterFlows.wl`:
- All 9 shallow-water equations equal `flowsFAST/*.wdx` exactly, at random rationals for every symbol and every shape atom.
  This holds both with the abstract dsR mapped to dsNh, … (E10) and with the derived ∂sR substituted on both sides.
- HEQ ≡ 0.
- Every Protected, DeclaredNumeric and Causality flow samples to exactly 0.
- The three frequency terms of script §6b are generated.
- The full suite runs in about 26 s.

---

## Closure (`Kernel/Closure.wl`) — `CheckClosure`

**Objective.** Decide, monomial by monomial, whether the integrated flows close on the ansatz. Then turn each declared
monomial into a flow equation.

**Input / output.** Input: the Integrated stage, plus the Symmetries stage it was projected with. Using any other Symmetries
stage aborts with `Closure.sameSymmetries`. Output: the "Closure" stage, with these keys:
- `Equations`: the declared records, each with `Normalisation`, `Flowing` and `RHS`;
- `Report`: every classified monomial, including those that vanish by structure.

**What happens to each class.**
- **Protected**, **DeclaredNumeric** and **Causality**: all 3 exact samples must be 0. If one is not, the module aborts
  with `Closure.protected`, `Closure.declaredNumeric` or `Closure.causality`; the hint lists the offending monomials and the
  reasons of their tags. A nonzero Protected flow means one of two things: the declared symmetry is not a symmetry of the
  ansatz with this regulator, or the ansatz breaks it.
- **Neglected**: one WARN line per monomial, with its value at a valid random point (see FrequencyIntegration). This value
  is one sample, not a bound: it says only that the term is generated.
- **Structure**: an INFO count. These monomials are not computed.
- **Declared**: the Γ_k coefficient is split into its pure number (sign included) and the rest, the *flowing combination*
  (Q10, agreed with the user). For example, −2 Dhk aDhk gives normalisation −2 and combination Dhk aDhk. The equation is
  ∂s(combination) = Value/normalisation.
  - The (2d)ⁿ of the p² projection is already in Value. The normalisations therefore equal the LHS factors in the script's
    `jobs` table divided by (2 Dim)ⁿ: −4Dim → −2, −2Dim → −1.

**Assumptions printed.** Every `NonRenormalised` tag with its reason (not checkable in v1, because that needs ∂sΓ⁽³⁾), and
the banner "Ward identities are assumed (the Protected tags), not derived".

**Show.** `ShowClosureReport[closure]` prints a table pair | monomial | class | outcome, preceded by a plain-words
explanation of each class.

**Acceptance.** `tests/acceptance/ShallowWaterClosure.wl`:
- the 9 combinations and normalisations of the script;
- script §6b: iω and ω² in Γ_{HΨ̄} and ω² in Γ_{Ψ̄Ψ̄} are reported as generated.

`tests/unit/ClosureTest.wl` (KPZ) covers the split and the negative test of design §10: a protected term with nonzero flow
(ω² in Γ_{h̄h̄} tagged Protected) aborts. The untagged-monomial negative test is in `SymmetriesTest.wl`, because the
monomial is caught there, before any flow is computed.

---

## Assemble (`Kernel/Assemble.wl`) — `AssembleEquations`

**Objective.** Write every equation in two forms:
- **compact**: scale-stripped atoms and abstract ∂sR entries, as the script's `compact.wdx`;
- **explicit**: the q-integrand in couplings, `dsCoupling[c]`, shapes and `Dim`, ready for FlowRG (design §6 item 10). For a
  declared d = 2 the integrand has the number 2 in place of `Dim`, and the card says "d = 2 (declared)".

**Scale-stripped atoms.** These are generated automatically, replacing the script's hand-written `toAtoms`/`atomVals`.
- Each propagator atom v(s) = GammaAtom[a,b,m](s), with s = q², is written v = c·sⁿ·X. Here c is the number and coupling
  factor common to every term of v, n is the smallest power of s among those terms, and X is what is left.
- When X is a sum, it becomes the atom `ScaleAtom[a,b,m,0]`; when it is a single term, it stays explicit.
- The j-th s-derivative is written c·X_j with X_j = dʲ(sⁿX)/dsʲ. It becomes `ScaleAtom[a,b,m,j]` when it is a sum.
- In shallow water this reproduces the script's atoms exactly: xE = Ak + r, xE1 = Ak + r + y r′, xE2 = (2r′ + y r″)/k², and
  so on, with xL = gk + g2k q² and fl′ = g2k left explicit.
- Check `Assemble.atomSync`: atom form = explicit value, exactly, for every atom and every derivative used.

**Manipulations.**
1. The Δ coefficients and their s-derivatives are replaced by their expressions in GammaAtom atoms.
2. GammaAtom atoms are replaced by their atom forms, with q² = the user's `q`^2 (the regulator's loop-momentum symbol).
3. The result is collected over the dsR entries and brought to canonical form with `Together`. It is checked to be linear
   and homogeneous in dsR (`Assemble.linearInDsR`): the loop is proportional to ∂sR.
4. `compactify`, ported from script §5. The candidates are per-dsR or whole-fraction forms, with the numerator content
   pulled out and nested `Collect` over the atoms: derivative atoms first, then the others in two orders. Small results also
   get a time-limited `Simplify` per dsR coefficient. The smallest explicit LeafCount wins, and it is checked exactly against
   the canonical form (`Assemble.compactExact`).
5. Explicit form: atoms → values, then the derived ∂sR (`DsRRules`) is substituted only here. `Assemble.dsRLinear` checks
   that ∂sR is linear in the `dsCoupling` symbols. `Assemble.noIntermediates` checks that no GammaAtom, ScaleAtom, Δ,
   DeltaCoefficient, dsR, `LoopMomentumSquared` or private symbol remains.

**Output contract per equation.**
- `Flowing`, `Normalisation`, `Pair`, `Frequency`, `Momentum`;
- `Compact`, `CompactForm` (which candidate won), `DsR`;
- `Integrand`: the ∫dᵈq/(2π)ᵈ integrand of ∂s(Flowing);
- `DsCouplings`, `IdenticallyZero`, `LeafCount`. The compact LeafCount counts each ScaleAtom and dsR as one leaf, as the
  script's xE and dsNh do.

The stage also carries:
- `Atoms`: ScaleAtom → explicit value;
- `AtomForms`;
- `ValidityConditions`: `Atoms` and `Explicit`; they must hold at every q, which is checked at run time downstream;
- `DsRRules`, `LoopMomentum`, `Scale`, `Shapes`, `Couplings`, `Dimension`.

Shapes appear as r[q²/k²], r′ and r″. FlowRG supports only r′ today (docs/FUTURE.md).

**Show.** `ShowEquation[equations, comb]` prints the compact form, the atoms it uses, the ∂sR entries and the validity
conditions. `EquationOf[equations, comb]` returns the equation record.

**Acceptance.** `tests/acceptance/ShallowWaterEquations.wl`:
- all 9 `Integrand`s equal `flowsFAST/*.wdx` exactly, with the derived ∂sR on both sides;
- compact LeafCounts are no larger than 1.5× the script's `compact.wdx`. In practice they are close or smaller: g2 28 306 vs
  30 130, D2phi 111 465 vs 115 792, Dphi·aDphi 2 141 vs 1 572;
- HEQ ≡ 0;
- assembly takes about 35 s, most of it (20 s) on D2phi.

**LPA acceptance** (`tests/acceptance/ShallowWaterLPA.wl`). The LPA ansatz, DE4plus without the Kk, g2k and D2phik terms,
is derived from scratch, with those three flows reported as Neglected. Its equations equal the old notebooks' exactly:
- `flowsDE4plus_agreed` /. {K, D2ϕ, g2 → 0}, for g, Dϕ, N, Dh and A. This is the input of
  `LPA_KPZ_thermal_forcing/AdimensionalFlowEQGenericReg.nb`.
- `flowsLPAold` with λ = 1, A = 0 and no H-noise regulator, for g, Dϕ, N and Dh, with the LHS factors of
  `shallow_water_DE4plus.nb` §14.

The old-notation map and the provenance are in `tests/reference/legacy/README.md`.

**DE4 against the old notebooks** (`tests/acceptance/ShallowWaterDE4Legacy.wl`).
- Exact agreement with `flowsDE4plus_agreed` for g, Dϕ, N, Dh, A and K.
- Exact agreement with the notebook's `flowsNEW` for g, Dϕ, N, Dh, A, K, H and g2.
- Not compared:
  - the agreed HEQ, which holds only a placeholder symbol;
  - the flowsNEW D2ϕ file, which is corrupt;
  - the agreed g2 and D2ϕ. These contain unsubstituted f″ and, even when completed, disagree with this library, which
    agrees with flowsNEW, flowsFAST and the oracle. See `tests/reference/legacy/README.md`.

---

## Oracle (`Kernel/Oracle.wl`) — `RunOracle`

**Objective.** A numerical cross-check of the assembled equations, by a route independent of Propagator, Projection,
FrequencyIntegration and Assemble. It generalises script §6c.

**Method.** At one random point:
- couplings and `dsCoupling[c]` are random rationals and k is random;
- each shape is a concrete function (1/(1+y²), e^{−y}, 2/(2+y), …);
- d = 2 and d = 3 (`"Dimensions" -> Automatic`); for a declared d = 2, only d = 2. An explicit list must be a subset of
  {2, 3} for a symbolic d, and exactly {2} for a declared d = 2; anything else aborts with `Oracle.declaredDimension`;
- the loop momentum is an explicit d-vector; `wedge[a, b]` becomes a₁b₂ − a₂b₁ once p is explicit.

The point is redrawn until the validity conditions hold, because the symbolic side assumes causal poles. Then, for each
equation:
1. Γ2 + R is built numerically as a matrix polynomial in ω and inverted as a matrix. Γ3, Γ4 and ∂sR are evaluated
   numerically.
2. The loop trace (the formula of Loops) is a rational function of ω with numeric coefficients. It is integrated by
   `NIntegrate` over the real axis, with 40 digits. An absolute accuracy goal stops the search on integrands that vanish.
3. The projections are taken by central finite differences with h = 10⁻⁵:
   - the p-Laplacian over the d unit vectors, divided by 2d;
   - d^m/dν^m, divided by iᵐ m!;
   - the normalisation is applied.
4. The result is compared with `Integrand` at the same point: relative error < 10⁻⁸, or absolute for zero flows.
   `"ZeroFlows" -> True` also checks every Protected, DeclaredNumeric and Causality flow. This exercises the ν
   finite differences and gives an independent numerical check of the Ward identities.

**What it does not test.** It shares the vertices and the loop formula (prefactor, −½ tadpole, 1/(2π)) with the symbolic
chain, so it does not test the normalisation (E3/E4). The card says so. The bubble normalisation is tested by the KPZ comparison
with Canet's flow instead (see Loops).

**Acceptance.**
- Shallow water, 9 equations in d = 2 and 3: relative errors ≤ 2·10⁻⁹. This is the finite-difference error, about h²; the
  script used a tolerance of 10⁻⁶. HEQ is numerically 5·10⁻¹³. Runs in 11 s.
- KPZ with zero flows (`AssembleTest.wl`): 20 comparisons in 2 s.
- `jac` in d = 2 (`tests/acceptance/Jacobian2D.wl`): 5 equations, relative errors ≤ 10⁻⁹.

---

## Stages (`Kernel/Stages.wl`) — `SaveStage`, `LoadStage`

**Objective.** Optional persistence of any stage (Q13).
- `SaveStage[stage, file]` writes the stage, with every stage it embeds, as WDX.
- `LoadStage[file, "Expect" -> {stages…}]` reads it back. It recomputes the fingerprint (a modified file aborts with
  `Stage.fingerprint`), then checks that the fingerprint of every expected stage appears among the saved stage's
  fingerprints and upstream fingerprints (`StageFingerprints`). Pass, for example, the current theory, ansatz, regulator and
  symmetries; a mismatch aborts with `Stages.stale`.

**Tests.** `tests/unit/StagesTest.wl`:
- a round trip gives the identical stage, which is usable downstream;
- an ansatz with an extra coupling makes the saved stage stale;
- a modified file and a missing file abort.

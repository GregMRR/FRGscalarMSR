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

**Assumption.** An identity that holds at 3 random rational points is treated as exact. A rational function of degree D in
its variables that vanishes at a random point drawn from M values does so by accident with probability ≤ D/M (the
Schwartz–Zippel bound), so three independent points make an accidental pass negligible.

---

## Theory (`Kernel/Theory.wl`) — `DefineTheory`

**Objective.** Declare the physical scalars, their response fields, the field order and the space dimension.

**Input / output.** Options `"Scalars" -> {φ1, …}` and `"Dimension" -> Dim`. Output: the `"Theory"` stage, with keys
`Scalars`, `ResponseFields`, `Fields` (in field order), `ResponseOf`, `Dimension`.

**Assumptions.** N = 1 or 2. With 2 scalars, the propagator denominator Δ(ω) has degree 2 in ω, which is the largest degree
v1 integrates in closed form.

**Manipulations.** Creates the symbol `φbar`, in the same context as φ, for each scalar.

**Checks.**
- `Theory.N`: N ∈ {1, 2}.
- `Theory.scalars.*`, `Theory.responseFields.*`, `Theory.dimension.*`: the symbols are plain, not System symbols, not reserved
  (`k q p y t x ω ν`), do not collide with library names, have no definitions, and are not duplicated.
- `Theory.dimensionDistinct`: `Dim` is not a field.

**Show.** The card prints the field order.

---

## Ansatz (`Kernel/Ansatz.wl`) — `DefineAnsatz`

**Objective.** Turn the action density of Γ_k into a table of monomials, coefficient × Π(derivatives of single fields), that the
Vertices module Fourier-transforms (milestone 2).

**Input / output.** The Theory stage, plus `"Density" -> expr`, `"Couplings" -> {…}` and `"Pipeline" -> "DerivativeExpansion"`.
Output: the `"Ansatz"` stage. Its key `"Monomials"` is a list of
`<|"Coefficient", "Factors" -> {FieldFactor[φ, n_t, {indices}], …}, "Fields", "FieldCount", "ResponseCount", "TimeOrder",
"MomentumOrder"|>`.

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

**Checks** (named as in the FAIL block).
- `Ansatz.noSystemDerivatives` (E7): no `Dt`, `D`, `Grad`, `Div`, `Laplacian` or `Derivative`.
  - Caveat: System operators applied to explicit coordinates (`Laplacian[H, {x, y}]`) evaluate to 0 before the library sees
    them. Only the package operators are safe.
- `Ansatz.knownSymbols`: every symbol is a field, a declared coupling or `Dim`. An explicit `t` or `x` gets its own hint.
- `Ansatz.couplingsOutsideOperators`: no `dt[nk H]`.
- `Ansatz.operatorArity`, `Ansatz.operatorOnField`.
- `Ansatz.vectorRank`, `Ansatz.vectorSum`, `Ansatz.vectorProduct`, `Ansatz.scalarDensity`: the vector algebra is well formed.
- `Ansatz.polynomial`, `Ansatz.unknownStructure`: the density is polynomial in the fields.
- `Ansatz.MSR`, `Ansatz.noLinear`, `Ansatz.quartic`, `Ansatz.momentumOrder` (E2).
- `Ansatz.pipeline`: only `"DerivativeExpansion"` exists in v1.
- Warning `Ansatz.unusedCouplings`: a declared coupling is absent from the density.

**Show.** `ShowMonomials[ansatz]` prints the table coefficient | fields | operators | order. It is printed automatically by
`DefineAnsatz`. `MonomialString[m]` writes a monomial back in the DSL.

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
   `scal[q_l, q_m]`, and within one factor `scal[q_l, q_l]`.
2. The n-th functional derivative of c ∫dt dᵈx Π φ is the sum, over every assignment of the monomial's factors to legs that carry
   the same field, of c Π(symbols). This sum is the full symmetrisation: identical fields give the multiplicity factors (2 for
   φ̄², 3! for φ̄φ³).
3. The prefactors (2π)^{d+1} δ(Σω) δ(Σq) from ∫dt dᵈx and from the derivative normalisation cancel against those in the
   definition of Γ⁽ⁿ⁾. Nothing else remains: there is no 2π factor at vertex level.

**Checks.**
- `Vertices.planeWave` (independent). For every field multiset with n = 2..4, including those with no vertex, set
  φ_a(t,x) = Σ_l ε_l e^{i(ω_l t − q_l·x)}, with explicit 3-component vectors, in the **raw density** using System `D`, `Grad`,
  `Div` and `Laplacian`. The coefficient of ε_1⋯ε_n must equal the derived vertex exactly. This tests the parser and the Fourier
  transform together, independently of their code.
- `Vertices.gamma2Symmetry`, `Vertices.causality` (Γ_φφ = 0), `Vertices.gamma2Polynomial`. These are implementation checks (E11):
  near-tautological by construction and labelled as such on the card.
- `Vertices.retardedPole[φ]` (E1) is the convention check.
  - It takes the free response G_{φ̄φ}(ω) = 1/Γ2_{φφ̄}(ω), with couplings at random positive values and q² = 1/2, and requires
    its pole in the upper half plane.
  - The time-mirrored convention or a damping term of the wrong sign aborts.
  - With no damping term (pole on the real axis) or no dt term, the check is inconclusive and reported as a WARN.
- WARN `Vertices.gamma4`: Γ⁽⁴⁾ ≠ 0. The tadpole term of the flow is hand-derived and untested.

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
  scale symbol in every later stage.
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

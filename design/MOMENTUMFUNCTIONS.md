# Momentum-dependent couplings (functional Γ2 slots) — design

**Status: implemented** (except the KPZ-LO literature benchmark, see "Implementation notes"). Code: every module of the cascade;
tests `tests/unit/MomentumFunctionsTest.wl`, setup `tests/reference/KPZLOSetup.wl`. Decisions settled with the user in a grilling session (2026-10-23) and four
adversarial rounds by the-enemy: `design/reviews/MomentumFunctions-r1.md` (2 blocker, 5 major, 6 minor, 1 nit), `-r2.md`
(1 blocker, 2 major, 6 minor, 1 nit), `-r3.md` (1 blocker, 3 major, 4 minor), `-r4.md` (0 blocker, 2 major; minors not sought).
Probes: `design/reviews/MomentumFunctions-r2-probe-*.wl`, `-r3-probe-*.wl`. Labels [E#] refer to those findings; each is
resolved in the item that carries it. Do not reopen the decisions.

## Context
Today every coupling depends on k only and flows are projected at ν = p = 0 (derivative expansion). We add couplings that are
functions of p² = |p|² in the quadratic (Γ2) part of the ansatz, as in the LO/BMW-type KPZ scheme (f_ν(p), f_D(p)). Their flow is
computed at generic external p: the loop carries propagators on two momentum lines, so the frequency integral has two pole
families and the angular integral stays explicit. Output: dimensionless equations; there is no FlowRG export (another library
integrates them), so `ExportJulia` refuses functional couplings.

## Reasons for the main choices
- Γ2-only functions (A1): a function in a local cubic term has no unique momentum argument; Γ2 legs carry ±p.
- Polynomial frequency (A1): keeps the residue integration exact (FrequencyIntegration needs rational ω-integrands).
- Function × reference (A3): the dimensionless function is then F̂ = F, and η stays a number attached to a constant reference.
- Lower orders as ordinary monomials (A4) [E2]: a generated mass m² would otherwise hide in F as m²/p²; today's tags keep meaning.
- Midpoint routing + r → −r symmetrisation (C1) [E15, E25]: the only routing where the subtracted integrand is O(p^{2n₀})
  pointwise for every n₀; ∂sR sits on one line only, so symmetrisation is still needed.
- ŷ = 0 objects by the p → 0 operator (C4) [E27]: same integral as the midpoint row at ap → 0 but a different integrand, one fewer
  shape derivative (Litim δ terms survive), and no 0/0 (E33).
- `Related` on functions, η balance derived (D7) [E16, E28]: the user states the relation once; the library derives what it
  implies and imposes it, instead of asking the user to retype the regulator.

## Implementation notes (deviations from the text below, and why)
- **Loop variable.** The midpoint loop momentum is the regulator's `q` (not a new `r`): q → k aq in Adimensional is unchanged. Below,
  read r, ar as q, aq. The cosine is `c`, the dimensionless external momentum `ap` (both reserved names now).
- **Record keys.** A functional record keeps `"Momentum" -> n₀` (an integer: the dimension checks use it) and carries `"Functional" -> F`,
  `"Argument" -> p^2`, `"Variables" -> {q, p, c}`, `"AtZero"` (the exact p → 0 row). `EquationOf[eqs, Fnu]` and `ShowBeta[adim, Fnu]` work.
- **Compact form.** Generic-p rows keep the atoms applied to their momentum (`ScaleAtom[a, b, m, j][x]`, `AdimAtom[a, b, m, j][x]`,
  `dsR[a, b][x]`) and are collected over the ∂sR entries only (no `compactify` candidates).
- **Vertices.planeWave** replaces each function by a constant on both sides (the argument F(q²) is a definition), instead of a
  polynomial test operator: with independent leg momenta the operand choice would matter off-shell.
- **d = 1 measure.** `FunctionalMeasure[1, …]` is not defined (0·∞); the stage stores `FunctionalMeasure -> 1/(2 Pi)` and
  `FunctionalCosines -> {-1, 1}` when d = 1 (e.g. after `Related[…, "AtDimension" -> 1]`).
- **All-Fixed η system.** When no reference is Running (the default with functions), the η system has no linear part; Adimensional now
  pads it with a zero matrix (it aborted with `Adimensional.etaLinear` before).
- **Benchmark (D15/E3) not done.** The general-d LO equations of Canet, Chaté, Delamotte, Wschebor, PRL 104, 150601 (2010)
  (arXiv:0905.1025, Eqs. 14–17; the 2011 PRE defers to them) as printed do not reproduce the p → 0 limit: their η_D integrand at p = 0
  differs from the derivative-expansion one (which the KPZ acceptance test fixes against Canet 2025, Eq. 20) by its structure (powers of
  γ̃11 + r and of q̃²). Eq. 15 also lacks the 1/p̃² of the γ11 = Γ11/p² projection. A pointwise or integrated comparison therefore fails for
  reasons in the reference, not in the derivation; the library side is checked by the oracle at generic p and the p → 0 regression.

## Notation
- Dimensionful p (external), r (loop, midpoint routing: lines at r − p/2 and r + p/2), s± = |r ± p/2|².
- Dimensionless ap = p/k, ar = r/k, ŷ = ap², c = cosθ = p·r/(p r). F̂(ŷ) = dimensionless function.
- Slot = (Γ2 entry (a,b), frequency order m). n₀ = power of p² multiplying F in its slot.
- "Generic-p rows": functional slots and `"AllMomenta"`-tagged slots. "p → 0 rows": everything else, and every ŷ = 0 object.

## A. Writing the ansatz
- A1 Only couplings of quadratic monomials may be functions of p² (Γ3/Γ4 constant). Frequency stays polynomial (iν)^m, projected
  at ν = 0. Mixed ansätze allowed; no functions ⇒ pipeline identical to today.
- A2 `DefineAnsatz[th, "Density" -> hbar (dt[h] - nk Fnu lap[h] - lamk/2 grad[h].grad[h]) - Dk FD hbar^2,
  "Couplings" -> {nk, lamk, Dk}, "MomentumFunctions" -> {Fnu, FD}]`. Functions are declared only in `"MomentumFunctions"` [E6]
  (fresh-symbol / not-a-field checks as couplings); `"Pipeline"` set automatically, a contradicting user value aborts with hint
  "drop "Pipeline"" [E11]. Functions appear bare (`Ansatz.functionBare`) [E8], only in quadratic terms
  (`Ansatz.functionQuadraticOnly`), outside operators (`Ansatz.functionNotInOperator`).
- A3 Functional term shape [E4, E21]: checked on the Γ2 slot coefficient after Vertices (`Vertices.functionalTerm`, slot quoted):
  number × (at most one constant coupling) × F¹ × (p²)^{n₀}. The constant coupling must be a reference
  (`Dimensions.functionPrefactor`); otherwise F is dimensionful, or declared in `"Dimensionless"` (accepted when consistent with
  the derived dimension) [E17]. `Dimensions.freedoms` hint names MomentumFunctions.
- A4 Slot structure [E2]: Σ_{n<n₀} c_n p^{2n} + p^{2n₀}·[ref]·F(p²). Lower orders are ordinary monomials, classified at p → 0 as
  today. A constant at power ≥ n₀, or a second function, aborts (`Symmetries.functionalSlot`). Each function sits in exactly one
  slot with nonzero Γ2 (`Symmetries.functionSlot`: two slots → hint `Related`; zero slots, e.g. `Fnu div[hbar grad[h]]`) [E5].
- A5 Every option that takes couplings rejects functions with a named check: `Regulator.noFunctions`,
  `Symmetries.nonRenormalisedFunction`, `Dimensions.functionReference`, `Adimensional.functionEquation` (except the `"Equation"`
  of A7) [E6, E22].

## B. Symmetries and tags
- B1 Orders ≥ n₀ of a functional slot: one row of class `Functional` on the card; TruncationOrders "Momentum" ignores it (card line)
  and applies to lower orders only (`truncationCoversAnsatz` adapted) [E20, E11].
- B2 Tags at ≥ n₀ [E29]: only an explicit integer `"Momentum" -> n₀` matches order n₀ (`All` and patterns never do); meaning: the
  loop part of the slot's flow at p = 0 vanishes (∂s(ref·F(0)) has no loop contribution), checked with the p → 0 operator. Any
  other tag at ≥ n₀ aborts `Symmetries.functionalTag`. The Galilean tag keeps today's spelling when λ → f^λ.
- B3 `Protected[pair, "Frequency" -> m, "Momentum" -> "AllMomenta"]` / `Neglected[…, "AllMomenta"]` on function-free slots act at
  generic p (exact zero / WARN with size). On a functional slot → `Symmetries.contradiction` (hint updated). [E18]
- B4 Sibling slot (same pair, other m) of a functional slot without an `"AllMomenta"` tag: ASSUMPTION
  `Symmetries.siblingAtZeroMomentum` ("computed at p → 0 only") — makes LO's f^λ ≡ 1 visible [E20].
- B5 Constant numeric slots without tags: exact zero at p → 0 only (today); no generic-p report [E12].

## C. Computing generic-p rows
- C1 Midpoint routing r (lines r ∓ p/2) [E15]; ∂sR carries its line's momentum before anything else, `dsR[c,d][s]` [E26];
  the integrand is symmetrised under r → −r (≡ p → −p at fixed r) [E25], making it even in p at fixed r. Constant slots keep
  today's routing.
- C2 Lower orders subtracted pointwise with their un-averaged Taylor coefficients in p at fixed r (they integrate to flow_n);
  divide by p^{2n₀}. `Projection.functionalOrder`: the result is O(1) pointwise as p → 0 (quotes the first offending order) [E1].
- C3 Frequency integration with two pole families Δ(ω, s−), Δ(ω, s+) (X±, Y±); residues of both families combined and cancelled
  before storing [E3]; `.evenInBeta`/`.rational` per family; `FrequencyIntegration.coincidentFinite` (finite at s+ → s−, which
  happens inside the domain at c = 0). Unit tests vs NIntegrate incl. s+ = s−, overdamped, underdamped, degree 1.
- C4 Every ŷ = 0 object (η rows of A7, B2 checks, the p → 0 regression) uses today's p → 0 operator, not a midpoint row at ap → 0
  (same integral, different integrand and shape-derivative order) [E27].

## D. Equations and dimensionless form
- D1 Closure: functional equations ∂s([ref]·F(p²)) = RHS(p, r, c); `Closure.protectedFunctional` (B3), B2 checks.
- D2 Assemble: `ScaleAtom[a,b,m,j][s]` at s±; integrand in (r, p, c); `noIntermediates` extended. Records: `"Momentum" ->
  "Functional"`, `"Argument"`; `EquationOf[eqs, Fnu]` [E23].
- D3 Adimensional: p → k ap, r → k ar, F(p²) → [F] F̂(ŷ); BetaL gains 2ŷ ∂ŷ F̂ (`linearNoMomentum` → free of ar; ap allowed in
  functional rows) [E11]; `AdimAtom[a,b,m,j][arg]`; F̂ named `aFnu` (override `"Names"`); functions in `"State"` by default
  (declaration order); `ShowBeta[adim, Fnu]`; homogeneity/unit-invariance extended to function arguments; causal poles with F̂.
- D4 Measure stored per row [E10]: generic-p rows K_d ar^{d−1}(1−c²)^{(d−3)/2} dc, K_d = S_{d−1}/(2π)^d; d = 1 stored as the
  separate case Σ_{c=±1}; p → 0 rows keep today's LoopMeasure.
- D5 Normalisation: F̂(0) free by default. References are constant; a reference appearing only times functions is Fixed (η at
  run time, `MustFix`). `nk -> {"nu", "Equation" -> nk Fnu}` imposes F̂(0) = 1, η from the p → 0 operator (C4) [E1, E19].
- D6 `"Equation"` always means Running, for `nk ank` and `nk Fnu` alike; `"Mode" -> "Fixed"` with `"Equation"` aborts
  (`Dimensions.equationModeFixed`); `"Mode" -> "Running"` accepted, redundant [E30]. (Fixes today's silent ignore.)
- D7 `Related` (an `AdimensionalEquations` option) [E16, E28]:
  `"Related" -> {Related[FD -> Fnu, "AtDimension" -> 1, "Reason" -> "time reversal (d = 1)"]}`.
  - FD eliminated (Rule semantics); optional number α̂ (`FD -> 2 Fnu`), dimensionless (`Adimensional.relatedDimension`).
  - Derived, imposed, printed as ASSUMPTION: η balance of the two slots' references (η_D = η_ν; η_D = 0 if one side has none);
    the shape relation from the Regulator stage (imposed; abort `Adimensional.relatedRegulator` only on a prefactor mismatch);
    Dim → AtDimension in every row when given (abort if Theory fixed another d; see D8); initial condition FD̂ = α̂ Fnû for the downstream integrator.
  - Check `Adimensional.related`: β_FD̂ = α̂ β_Fnû pointwise (c sampling per D8/E34).
  - Named aborts: function related to a constant coupling, both sides in one slot, chains/cycles, eliminated function in "State".

- D8 Round-4 resolutions (`design/reviews/MomentumFunctions-r4.md`):
  - [E33] The functional row (I − lower orders)/ap^{2n₀} is exact for ŷ > 0 but 0/0 at ŷ = 0 and loses float precision near it
    (≈ 10⁻¹⁶/ap^{2n₀}; matters for n₀ ≥ 2). The library never evaluates it there (C4), so equations and checks are unaffected:
    INFO line `Adimensional.functionalSmallMomentum` on the card for n₀ ≥ 1, and an entry in `docs/FUTURE.md` (possible fix:
    store an "AtZero" row + next-order series for downstream use). No other change.
  - [E34] `"AtDimension"` semantics: = 1 → check at c = ±1; = n ≠ 1 → Dim → n, random rational c ∈ (−1, 1) plus ±1; omitted →
    every d (symbolic Dim, random c). A relation holding only after the angular integral aborts with a hint saying so. The
    ∂sR line attribution (E26) is guarded by the Oracle only (stated in the design record).

## E. Export, Oracle, benchmark
- E1 `ExportJulia` aborts `JuliaExport.functionalCouplings` on any momentum/frequency-dependent coupling.
- E2 Oracle at generic (r, p, c): numeric G on both lines, ∂sR on its own line, r → −r symmetrisation, NIntegrate over ω, explicit
  smooth F̂ (e.g. 1/(1+y)); includes a point with s+ = s−.
- E3 `benchmarks/KPZ_LO/` vs Canet–Chaté–Delamotte–Wschebor LO equations (f^λ ≡ 1), transcribed from arXiv (PRE 84, 061128 /
  PRL 104, 150601) citing equation numbers and errata [E13]; reference shifted to midpoint routing and r → −r symmetrised before
  pointwise comparison. p → 0 regression: functional KPZ's ŷ = 0 rows (C4) equal today's KPZ equations exactly.

## Check → hint table (goes into the design record) [E22, E32]
| check | trigger | hint |
|---|---|---|
| Ansatz.functionQuadraticOnly | `Fnu hbar h^2` | functions multiply quadratic terms only (Γ3/Γ4 couplings stay constant) |
| Ansatz.functionNotInOperator | `lap[Fnu h]` | write the function outside the operators: `Fnu lap[h]` |
| Ansatz.functionBare | `Fnu[0]` | write `Fnu` bare: it already means Fnu(p²) |
| Ansatz.pipeline (adapted) | contradicting `"Pipeline"` | drop "Pipeline": it is set from "MomentumFunctions" |
| Vertices.functionalTerm | `nk Ak Fnu`, `Fnu^2`, `1/Fnu`, F at two powers | one term ref × Fnu × (p²)^n₀: absorb the extra factor into the function |
| Symmetries.functionalSlot | constant at ≥ n₀, two functions | `-nk hbar lap[h]` is F(0)'s job: drop it / keep one function |
| Symmetries.functionSlot | function in 0 or 2 slots | term vanishes in Γ2: drop it / one function per slot, use `Related` |
| Symmetries.functionalTag | tag at > n₀, Neglected at n₀ | only `"Momentum" -> n₀` (Protected) is allowed on a functional slot |
| Symmetries.contradiction (adapted) | AllMomenta tag on a functional slot | use `"Momentum" -> n₀` for F(0), or `Related` |
| Symmetries.nonRenormalisedFunction | `NonRenormalised[Fnu]` | NonRenormalised takes constant couplings |
| Regulator.noFunctions | `nk Fnu q^2 rn` | regulators take constant couplings: use the reference `nk` |
| Projection.functionalOrder | lower orders not cancelled | quotes slot and first non-vanishing order p^{2j}, j < n₀: tag or declare it |
| FrequencyIntegration.coincidentFinite | residues not cancelled at s+ = s− | internal: report |
| Dimensions.freedoms (adapted) | undetermined function | …or declare the MomentumFunction in "Dimensionless" |
| Dimensions.functionPrefactor | function times a non-reference | make that coupling a reference or absorb it into the function |
| Dimensions.functionReference | function in "References" | references are constant couplings |
| Dimensions.equationModeFixed | Mode Fixed + Equation | "Equation" makes the reference Running: drop one |
| Adimensional.functionEquation | function as a bare "Equation" | write `ref Fnu` |
| Adimensional.relatedDimension | non-dimensionless α̂, or Theory d ≠ AtDimension | … |
| Adimensional.relatedRegulator | regulator prefactors not related | write R with the same prefactor structure in both slots |
| Adimensional.related | β identity fails at c = ±1 | the relation is not a symmetry of this ansatz/regulator |
| JuliaExport.functionalCouplings | any function | functional equations are not exportable to FlowRG |

## Cascade per module
Theory (reserve ap, ar, c) · Ansatz (A2) · Vertices (F → F[scal[q,q]]; `gamma2Polynomial` → polynomial in iω with function-valued
coefficients; planeWave with independent legs, F as polynomial test operator on either operand; A3) [E11] · Regulator (A5) ·
Propagator (atom extraction with function-valued coefficients; validity conditions with F, ASSUMPTION) · Symmetries (A4, B) ·
Loops (∂sR line argument, C1; midpoint routing for generic-p rows) · Projection (C1, C2, C4; `LHS`, `complete` adapted) ·
FrequencyIntegration (C3) · Closure (D1) · Assemble (D2) · Dimensions (A3, A5, D6) · Adimensional (D3–D7) · JuliaExport (E1) ·
Oracle (E2) · Stages (nothing).

## Docs, tests
- `design/MOMENTUMFUNCTIONS.md` (this content + reasons + hint table) and `design/reviews/MOMENTUMFUNCTIONS/`.
- `docs/MODULES.md`, `docs/CONVENTIONS.md` (ap, ar, c, ŷ, F̂, routing, measure), `docs/FUTURE.md`, CLAUDE.md map.
- Lessons for Ansatz, Symmetries, Loops, Projection, FrequencyIntegration, Adimensional, JuliaExport; rebuild.
- Unit tests per module, `expectFail` for every new check; `tests/acceptance/KPZ_LO.wl`; p → 0 regression; D6 behaviour change test.

## Order of work
1. Design record + reviews. 2. Ansatz/Vertices/Symmetries/Regulator checks. 3. FrequencyIntegration two families. 4. Loops line
argument + midpoint, Projection. 5. Closure. 6. Assemble, Oracle. 7. Dimensions/Adimensional (incl. D6, D7). 8. JuliaExport.
9. KPZ_LO benchmark + regression. 10. Docs, lessons.

## Verification
- `math -script tests/runTests.wl` exit 0.
- `math -script benchmarks/KPZ_LO/run.wl`, `benchmarks/KuramotoSivashinsky/run.wl` exit 0.
- `RunOracle` at generic p and at s+ = s− on KPZ_LO.
- `math -script lessons/makeLessons.wl`; touched lessons headless.

# FRGscalarMSR

A Mathematica library that takes an ansatz for Γ_k of N scalar fields in the MSR (response-field) formalism and produces a
closed, checked set of Wetterich flow equations, using the derivative-expansion pipeline.
- The design is in `design/FRGDESIGN.md` (decision records for later additions: `design/LEVICIVITA.md`, `design/ADIMENSIONAL.md`).
- Each module is described in `docs/MODULES.md`: what it does, its assumptions and its checks.
- The conventions are in `docs/CONVENTIONS.md`: Fourier transform, vertex legs, the density language.
- Planned extensions are in `docs/FUTURE.md`.

**Status: v1 complete (milestones 1–6).**
- Shallow water: the 9 DE4plus equations reproduce `DE4plus_fast.wl`'s `flowsFAST` exactly, and a numerical oracle confirms them
  in d = 2 and 3.
- LPA: the LPA ansatz (no K, g2, D2phi) reproduces the old notebooks' LPA equations exactly (`flowsDE4plus_agreed` reduced to
  LPA, and the original `flowsLPAold`).
- DE4 vs the old notebooks: exact agreement with `flowsDE4plus_agreed` (g, Dϕ, N, Dh, A, K) and `flowsNEW` (g, Dϕ, N, Dh, A, K, H, g2).
  The old agreed g2 and D2ϕ disagree and are considered wrong (`tests/reference/legacy/README.md`).
- KPZ: the one-loop flow reproduces L. Canet's closed-form flow (arXiv:2509.01472, Eq. 20) exactly, for any d.
- 2D Levi-Civita term `jac` (d = 2 declared): the equations agree with the numerical oracle in d = 2 (`tests/acceptance/Jacobian2D.wl`).
- Dimensionless flows and FlowRG export (`Dimensions`, `Adimensional`, `JuliaExport`): dimensions derived from the ansatz, KPZ equal
  to Canet's Eq. 20 after adimensionalisation, generated files evaluated in Julia and loaded by FlowRG.

## Install

- Mathematica (Wolfram Language 14 or later; the tests use the `math` command-line kernel, Wolfram 15).
- Nothing to build: load the library from a notebook or a script with `Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"]`.
- Optional, for running the exported flows: Julia and the FlowRG package (`~/Physics/Sources/FlowRG.jl`).

## Usage

The library has two parts:
1. **Flow equations**: from an ansatz for Γ_k to the dimensionful equations ∂s(coupling) = ∫dᵈq/(2π)ᵈ (integrand).
2. **Dimensionless flows and export**: dimensions derived from the ansatz, dimensionless β functions and anomalous dimensions,
   and a FlowRG equations file.

Every step returns an immutable *stage* (an Association) and prints a summary card: one PASS line per check, plus warnings and
assumptions. A failed check prints a red FAIL block with a hint and stops the evaluation; headless runs exit with a non-zero code.
Full expressions are printed only on request, with the `Show…` commands.

### Quickstart: KPZ from the ansatz to the FlowRG file

The whole chain for the KPZ equation ∂t h = ν ∇²h + λ/2 (∇h)² + noise (variance 2D), about 5 s:

```mathematica
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"]

(* part 1: flow equations *)
theory     = DefineTheory["Scalars" -> {h}, "Dimension" -> Dim];
ansatz     = DefineAnsatz[theory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2,
               "Couplings" -> {nuk, lamk, Dk}];
vertices   = BuildVertices[ansatz];
regulator  = DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
prop       = BuildPropagator[vertices, regulator];
symmetries = DeclareSymmetries[ansatz,
               Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
               Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
               Neglected[_, "Reason" -> "beyond the ansatz"],
               NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
               "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>];
loops      = BuildLoopIntegrands[vertices, prop, regulator];
projected  = ProjectOnMonomials[loops, symmetries];
integrated = IntegrateFrequencies[projected];
flows      = CheckClosure[integrated, symmetries];
equations  = AssembleEquations[flows];
ShowEquation[equations, nuk]

(* part 2: dimensionless flows and the FlowRG file *)
dims = DeriveDimensions[ansatz, "References" -> {nuk -> "nu", Dk -> "D"}];
adim = AdimensionalEquations[equations, dims];
ShowBeta[adim, lamk]
ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{$TemporaryDirectory, "kpz", "Equations2D.jl"}]]
```

and in Julia:

```julia
using FlowRG
m = load_model("/tmp/kpz/Equations2D.jl"; regulator = Litim)
etas(m, [0.7]), beta(m, [0.7])          # (nu = 9.5e-5, D = 0.0197), [-0.00679]
```

The two-field shallow-water case is in `examples/ShallowWaterDE4plus.nb`.

### Part 1, step by step: the flow equations

1. **Fields.** `DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]` declares N = 1 or 2 physical scalars and creates the
   response field `<name>bar` of each (`Hbar`, `Psibar`). `Dim` is the space dimension, a symbol kept symbolic to the end (or
   `FixedDimension[2]`, see below).
2. **Ansatz.** `DefineAnsatz[theory, "Density" -> …, "Couplings" -> {…}]`: Γ_k = ∫dt dᵈx (density). Write the density with the
   library's own operators `dt[…]`, `grad[…]`, `div[…]`, `lap[…]`, `jac[A, B]` and `grad[a].grad[b]` for a scalar product (not the
   System `Dt`, `Grad`: `Dt[nk H]` would differentiate the coupling). Couplings depend on k only and stand outside the operators
   (`nk dt[H]`, never `dt[nk H]`). Every monomial needs a response field; at most quartic in the fields, at most p² in momentum.
   The parsed monomials are printed: check them.
3. **Vertices.** `BuildVertices[ansatz]` Fourier-transforms every monomial (∂t → +iω, ∇ → −iq, all momenta incoming) and prints
   every vertex under that convention, so you can check it by hand (`docs/CONVENTIONS.md`).
4. **Regulator.** `DefineRegulator[ansatz, {a, b} -> R_ab, …, "Shapes" -> {r1, …}]`. Each entry is the amount **added to Γ2_ab**,
   including the factor 2 of a diagonal pair: `Dk hbar^2` in the density gives Γ2_{hbar hbar} = −2 Dk, so its regulator is
   `-2 Dk rd[q^2/k^2]`. Shapes are abstract functions of q²/k²; they stay abstract up to the Julia file, where the regulator is
   chosen at run time. Only φ̄φ and φ̄φ̄ pairs.
5. **Propagator.** `BuildPropagator[vertices, regulator]` inverts Γ2 + R exactly and prints the conditions for causal poles.
6. **Symmetries.** `DeclareSymmetries[ansatz, tags…, "TruncationOrders" -> <|"Frequency" -> mf, "Momentum" -> mp|>]` gives every
   monomial (iω)^m (p²)^n of every Γ2 entry, up to the truncation orders, one class:
   - declared in the ansatz: it gets a flow equation;
   - `Protected[{a, b}, "Frequency" -> m, "Momentum" -> n, "Reason" -> "…"]`: a symmetry forbids it; its flow is checked to be
     exactly 0;
   - `Neglected[…, "Reason" -> "…"]`: generated but dropped; its size is reported;
   - zero by structure (e.g. odd ω powers of a diagonal entry): detected, no tag needed.
   `{a, b}` may be a pattern (`{Hbar, _}`, `_`); "Frequency"/"Momentum" take a number, a list or `All` (default). Precedence:
   Protected > declared > Neglected. An untagged monomial stops the evaluation. `NonRenormalised[c, "Reason" -> "…"]` declares a
   vertex coupling that does not flow (not checkable here: it would need ∂sΓ⁽³⁾).
7. **Loops to equations.** `BuildLoopIntegrands` → `ProjectOnMonomials` (external frequency and momentum → 0) →
   `IntegrateFrequencies` (exact residues) → `CheckClosure` (protected flows vanish, dropped terms reported) →
   `AssembleEquations`.
8. **Inspect and check.** `ShowEquation[equations, comb]`, `EquationOf[equations, comb]`, `ShowMonomials[ansatz]`,
   `ShowVertex[vertices, {Hbar, H, Psi}]`, `ShowGamma2[vertices]`, `ShowRegulator[regulator]`, `ShowPropagator[prop, {H, Psibar}]`,
   `ShowClassification[symmetries]`, `ShowIntegrand[projected, nk Ak]`, `ShowClosureReport[flows]`. `RunOracle[equations]`
   recomputes every equation numerically (matrix inverse, numerical frequency integral, finite differences) in d = 2 and 3.
9. **Save.** `SaveStage[equations, file]`, `LoadStage[file, "Expect" -> {theory, ansatz, regulator, symmetries}]` (the
   fingerprints must match, else it aborts).

The output equations are ∂s(coefficient as written in the ansatz) = ∫dᵈq/(2π)ᵈ (integrand), e.g. ∂s(nk Ak). In each integrand:
- `dsCoupling[c]` = k dc/dk enters linearly;
- the regulator shapes r(q²/k²) and their derivatives stay abstract;
- `Dim` stays symbolic, except with a declared d = 2 (next section).

### 2D Levi-Civita term: `jac` and `FixedDimension[2]`

ε_ij ∂_iA ∂_jB (ε₁₂ = +1) is written `jac[A, B]`, with A and B scalar expressions of the fields. ε exists only in d = 2, so the
dimension must be declared; it is never fixed silently (a literal `2`, or a `Dim` that has the value 2, aborts):

```mathematica
theory = DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> FixedDimension[2]]      (* card: "dimension 2 (declared)" *)
ansatz = DefineAnsatz[theory, "Density" ->
    Hbar (dt[H] + lk jac[Psi, H] + kap grad[Psi].grad[H] - nk lap[H]) - dk Hbar^2 + Psibar (dt[Psi] - mk lap[Psi]) - ek Psibar^2,
  "Couplings" -> {lk, kap, nk, dk, mk, ek}]
```

- `jac[A, B]` = −`jac[B, A]`; terms that cancel this way (e.g. `jac[A, A]`) are dropped and reported.
- At most one `jac` per monomial (two would be order p⁴), and no quadratic `jac` term (it is a total derivative): both abort.
- In the vertices, `jac[A, B]` with A at leg l and B at leg m becomes −`wedge[q_l, q_m]`, where `wedge[a, b]` = ε_ij a_i b_j
  (printed `wedge(a,b)`). Γ2 never carries a wedge.
- The rest of the pipeline is unchanged. The equations carry the number 2 in place of `Dim`, and `RunOracle` defaults to
  `"Dimensions" -> Automatic`: {2, 3} for a symbolic d, {2} for a declared d = 2. An explicit list must be a subset of {2, 3}
  (symbolic d) or exactly {2} (declared d = 2), e.g. `RunOracle[equations, "Dimensions" -> {3}]` for a symbolic d only.
- Details: `docs/CONVENTIONS.md`, `docs/MODULES.md` (Ansatz, Projection, Oracle) and the design note `design/LEVICIVITA.md`.

### Part 2, step by step: dimensionless flows and the FlowRG file

1. **Dimensions.** `DeriveDimensions[ansatz, "References" -> {…}, "Dimensionless" -> {…}]` derives the dimension of every
   coupling, field and of time from the ansatz (Γ_k is dimensionless, [x] = 1/k), written as k^x ref₁^y₁ … in the basis of k and
   **N+1 reference couplings** that you choose (N for the fields, one for time; for shallow water `nk`, `Dhk`, `Dphik`). Each
   reference names its anomalous dimension: `nk -> "nu"` gives η_ν (symbol `etanu`), with ∂s nk = −η_ν nk.
   - Declare in `"Dimensionless"` the amplitudes that share a monomial with another coupling (`Ak` in `nk Ak`, `aDphik` in
     `Dphik aDphik`): one monomial cannot split the dimension between them. If you forget, the error names them.
   - Checked: exactly N+1 free units, independent references, every monomial consistent, and the action invariant under the
     derived rescaling.
   - `ShowDimensions[dims]` prints coupling | monomials | dimension | dimensionless name | linear part of its flow, then the field
     dimensions, [ω] and z (e.g. [ω] = k² nk, z = 2 − η_ν). `CheckDimensions[dims, <|gk -> k Dphik^(1/2) Dhk^(-1/2) nk, …|>]`
     compares your hand values exactly.
2. **Dimensionless flows.** `AdimensionalEquations[equations, dims, "State" -> {gk, lk, …}]`: the state couplings (dimensionful
   names, in the order you want in Julia) become `agk`, `alk`, … (declared dimensionless couplings keep their name), with
   ∂s ĉ = BetaL + ∫daq LoopMeasure BetaNL: BetaL is the canonical part ĉ(−x + Σ y_j η_j) plus local terms, **not integrated**;
   BetaNL is the loop integrand, with LoopMeasure = S_d aq^{d−1}/(2π)^d and q = k aq. The η's solve a linear system with the same
   split. Every β is checked to be free of k and of the references.
   - A reference that appears **alone** in an equation (`nk` in `- nk lap[Psi]`) is *running*: its η comes from that equation.
   - A reference that appears **only times an amplitude** (`Dphik aDphik`) is *fixed*: you fix its η at run time (e.g. a forcing
     whose energy injection must stay finite); the equation becomes the flow of the amplitude, and the file lists the η in
     `MUST_FIX`. `"References" -> {Dphik -> {"phi", "Mode" -> "Running"}}` instead takes η from the amplitude equation with the
     amplitude counted as 1 (the convention of the old hand-built exports; add `"Equation" -> nk ank` when several qualify).
   - `NonRenormalised` couplings keep their canonical part only.
   - `ShowBeta[adim, gk]` and `ShowEtaSystem[adim]` print them.
3. **Export.** `ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{dir, "eqs", "Equations2D.jl"}]]` writes one file per
   space dimension, following FlowRG's `docs/equations-interface.md` (COUPLINGS, ETAS, SHAPES, `betal`, `betanl`, `etaMl`, …). Load
   it with `load_model(path; regulator = Litim, fixed_etas = (phi = 0.0,))`, giving every `MUST_FIX` eta.

From a stage saved by part 1 (e.g. in another notebook):

```mathematica
equations = LoadStage[FileNameJoin[{Dir, "swLPAplusequations.wdx"}]];
dims = DeriveDimensions[AnsatzOf[equations], "References" -> {nk -> "nu", Dhk -> "h", Dphik -> "phi"}, "Dimensionless" -> {Ak, aDphik}];
adim = AdimensionalEquations[equations, dims, "State" -> {gk, lk, Hk, Ak, mk, Dmk, aDphik}];
ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{Dir, "eqs", "Equations2D.jl"}]]
```

Limits of part 2:
- Regulator prefactors may contain reference couplings only (`nk q^2 arnh`, not `nk (1 + Ak) q^2 arnh`): ∂sR would otherwise
  contain a β itself. Lifting this is planned (`docs/FUTURE.md`).
- Flows of momentum-derivative couplings (DE4plus g2, D2φ) need r″, written `D2armh(aq^2)`; FlowRG does not bind it yet
  (`docs/handoff/FlowRG-MUST_FIX.md`, which also specifies the `MUST_FIX` guard).
- Large files (DE4plus, ~1 MB) take minutes for Julia to compile on first use.

Checks: KPZ reproduces Canet's Eq. (20) exactly after adimensionalisation (`tests/acceptance/KPZAdimensional.wl`); DE4plus agrees
with the old hand-built `Equations2D.jl` up to its two known errors (`tests/acceptance/ShallowWaterDE4Adimensional.wl`). Design and
decisions: `design/ADIMENSIONAL.md`.

## Examples

`examples/ShallowWaterDE4plus.nb` and `examples/KPZ.nb` are generated by `examples/makeNotebooks.wl`. Edit the cells there and
regenerate. To check that a notebook runs without a front end:

```
math -script examples/runNotebook.wl examples/KPZ.nb                  # ~10 s
math -script examples/runNotebook.wl examples/ShallowWaterDE4plus.nb  # ~2.5 min
```

## Tests

```
math -script tests/runTests.wl                          # unit + acceptance, ~6 min; exit code 0 means every test passed
math -script tests/runTests.wl --only DimensionsTest    # only the test files whose name contains DimensionsTest
math -script tests/runTests.wl --verbose                # also print the cards and FAIL blocks
```

Tests that need `julia` or FlowRG (`~/Physics/Sources/FlowRG.jl`, or `FLOWRG_PATH`) print SKIP when those are missing;
`--require-julia` turns a SKIP into a failure.

## Known limitations

- **Ward identities** are assumed (declared by the user as `Protected`), not derived. `NonRenormalised` vertex couplings are
  not checkable in v1, because that needs ∂sΓ⁽³⁾.
- **Loop normalisation (E3/E4).** The trace formula is hand-derived in `docs/MODULES.md` (Loops).
  - Its bubble part (prefactor, sign, 2π) is checked independently against Canet's KPZ flow (`tests/acceptance/KPZ.wl`).
  - The −½ Γ⁽⁴⁾ tadpole is **untested**, and the cards say so whenever Γ⁽⁴⁾ ≠ 0.
  - Agreement with `flowsFAST` and with the oracle does not test the normalisation, because they share the formula.
- **Truncation.** Momentum projections only up to p² (E2). Δ (the determinant whose roots are the propagator poles) must be of
  degree ≤ 2 in ω, so at most 2 scalars. Both limits abort clearly; the extensions are in `docs/FUTURE.md`.
- **Validity conditions.** The causal-pole conditions must hold at every loop momentum; they are emitted with the equations
  and are to be checked at run time.

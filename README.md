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
5. **Propagator.** `BuildPropagator[vertices, regulator]` inverts Γ2 + R exactly and prints the conditions for causal poles
   (exported later, dimensionless, as `causal_poles`).
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
   - **Causal poles.** The propagator's validity conditions (all roots of Δ(ω) in the upper half plane) are rewritten in
     minimal dimensionless form, in `adim["CausalPoles"]`: denominators cleared, positive factors (numbers, powers of k, of the
     references, of aq) stripped, assuming the references are > 0. Shallow water gives `Ak + ank + arnh(aq²) + arnphi(aq²) > 0`
     and `agk aHk alk + aq² (ag2k aHk alk + (Ak + arnh)(ank + arnphi)) > 0`; KPZ gives `1 + rn(aq²) > 0` (shape positivity is
     not used: the shapes are chosen in Julia). Each must hold at every aq. `ShowCausalPoles[adim]` prints them
     (design/CAUSALITY.md).
3. **Export.** `ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{dir, "eqs", "Equations2D.jl"}]]` writes one file per
   space dimension, following FlowRG's `docs/equations-interface.md` (COUPLINGS, ETAS, SHAPES, `betal`, `betanl`, `etaMl`, …). Load
   it with `load_model(path; regulator = Litim, fixed_etas = (phi = 0.0,))`, giving every `MUST_FIX` eta.
   - The file also defines `causal_poles(couplings…, aq)`, the causal-pole conditions: each entry must be > 0 at every aq (only
     the sign means something). FlowRG does not call it: evaluate it yourself, e.g. `M.causal_poles(u..., aq)` on a grid of aq,
     with the shapes bound as in FlowRG. Omitted when no condition is left.

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

## Function reference

Exported functions per module, in load order. `docs/MODULES.md` gives each module's assumptions, manipulations and checks. Vectors
in `VertexValue`/`LoopIntegrand` are symbols, scalar products are `scal[a, b]`. `s` stands for q², `w` for the frequency.

**Theory**

| Call | Arguments | Returns / prints |
|---|---|---|
| `DefineTheory[…]` | `"Scalars" -> {φ1, …}`, `"Dimension" -> Dim` or `FixedDimension[2]` | Theory stage; field order (φ1, φ1bar, φ2, …) |
| `ResponseField[theory, φ]` | theory stage, a scalar | its response field `φbar` |

**Ansatz** (operators for the density: `dt`, `grad`, `div`, `lap`, `jac`)

| Call | Arguments | Returns / prints |
|---|---|---|
| `DefineAnsatz[theory, …]` | `"Density" -> density`, `"Couplings" -> {c1, …}` | Ansatz stage; prints the parsed monomials |
| `ShowMonomials[ansatz]` | ansatz stage | table coefficient / fields / operators |
| `MonomialString[monomial]` | one parsed monomial | its field part written in the DSL |

**Vertices**

| Call | Arguments | Returns / prints |
|---|---|---|
| `BuildVertices[ansatz]` | ansatz stage | Vertices stage; prints the vertex table |
| `VertexValue[vertices, {a1, …, an}, {{w1, q1}, …, {w(n-1), q(n-1)}}]` | fields; the n-1 independent (ω, q) legs | Γ⁽ⁿ⁾, last leg = −(sum) |
| `Gamma2Matrix[vertices, w, q]` | frequency, momentum symbol | Γ2_ab(w, q), field order, no regulator |
| `Gamma3Matrix[vertices, i, {nu, p}, {w, q}]` | field `i` (or its position), two (ω, q) legs | matrix Γ⁽³⁾_{i a b} |
| `ShowVertex[vertices, {a1, …, an}]`, `ShowGamma2[vertices]`, `ShowVertexTable[vertices]` | | print |

**Regulator**

| Call | Arguments | Returns / prints |
|---|---|---|
| `DefineRegulator[ansatz, …]` | `{a, b} -> R_ab, …`, `"Shapes" -> {r1, …}` | Regulator stage |
| `RegulatorMatrix[regulator, s]` | regulator stage, value of q² | R_ab(s), field order |
| `DsRegulatorMatrix[regulator, s]` | same | derived ∂s R_ab(s) (`dsCoupling[c]`, shape derivatives) |
| `DsRAbstractMatrix[regulator]` | | matrix of the symbols `dsR[a, b]` |
| `DsRRules[regulator, s]` | | rules `dsR[a, b] -> ∂s R_ab(s)` |
| `ShowRegulator[regulator]` | | prints R_ab and ∂s R_ab |
| `dsCoupling[c]`, `dsR[a, b]` | symbols appearing in results | k dc/dk; abstract entry of ∂s R |

**Propagator**

| Call | Arguments | Returns / prints |
|---|---|---|
| `BuildPropagator[vertices, regulator]` | both stages | Propagator stage |
| `PropagatorMatrix[prop, w, s]` | | G(w) in `GammaAtom` atoms and `Δ[±w, s]` |
| `FullMatrix[prop, w, s]` | | (Γ2 + R)(w) in `GammaAtom` atoms |
| `AtomRules[prop]` | | rules `GammaAtom[a, b, m] -> Function[s, value]` |
| `DeltaRules[prop]` | | rules for `Δ` and `DeltaCoefficient` in atoms |
| `ValidityConditions[prop]` | | inequalities for causal poles (a function of `s`) |
| `ShowPropagator[prop]`, `ShowPropagator[prop, {a, b}]` | | atoms, Δ, validity; or G_ab |
| `GammaAtom[a, b, m][s]`, `Δ[w, s]`, `DeltaCoefficient[m][s]` | symbols appearing in results | coefficient of (iω)^m; det of the (φ, φbar) block; coefficient c_m of Δ |

**Symmetries**

| Call | Arguments | Returns / prints |
|---|---|---|
| `DeclareSymmetries[ansatz, tags…, …]` | `Protected[…]`, `Neglected[…]`, `NonRenormalised[c, "Reason" -> …]`, `"TruncationOrders" -> <\|"Frequency" -> mf, "Momentum" -> mp\|>` | Symmetries stage |
| `ShowClassification[symmetries]` | | classification table of all monomials |
| `MonomialClasses` | | the class names (`"Declared"`, `"Protected"`, `"Neglected"`, …) |

**Loops**

| Call | Arguments | Returns / prints |
|---|---|---|
| `BuildLoopIntegrands[vertices, prop, regulator]` | three stages | Loops stage |
| `LoopIntegrand[loops, {a, b}, nu, p]` | field pair, external frequency and momentum | integrand of ∂s Γ2_ab (with `dsR`) |
| `LoopVariables[loops]` | | `{w, q}`, the loop variables |

**Projection**

| Call | Arguments | Returns / prints |
|---|---|---|
| `ProjectOnMonomials[loops, symmetries]` | | Projection stage |
| `ShowIntegrand[projected, coefficient]`, `ShowIntegrand[projected, {a, b}, {m, n}]` | a coupling combination, or pair and powers of (iω)^m (p²)^n | prints the ω-integrand |
| `LoopMomentumSquared` | | the symbol for q² in integrands and equations |

**FrequencyIntegration**

| Call | Arguments | Returns / prints |
|---|---|---|
| `IntegrateFrequencies[projected]` | | Integrated stage |
| `ShowFlow[integrated, {a, b}, {m, n}]` | pair, powers | integrated flow of a monomial |
| `FrequencyIntegral[n, a, b, deg]`, `FrequencyIntegralC[m]` | | ∫dω ωⁿ/(Δ(ω)^a Δ(−ω)^b) in the coefficients c_m |

**Closure**

| Call | Arguments | Returns / prints |
|---|---|---|
| `CheckClosure[integrated, symmetries]` | | Closure stage |
| `ShowClosureReport[closure]` | | class and outcome per monomial |

**Assemble**

| Call | Arguments | Returns / prints |
|---|---|---|
| `AssembleEquations[closure]` | | Equations stage |
| `ShowEquation[equations, comb]` | `comb` as written in the ansatz, e.g. `nk Ak` | prints ∂s(comb) = RHS, atoms, validity |
| `EquationOf[equations, comb]` | | the equation record (keys include `"Integrand"`) |
| `ScaleAtom[a, b, m, j]` | symbol in results | scale-stripped atom; values in `equations["Atoms"]` |

**Oracle**

| Call | Arguments | Returns / prints |
|---|---|---|
| `RunOracle[equations, opts]` | `"Dimensions"`, `"Tolerance" -> 10^-8`, `"Seed"`, `"ZeroFlows"` | Oracle stage; numerical cross-check |

**Stages**

| Call | Arguments | Returns / prints |
|---|---|---|
| `SaveStage[stage, file]` | any stage, `.wdx` path | the file name |
| `LoadStage[file, "Expect" -> {stages…}]` | the current theory, ansatz, regulator, symmetries | the stage; aborts if it was not built from them |
| `StageFingerprints[stage]` | | fingerprints of the stage and everything it embeds |

**Dimensions**

| Call | Arguments | Returns / prints |
|---|---|---|
| `DeriveDimensions[ansatz, …]` | `"References" -> {nk -> "nu", …}`, `"Dimensionless" -> {…}`, `"Names"`, `"Scale"` | Dimensions stage |
| `ShowDimensions[dims]` | | dimension table |
| `CheckDimensions[dims, <\|gk -> …, "Frequency" -> …\|>]` | hand-derived dimensions | compares exactly, aborts on mismatch |
| `DimensionOf[dims, x]` | coupling, field, `"Frequency"` or `"Time"` | k^x ref^y |
| `LinearPart[dims, c]` | coupling | canonical part of the dimensionless flow |
| `EmbeddedStage[stage, name]` | any stage, a stage name (`"Propagator"`) | the stage embedded at any depth, or `Missing` |
| `AnsatzOf[stage]` | any later stage | its Ansatz stage |

**Adimensional**

| Call | Arguments | Returns / prints |
|---|---|---|
| `AdimensionalEquations[equations, dims, "State" -> {…}]` | | Adimensional stage |
| `ShowBeta[adim, c]` | state coupling | linear part, loop integrand, atoms |
| `ShowEtaSystem[adim]` | | the η system row by row |
| `ShowCausalPoles[adim]` | | causal-pole conditions, minimal and dimensionless (`adim["CausalPoles"]`: `Compact`, `Explicit`, `Stripped`) |
| `AdimAtom[a, b, m, j]`, `LoopMeasure[d, aq]` | symbols appearing in results | dimensionless atom (values in `adim["Atoms"]`); angular measure |

**JuliaExport**

| Call | Arguments | Returns / prints |
|---|---|---|
| `ExportJulia[adim, "Dimension" -> d, "File" -> path]` | `"Form" -> "Compact"` (default) or `"Explicit"` | the path of the FlowRG file (FlowRG parts + `causal_poles`) |
| `JuliaParts[adim, d]` | | the parts of the file as Mathematica expressions (`"causalPoles"` included) |

**Core**: `FRGAssert[name, cond, details, hint]` raises a named check, `ExactlyZeroQ[expr]` tests exact zero at random rationals,
`PrintCard[stage]` reprints a stage card.

### Inspecting a loaded stage

A stage is an Association: `Keys[stage]` lists its contents, `stage["Stage"]` its name, `PrintCard[stage]` its card. A stage embeds
the stages it was built from, so one loaded `Equations` file gives access to the whole chain:

```
Equations ⊃ Closure ⊃ Integrated ⊃ Projection ⊃ Loops ⊃ { Vertices, Propagator, Regulator } ⊃ Ansatz ⊃ Theory     (Symmetries beside Closure)
Adimensional ⊃ { Equations, Dimensions }
```

`EmbeddedStage[stage, "Name"]` returns the embedded stage (or `Missing`); then use that module's accessors from the tables above.

| I want | Call (`eqs` = a loaded Equations stage) |
|---|---|
| the ansatz, its monomials | `ShowMonomials[AnsatzOf[eqs]]` |
| vertices, Γ2 | `v = EmbeddedStage[eqs, "Vertices"]`, `ShowVertex[v, {…}]`, `ShowGamma2[v]` |
| regulator | `ShowRegulator[EmbeddedStage[eqs, "Regulator"]]` |
| propagator, Δ | `prop = EmbeddedStage[eqs, "Propagator"]`, see below |
| symmetry classes | `ShowClassification[EmbeddedStage[eqs, "Symmetries"]]` |
| closure report | `ShowClosureReport[EmbeddedStage[eqs, "Closure"]]` |
| the flow of a coupling | `ShowEquation[eqs, comb]`, `EquationOf[eqs, comb]` |
| dimensions, β function (from an Adimensional stage) | `ShowDimensions[adim["Dimensions"]]`, `ShowBeta[adim, c]` |
| causal-pole conditions, dimensionless | `ShowCausalPoles[adim]` |

**Δ of the propagator matrix.** Δ(ω) = det of the (φ, φbar) block of Γ2 + R, whose roots are the poles of the retarded propagator.
It is stored in the Propagator stage as a function of the frequency and q², written in `GammaAtom` atoms:

```mathematica
eqs  = LoadStage["kpz.wdx", "Expect" -> {theory, ansatz, regulator, symmetries}];
prop = EmbeddedStage[eqs, "Propagator"];

prop["Delta"][w, s]                  (* I w + GammaAtom[h, hbar, 0][s]   for KPZ *)
prop["DeltaCoefficients"][s]         (* {c0, c1, c2}: Δ = c0 + c1 (I w) + c2 (I w)^2 *)
prop["DeltaDegree"]                  (* 1 or 2 *)
ValidityConditions[prop][s]          (* conditions for all roots in the upper half plane *)

(* with the couplings and regulator shapes written out *)
prop["Delta"][w, s] /. (g : GammaAtom[__])[x_] :> prop["AtomValues"][g][x]
                                     (* nuk s + I w + nuk s rn[s/k^2]   for KPZ *)
ShowPropagator[prop]                 (* prints atoms, Δ, coefficients, validity *)
```

**The propagator in the couplings.** `prop["G"][w, s]` is written in atoms and `Δ`; two substitutions make it explicit:

```mathematica
atoms = (g : GammaAtom[__])[x_] :> prop["AtomValues"][g][x];
G = prop["G"][w, s] /. DeltaRules[prop] /. atoms;      (* matrix in field order (h, hbar), no GammaAtom or Δ left *)
G[[1, 2]]                                              (* G_{h hbar}: 1/(nuk s - I w + nuk s rn[s/k^2]) for KPZ *)
G . (prop["M"][w, s] /. atoms) // Simplify             (* identity matrix: check *)
```

The same stage holds `prop["G"][w, s]` (propagator matrix, same as `PropagatorMatrix[prop, w, s]`) and `prop["M"][w, s]`
(Γ2 + R, `FullMatrix`). From a Loops stage, `loops["Propagator"]` gives it directly.

## Examples

`examples/ShallowWaterDE4plus.nb`, `examples/KPZ.nb` and `examples/KPZLO.nb` (the LO pipeline: functions of p², flows at
generic p, the Python module for BMWflow) are generated by `examples/makeNotebooks.wl`. Edit the cells there and regenerate. To check that a notebook runs without a front end:

```
math -script examples/runNotebook.wl examples/KPZ.nb                  # ~10 s
math -script examples/runNotebook.wl examples/KPZLO.nb                # ~10 s
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
- **Validity conditions.** The causal-pole conditions must hold at every loop momentum. The library does not check them: they
  are exported, minimal and dimensionless, as `causal_poles` in the equations file, and FlowRG does not call it yet
  (`docs/FUTURE.md`). Stripping their positive factors assumes the reference couplings are > 0.

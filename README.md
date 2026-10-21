# FRGscalarMSR

A Mathematica library that takes an ansatz for Γ_k of N scalar fields in the MSR (response-field) formalism and produces a
closed, checked set of Wetterich flow equations, using the derivative-expansion pipeline.
- The design is in `design/FRGDESIGN.md`.
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

## Workflow

Every step returns an immutable *stage* and prints a summary card: one PASS line per check, plus warnings and assumptions. A
failed check prints a red FAIL block and stops the evaluation; headless runs exit with a non-zero code.

```mathematica
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"]
theory     = DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]           (* Hbar, Psibar are created *)
ansatz     = DefineAnsatz[theory, "Density" -> Hbar (dt[H] - nk Ak lap[H] + ...) + ..., "Couplings" -> {nk, Ak, ...}]
vertices   = BuildVertices[ansatz]                    (* prints every Fourier vertex *)
regulator  = DefineRegulator[ansatz, {Hbar, H} -> nk q^2 arnh[q^2/k^2], ..., "Shapes" -> {arnh, ...}]
prop       = BuildPropagator[vertices, regulator]
symmetries = DeclareSymmetries[ansatz, Protected[...], Neglected[...], NonRenormalised[...], "TruncationOrders" -> <|...|>]
loops      = BuildLoopIntegrands[vertices, prop, regulator]
projected  = ProjectOnMonomials[loops, symmetries]
integrated = IntegrateFrequencies[projected]
flows      = CheckClosure[integrated, symmetries]     (* vanishing flows asserted, dropped terms reported *)
equations  = AssembleEquations[flows]                 (* compact and explicit forms *)
RunOracle[equations]                                  (* numerical cross-check in d = 2, 3 (d = 2 only if declared) *)
ShowEquation[equations, gk]   ShowClosureReport[flows]   SaveStage[equations, file]   LoadStage[file, "Expect" -> {...}]
```

The output equations are ∂s(coefficient as written in the ansatz) = ∫dᵈq/(2π)ᵈ (integrand). In each integrand:
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

### Dimensionless flows and the FlowRG file

The dimensions of every coupling are derived from the ansatz (Γ_k dimensionless), in the basis of k and N+1 reference couplings
you choose; then the equations are made dimensionless and written as a FlowRG equations file:

```mathematica
dims = DeriveDimensions[ansatz, "References" -> {nk -> "nu", Dhk -> "h", Dphik -> "phi"}, "Dimensionless" -> {Ak, aDphik}]
ShowDimensions[dims]                         (* coupling | monomials | dimension | name | linear part; fields, [w], z *)
CheckDimensions[dims, <|gk -> k Dphik^(1/2) Dhk^(-1/2) nk|>]   (* hand values, exact *)
adim = AdimensionalEquations[equations, dims, "State" -> {gk, lk, Hk, Ak, mk, Dmk, aDphik}]
ShowBeta[adim, gk]   ShowEtaSystem[adim]
ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{dir, "eqs", "Equations2D.jl"}]]
```

- ∂s ref = −η ref. A reference that appears alone in an equation (`nk`) gets its η from it; a reference that appears only times
  a dimensionless amplitude (`Dphik aDphik`) has its η fixed at run time (`MUST_FIX` in the file, FlowRG `fixed_etas`), and the
  equation becomes the flow of the amplitude. `"Mode" -> "Running"` restores the old convention (η from the amplitude equation).
- Regulator prefactors may contain reference couplings only (`nk q^2 arnh`, not `nk (1 + Ak) q^2 arnh`).
- From a saved stage: `DeriveDimensions[AnsatzOf[LoadStage[file]], …]`.
- KPZ reproduces Canet's Eq. (20) exactly after adimensionalisation (`tests/acceptance/KPZAdimensional.wl`); DE4plus agrees with
  the old hand-built `Equations2D.jl` up to its two known errors (`tests/acceptance/ShallowWaterDE4Adimensional.wl`).
- FlowRG extensions needed by some files (MUST_FIX guard, second shape derivatives `D2armh`) are specified in
  `docs/handoff/FlowRG-MUST_FIX.md`.

## Examples

`examples/ShallowWaterDE4plus.nb` and `examples/KPZ.nb` are generated by `examples/makeNotebooks.wl`. Edit the cells there and
regenerate. To check that a notebook runs without a front end:

```
math -script examples/runNotebook.wl examples/KPZ.nb                  # ~10 s
math -script examples/runNotebook.wl examples/ShallowWaterDE4plus.nb  # ~2.5 min
```

## Tests

```
math -script tests/runTests.wl [--verbose] [--require-julia]   # unit + acceptance; exit code 0 means every test passed
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

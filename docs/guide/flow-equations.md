# Flow equations, step by step

Part 1 of the derivative-expansion pipeline: from an ansatz for Γ_k to the dimensionful equations
∂s(coupling) = ∫dᵈq/(2π)ᵈ (integrand). The [README quickstart](../../README.md#quickstart-kpz-from-the-ansatz-to-the-flowrg-file)
runs the whole chain for KPZ; the two-field shallow-water case is in `examples/ShallowWaterDE4plus.nb`.

Every step returns an immutable *stage* (an Association) and prints a summary card: one PASS line per check, plus warnings and
assumptions. Full expressions are printed only on request, with the `Show…` commands. Each function is described in the
[API reference](../reference/README.md).

## Steps

1. **Fields.** `DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]` declares N = 1 or 2 physical scalars and creates the
   response field `<name>bar` of each (`Hbar`, `Psibar`). `Dim` is the space dimension, a symbol kept symbolic to the end (or
   `FixedDimension[2]`, see [below](#2d-levi-civita-term-jac-and-fixeddimension2)).
2. **Ansatz.** `DefineAnsatz[theory, "Density" -> …, "Couplings" -> {…}]`: Γ_k = ∫dt dᵈx (density). Write the density with the
   library's own operators `dt[…]`, `grad[…]`, `div[…]`, `lap[…]`, `jac[A, B]` and `grad[a].grad[b]` for a scalar product (not the
   System `Dt`, `Grad`: `Dt[nk H]` would differentiate the coupling). Couplings depend on k only and stand outside the operators
   (`nk dt[H]`, never `dt[nk H]`). Every monomial needs a response field; at most quartic in the fields.
   The parsed monomials are printed: check them.
3. **Vertices.** `BuildVertices[ansatz]` Fourier-transforms every monomial (∂t → +iω, ∇ → −iq, all momenta incoming) and prints
   every vertex under that convention, so you can check it by hand ([conventions](../CONVENTIONS.md)).
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
7. **Loops to equations.** `BuildLoopIntegrands` → `ProjectOnMonomials` (external frequency → 0, momentum projection on
   (p²)ⁿ) → `IntegrateFrequencies` (exact residues) → `CheckClosure` (protected flows vanish, dropped terms reported) →
   `AssembleEquations`.
8. **Inspect and check.** `ShowEquation[equations, comb]`, `EquationOf[equations, comb]`, `ShowMonomials[ansatz]`,
   `ShowVertex[vertices, {Hbar, H, Psi}]`, `ShowGamma2[vertices]`, `ShowRegulator[regulator]`, `ShowPropagator[prop, {H, Psibar}]`,
   `ShowClassification[symmetries]`, `ShowIntegrand[projected, nk Ak]`, `ShowClosureReport[flows]`. `RunOracle[equations]`
   recomputes every equation numerically (matrix inverse, numerical frequency integral, finite differences) in d = 2 and 3.
9. **Save.** `SaveStage[equations, file]`, `LoadStage[file, "Expect" -> {theory, ansatz, regulator, symmetries}]` (the
   fingerprints must match, else it aborts). See [inspecting stages](inspecting-stages.md).

## The output equations

The output equations are ∂s(coefficient as written in the ansatz) = ∫dᵈq/(2π)ᵈ (integrand), e.g. ∂s(nk Ak). In each integrand:

- `dsCoupling[c]` = k dc/dk enters linearly;
- the regulator shapes r(q²/k²) and their derivatives stay abstract;
- `Dim` stays symbolic, except with a declared d = 2 (next section).

## 2D Levi-Civita term: `jac` and `FixedDimension[2]`

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
- Details: [conventions](../CONVENTIONS.md), [modules](../MODULES.md) (Ansatz, Projection, Oracle) and the design record
  [design/LEVICIVITA.md](../../design/LEVICIVITA.md).

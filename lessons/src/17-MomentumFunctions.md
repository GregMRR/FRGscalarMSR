# Lesson 17 - Momentum functions

**Physics step.** Let a Γ2 coupling depend on the external momentum: ν_k → ν_k Fν(p²), D_k → D_k F_D(p²), the leading-order (LO) scheme of the KPZ literature. The flow of a function is a flow at **generic** external p: the loop carries propagators at two different momenta, the frequency integral has two pole families, and the angular integral stays in the equation. The decisions are recorded in `design/MOMENTUMFUNCTIONS.md`; this lesson follows the KPZ example through every stage.

**Stage.** The same chain as for KPZ (lessons 01 to 14), with `"MomentumFunctions"` declared in the ansatz. Output: functional equations in (q, p, c) and their dimensionless form; `ExportJulia` refuses them.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
theory = lessonUpTo["Theory"];
```

## 1. Writing the ansatz

A function is written bare, in a quadratic term, times its reference: `- nuk Fnu hbar lap[h]` means −ν_k Fν(p²) h̄ ∇²h. It is declared in `"MomentumFunctions"` only, not in `"Couplings"`.

```wl
loAnsatz = DefineAnsatz[theory, "Density" -> hbar (dt[h] - nuk Fnu lap[h] - lamk/2 grad[h].grad[h]) - Dk FD hbar^2,
   "Couplings" -> {nuk, lamk, Dk}, "MomentumFunctions" -> {Fnu, FD}];
```

Every input that makes no sense has a named check with a hint. A function under an operator, or in a cubic term:

```wl
showFailure[DefineAnsatz[theory, "Density" -> hbar dt[h] - nuk hbar lap[Fnu h] - Dk hbar^2, "Couplings" -> {nuk, Dk}, "MomentumFunctions" -> {Fnu}]]
```

## 2. Vertices and slots

In Γ2 the function becomes F(q²) at its leg momentum (both legs carry ±q):

```wl
loVertices = BuildVertices[loAnsatz];
Gamma2Matrix[loVertices, w, q]
```

A **slot** is a Γ2 entry at one power of iω. A slot with a function is Σ_{n<n₀} c_n p^{2n} + p^{2n₀} ref F(p²): its orders below n₀ are ordinary monomials, the rest is the function. `DeclareSymmetries` prints the orders ≥ n₀ as one row of class `Functional`. The KPZ tags keep their meaning: the shift symmetry still protects the p⁰ order of the Fν slot.

```wl
loSymmetries = DeclareSymmetries[loAnsatz,
   Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
   Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
   Neglected[_, "Reason" -> "beyond the LO ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
   "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>];
```

The card's ASSUMPTION `Symmetries.siblingAtZeroMomentum` says that the iω slot of Γ_{h̄h} is computed at p → 0 only: this is the LO approximation f^λ = 1, visible by default. `Neglected[{h, hbar}, "Frequency" -> 1, "Momentum" -> "AllMomenta", ...]` would compute it at generic p and report its size.

## 3. The loop at generic p

At generic p a constant-coupling ansatz gives a pointwise regular integrand, but with functions the quotient by p^{2n₀} is singular as p → 0 unless the integrand is even in p at fixed loop momentum. The library routes the loop at the midpoint (lines q − p/2 and q + p/2) and symmetrises under q → −q; the regulator derivative ∂sR sits on one line only, so the routing alone is not enough, and each `dsR` entry carries its line's momentum.

@quote Kernel/Loops.wl MidpointLoopIntegrand

```wl
loRegulator = DefineRegulator[loAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
loPropagator = BuildPropagator[loVertices, loRegulator];
loLoops = BuildLoopIntegrands[loVertices, loPropagator, loRegulator];
loProjected = ProjectOnMonomials[loLoops, loSymmetries];
```

The projection subtracts the orders below n₀ by their Taylor coefficients in p at fixed q (p → t p) and divides by p^{2n₀}; `Projection.functionalOrder` checks that the odd orders vanish, so the quotient is regular.

@quote Kernel/Projection.wl subtractLowerOrders

:::why Why the exact p = 0 row is computed separately
At ap = 0 the stored functional integrand is 0/0. Its value there (needed for η when F̂(0) = 1, and for a `Protected` tag at n₀) is the same integral as the derivative-expansion row of order n₀, but not the same integrand. Each functional row therefore also carries `"AtZero"`, built with the p → 0 operator of lesson 08.
:::

## 4. Several pole families

The generic-p integrand has Δ at |q − p/2|² and at |q + p/2|² (and at q² in the subtracted orders). `integrate` gives each distinct momentum argument its own family; a monomial with several families goes to `FrequencyIntegralFamilies`, which sums the residues of every family and cancels the factors that vanish where two families coincide (|q − p/2| = |q + p/2| happens inside the domain, at c = 0).

@quote Kernel/FrequencyIntegration.wl FrequencyIntegralFamilies

```wl
FrequencyIntegralFamilies[0, {{1, 1}, {1, 1}}, 1]
```

## 5. Equations in (q, p, c), and the oracle

```wl
loIntegrated = IntegrateFrequencies[loProjected];
loEquations = AssembleEquations[CheckClosure[loIntegrated, loSymmetries]];
Short[EquationOf[loEquations, FD]["Integrand"], 3]
```

`RunOracle` checks each functional equation numerically at an explicit (q, p, c): numeric inverse of Γ2 + R on both lines, ∂sR on its own line, the mirror q → −q, NIntegrate over ω, explicit smooth test functions; and its exact p = 0 row against the p → 0 path.

```wl
RunOracle[loEquations, "Dimensions" -> {2}];
```

## 6. Dimensionless form and Related

A function times a reference is dimensionless; its dimensionless version is named `a<>F` (F(p²) and F̂(p²/k²) are different functions). With `"Equation" -> nuk Fnu` the reference is Running and F̂ν(0) = 1.

```wl
loDims = DeriveDimensions[loAnsatz, "References" -> {nuk -> {"nu", "Equation" -> nuk Fnu}, Dk -> {"D", "Equation" -> Dk FD}}, "Scale" -> k];
loAdim = AdimensionalEquations[loEquations, loDims];
ShowBeta[loAdim, FD]
```

In d = 1 KPZ has the time-reversal symmetry: F̂_D = F̂ν. `Related` states it once; the library derives what it implies (η_D = η_ν, rd = rn, d = 1), checks that the two betas agree pointwise at c = ±1, imposes it and drops F_D:

```wl
loRel = AdimensionalEquations[loEquations, loDims, "Related" -> {Related[FD -> Fnu, "AtDimension" -> 1, "Reason" -> "time reversal (d = 1)"]}];
{loRel["State"], loRel["Dimension"]}
```

:::exercise Time reversal in d = 2
Ask for the same relation at d = 2. What happens, and why does the check need random cosines there?
:::

:::solution
```wl
showFailure[AdimensionalEquations[loEquations, loDims, "Related" -> {Related[FD -> Fnu, "AtDimension" -> 2, "Reason" -> "false in d = 2"]}]]
```
It aborts with `Adimensional.related`: the identity does not hold in d = 2. At c = ±1 (collinear p and q) the integrands reduce to the one-dimensional ones in every d, so a check at c = ±1 alone would accept it; at d ≠ 1 the library also draws random c in (−1, 1).
:::

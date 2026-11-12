# LO pipeline: couplings that are functions of p²

The second pipeline of the library. The Γ2 couplings become functions of the momentum, e.g. for KPZ
Γ^(1,1) = iω + ν p² f_ν(p²) and Γ^(0,2) = −2D f_D(p²), with f_ν(0) = f_D(0) = 1. The library derives the flow of each function
at generic p (one loop, two pole families: Δ at q and at Q = p + q), the anomalous dimensions from f(0) = 1, the dimensionless
form at fixed p/k, and writes a Python module for BMWflow.

The worked example is `examples/KPZLO.nb`; the exact comparison with the literature is the benchmark
[benchmarks/KPZLO](../../benchmarks/KPZLO/README.md). Design: [design/LOPIPELINE.md](../../design/LOPIPELINE.md); functions
at iω and in vertices (KPZ's f_λ): [design/LOFLAMBDA.md](../../design/LOFLAMBDA.md).

## 1. Ansatz with functions

`fn[f, X]` in the density is f(−Laplacian) acting on X: it multiplies the term by f(p²), p the momentum of X's leg.
`"Functions" -> {f -> reference}` states f(0) = 1 and names the k-dependent amplitude. The ansatz infers `"Pipeline" -> "LO"` from
the `fn`. The symmetries take the same tags as in the derivative expansion, but no truncation orders.

```mathematica
theory = DefineTheory["Scalars" -> {h}, "Dimension" -> Dim];
ansatz = DefineAnsatz[theory,
  "Density" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],
  "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}];
symmetries = DeclareSymmetries[ansatz,
  Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
  Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
  NonRenormalised[lamk, "Reason" -> "Galilean invariance"]];
```

## 2. Vertices, regulator, propagator, loops

The same stages as in the derivative expansion. The propagator atoms now contain `fnu[s]` and `fD[s]`.

```mathematica
vertices  = BuildVertices[ansatz];
regulator = DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
prop      = BuildPropagator[vertices, regulator];
loops     = BuildLoopIntegrands[vertices, prop, regulator];
```

## 3. Flows at generic p

`BuildLOEquations` replaces Projection … Assemble. For each function it takes the loop at external momentum p and zero external
frequency, integrates the loop frequency by residues at the poles of both families, and divides by the ansatz prefactor
(nuk p² for fnu, −2 Dk for fD). It checks that the entries outside the ansatz vanish at every p, the Ward identity, and that the
p⁰ order of the fnu row vanishes, and builds the η rows as the exact p → 0 limit. The integrand is explicit in p, q,
`CosTheta` = cos(p, q), f at q² and at Q² = p² + q² + 2 p q CosTheta, the shapes and `dsCoupling`.

```mathematica
equations = BuildLOEquations[loops, symmetries];
ShowLOEquation[equations, fD]
oracle = RunOracle[equations, "ZeroFlows" -> True];
```

## 4. Dimensionless flows

F_k(p²) = f(p²/k²), ap = p/k, aq = q/k. At fixed ap each function flows as ∂s f = η f + 2 ap² f′(ap²) + Loop. f(0) = 1 fixes
the η's. The constant coupling alamk = lamk Dk^(1/2) nuk^(−3/2) k^(d/2−1) has only a linear part.

```mathematica
dims = DeriveDimensions[ansatz, "References" -> {nuk -> "nu", Dk -> "D"}];
adim = AdimensionalEquations[equations, dims];
Thread[adim["Functions"] -> adim["Linear"]]
```

At f_D = f_ν = 1 the η system is the one of the derivative-expansion KPZ flow; `examples/KPZLO.nb` checks Canet's Eq. (20)
exactly at every d.

## 5. The BMWflow module

`ExportPython` writes the module of the contract [docs/handoff/BMWflow_LO.md](../handoff/BMWflow_LO.md): FUNCTIONS, ETAS,
COUPLINGS, P2_DIVIDED and the functions `linear`, `loop`, `eta_matrix`, `eta_vector`, `beta`. The loop measure is not in the
kernels (the BMWflow adapter applies it). `"G" -> alamk^2` writes every kernel in BMWflow's coupling g = alamk².

```mathematica
ExportPython[adim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{$TemporaryDirectory, "KPZLO", "KPZ_LO.py"}]]
```

## 6. NLO: functions of (ϖ, p)

The option `"ExternalFrequency" -> True` keeps the external frequency ϖ (`w`) in the loop: the functions on the left-hand side
depend on (p, ϖ), every function inside the loop stays at ϖ = 0 (Kloss, Canet, Wschebor 2012, Eq. 43). A function at (iϖ)⁰ flows
with the even part in ϖ of its Γ2 entry, a function at (iϖ)¹ (f_λ) with the odd part divided by iϖ. The η system and β do not
change. `AdimensionalEquations` adds (2 − η_ν) aw ∂f/∂aw to the linear part, and `ExportPython` writes a contract-v4 NLO module
for BMWflow's `ModelGenerated(..., "NLO", ...)` ([design/NLOPIPELINE.md](../../design/NLOPIPELINE.md), `benchmarks/KPZNLO`).

```mathematica
nloEquations = BuildLOEquations[loops, symmetries, "ExternalFrequency" -> True];
nloAdim = AdimensionalEquations[nloEquations, dims];
ExportPython[nloAdim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{$TemporaryDirectory, "KPZNLO", "KPZ_NLO.py"}]]
```

The limits of the LO pipeline are listed in [status](../STATUS.md#not-implemented-and-future-directions).

# FRGscalarMSR

A Mathematica library for the functional renormalization group (FRG) of N scalar fields in the MSR (response-field)
formalism. You write an ansatz for Γ_k; the library derives the one-loop Wetterich flow equations, checks them, makes them
dimensionless and exports them for numerical integration.

- **Derivative-expansion pipeline**: couplings that depend on k only → flow equations → dimensionless β functions → a
  FlowRG (Julia) equations file.
- **LO pipeline**: couplings that are functions of p² → flows at generic external momentum → a Python module for BMWflow.

Every step returns an immutable *stage* and prints a summary card with one PASS line per check. A failed check prints a red
FAIL block with a hint and stops the evaluation.

What is implemented, what is checked and what is planned: [docs/STATUS.md](docs/STATUS.md).

## Install

- Mathematica (Wolfram Language 14 or later). The tests use the `math` command-line kernel (Wolfram 15).
- Nothing to build. Clone the repository and load the library from a notebook or a script:
  `Get["/path/to/FRGscalarMSR/FRGscalarMSR.wl"]`.
- Optional: Julia and FlowRG to run the exported flows; Python 3 with NumPy and BMWflow for the LO modules.

## Quickstart: KPZ from the ansatz to the FlowRG file

The KPZ equation ∂t h = ν ∇²h + λ/2 (∇h)² + noise (variance 2D), about 5 s:

```mathematica
Get["/path/to/FRGscalarMSR/FRGscalarMSR.wl"]

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

## Documentation

All documentation is in [docs/](docs/README.md):

- Guides: [flow equations](docs/guide/flow-equations.md), [dimensionless flows and export](docs/guide/dimensionless-and-export.md),
  [LO pipeline](docs/guide/lo-pipeline.md), [inspecting stages](docs/guide/inspecting-stages.md),
  [development](docs/guide/development.md).
- [API reference](docs/reference/README.md), generated from the `::usage` strings of `Kernel/*.wl`.
- [Modules](docs/MODULES.md) (assumptions, manipulations and checks of each module), [conventions](docs/CONVENTIONS.md),
  [design records](design/).
- Worked examples: `examples/*.nb`. Teaching notebooks, one per module: [lessons/](lessons/README.md).
  Literature benchmarks: [benchmarks/](benchmarks/README.md).

## Tests

```
math -script tests/runTests.wl          # unit + acceptance, ~6 min; exit code 0 means every test passed
```

More commands: [docs/guide/development.md](docs/guide/development.md).

## License

GNU General Public License v3.0: see [LICENSE](LICENSE).

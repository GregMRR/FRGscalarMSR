# Benchmark: KPZ at LO (functions of p²)

**Reference.** The KPZ LO flow of BMWflow (`~/Physics/Sources/BMWflow`, commit `55108ad`): `notes/LO_flow_equations.tex` and
`notes/kpz_notation.tex` (git-ignored, local checkout only), coded by hand in `flow/models/model_kpz.py` (`Integral_pfixed_dD_LO`,
`eta_calc`, `g_rhs_calc`). Notation of Kloss, Canet, Wschebor, PRE 86, 051124 (2012), with f_λ = 1 and the functions at zero
frequency.

**Model**, symbolic d:

  Γ_k = ∫ h̄ { ∂t h − ν f_ν(−∇²) ∇²h − λ/2 (∇h)² } − D h̄ f_D(−∇²) h̄,  so Γ⁽¹¹⁾ = iω + ν p² f_ν(p²) and Γ⁽⁰²⁾ = −2 D f_D(p²),

with f_X(0) = 1. λ and the 1 in front of ∂t are protected (Galilean and shift symmetries). Regulators R_ν = ν q² r_n(q²/k²),
R_D = D r_d(q²/k²), two independent shapes; the notes use one shape r, substituted in the comparison.

**Why it matters.** First test of the LO pipeline (design/LOPIPELINE.md): functions of p² in the ansatz (`fn[f, X]`), the loop at
generic external momentum with two pole families (Δ at q and at Q = p + q), η from f(0) = 1 as the p → 0 limit of the same loop,
the flow at fixed p/k, and the Python export for BMWflow.

## Run
- `math -script benchmarks/KPZLO/run.wl`: a few seconds, exit 0 = all pass.
- Notebook: `math -script benchmarks/KPZLO/makeNotebook.wl` writes `KPZLO.nb`. Check it headless with
  `math -script examples/runNotebook.wl benchmarks/KPZLO/KPZLO.nb`.
- Suite: `math -script tests/runTests.wl --only KPZLO`.

## Files
- `Pipeline.wl`: the library stages, from `DefineTheory` to `AdimensionalEquations` (symbolic d).
- `Reference.wl`: the loop integrands of I_D and I_ν (kernels A_a, A_c, A_d), the η matrix and vector, the linear part and
  ∂s g, in the notes' variables ph, qh, ch = cos θ.
- `Checks.wl`: `kpzloLibrary[key]` rewrites the library's dimensionless stage in those variables (rn = rd = r, alamk² = g);
  `kpzloChecks[]` returns the records below.

## Results (all exact, pointwise in p, q, c, symbolic d)
Each comparison is tested with `ExactlyZeroQ`: random rationals for every symbol and for every atom f(x), f′(x), r(x), r′(x).

| Check | Reference | Result |
|---|---|---|
| `KPZLO.contract` | design/LOPIPELINE.md §6 | FUNCTIONS (fD, fnu), ETAS (etaD, etanu), COUPLINGS (alamk,), P2_DIVIDED (0, 1) |
| `KPZLO.loop.D` | I_D | library = notes |
| `KPZLO.loop.nu` | I_ν, with its 1/p² | library = notes |
| `KPZLO.eta.matrix` | I_DD, I_Dν, I_νD, I_νν | library = notes |
| `KPZLO.eta.vector` | I_D0, I_ν0 | library = notes |
| `KPZLO.linear` | η_X f_X + p f_X′ | library = notes |
| `KPZLO.coupling` | ∂s g = g(d − 2 − η_D + 3η_ν) | 2 alamk β_alamk = notes, g = alamk² |
| `KPZLO.symmetries` | Galilean, shift | Ward identity at p = 0, p⁰ order of Γ⁽¹¹⁾ = 0, Γ(h,h) = 0 |

The suite hook also checks that a 1/100 change of the A_d kernel of I_ν is detected. The numerical BMWflow acceptance (the
generated module run by BMWflow's `ModelGenerated` against `ModelKPZ`) is specified in `docs/handoff/BMWflow_LO.md`.

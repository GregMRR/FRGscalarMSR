# Benchmark: KPZ at LO with a flowing f_λ

**Model** (Kloss, Canet, Wschebor, PRE 86, 051124 (2012), LO, symbolic d; design/LOFLAMBDA.md):

  Γ_k = ∫ h̄ f_λ(−∇²) { ∂t h − λ/2 (∇h)² } − ν h̄ f_ν(−∇²) ∇²h − D h̄ f_D(−∇²) h̄,

so Γ⁽¹¹⁾ = iω f_λ(p²) + ν p² f_ν(p²), Γ⁽⁰²⁾ = −2 D f_D(p²), and the cubic vertex carries f_λ at the momentum of its h̄ leg
(Galilean invariance ties it to the f_λ of ∂t h). f_λ(0) = 1 is the Galilean Ward identity: f_λ has reference 1 and no η. Its
flow is ∂/∂(iν) of ∂sΓ⁽¹¹⁾(ν, p) at ν = 0, taken before the frequency integral (a double pole at Q = p + q).

**References.** No source gives the f_λ row in this form, so it is checked three ways:
- f_λ = 1 → the KPZLO benchmark (itself = BMWflow's LO notes), for f_D and f_ν;
- generic f_λ → the Kloss et al. kernels with f_λ slots kept in BMWflow `flow/models/model_kpz.py` (`Integral_pfixed_dD_LO`,
  set to 1 there), transcribed in `Reference.wl`, for f_D and f_ν;
- the derivative expansion of the same model at f_λ = 1 + x p², f_ν = f_D = 1: the DE flows of x, ν, D equal the p², p⁰ Taylor
  coefficients of the f_λ, f_ν, f_D rows (an independent route for the sign and normalisation of the ∂/∂(iν) row).

## Run
- `math -script benchmarks/KPZLOLambda/run.wl`: ~15 s, exit 0 = all pass. Suite: `math -script tests/runTests.wl --only KPZLOLambda`.
- Notebook: `math -script benchmarks/KPZLOLambda/makeNotebook.wl` writes `KPZLOLambda.nb`.
- BMWflow (contract v3, after BMWflow implements it; docs/handoff/BMWflow_LO.md): `math -script benchmarks/KPZLOLambda/export.wl`,
  then `~/sources/BMWflow/.venv/bin/python benchmarks/KPZLOLambda/runBMWflow.py` → z, χ in d = 1, 2, 3 with and without
  f_λ, plots in `out/` (git-ignored). In d = 1, f_λ ≠ 1 breaks the time-reversal symmetry of KPZ (∫ h″ f_λ(−∂²)(h′)² ≠ 0), so
  the exact d = 1 exponents are not protected in the f_λ model.

## Files
- `Pipeline.wl`: the LO stages (kpzlam*) and the DE stages at f_λ = 1 + x p² (kpzlamDE*).
- `Reference.wl`: the f_D, f_ν loop integrands with f_λ at p, q, Q.
- `Checks.wl`: `kpzlamChecks[]` (needs the KPZLO stages as well).
- `export.wl`, `runBMWflow.py`: acceptance 2.

## Results (exact, ExactlyZeroQ, symbolic d)
| Check | Result |
|---|---|
| `KPZLOLambda.contract` | FUNCTIONS (fD, fnu, flam), ETAS (etaD, etanu, None), P2_DIVIDED (0, 1, 0) |
| `KPZLOLambda.flam1.*` | f_λ = 1: integrands, η system, linear parts = KPZLO's (Equations and Adimensional stages) |
| `KPZLOLambda.coupling` | β of λ unchanged |
| `KPZLOLambda.reference.D`, `.nu` | generic f_λ: library = Kloss et al. kernels |
| `KPZLOLambda.de.flam`, `.de.nuD` | p² of the f_λ row = DE flow of x (sign flip rejected); p⁰ of f_ν, f_D rows = DE flows of ν, D |
| `KPZLOLambda.flam.linear`, `.flam.flows` | ∂s f_λ = p f_λ′ + I_λ, I_λ ≢ 0 |
| `KPZLOLambda.symmetries` | Ward check, unit reference, η rows hold f_λ(0) = 1 only |
| `KPZLOLambda.oracle` | residues (double poles) = numerical ω integral with exact ∂ν, d = 2, 3 |

The suite hook also checks that a 1 % change of f_λ(p) in the f_ν row is detected, the Python round trip of the module, and
that the module at f_λ ≡ 1 equals `KPZ_LO.py` (numpy, 1e-13).

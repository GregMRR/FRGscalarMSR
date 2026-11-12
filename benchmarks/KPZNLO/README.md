# Benchmark: KPZ at NLO (functions of p and of the external frequency)

**Model** (Kloss, Canet, Wschebor, PRE 86, 051124 (2012), arXiv:1209.4650, "KCW", Eqs. 43, 44, 49; design/NLOPIPELINE.md):
the ansatz of `KPZLOLambda` (and of `KPZLO` with f_λ ≡ 1),

  Γ_k = ∫ h̄ f_λ(−∇²) { ∂t h − λ/2 (∇h)² } − ν h̄ f_ν(−∇²) ∇²h − D h̄ f_D(−∇²) h̄,

with f_D, f_ν, f_λ functions of (ϖ, p) on the left-hand side. In the loop every function is at zero frequency (KCW Eq. 43); ϖ
enters only through G(ω + ϖ, Q). Rows: f_D = even part in ϖ of ∂sΓ⁽⁰²⁾/(−2D); f_ν = even part of ∂sΓ⁽¹¹⁾/(ν p²); f_λ = odd part of
∂sΓ⁽¹¹⁾ divided by iϖ (exact, finite at ϖ = 0). Linear part: KCW Eq. 56, with 2 − η_ν derived from the dimension of the frequency.

Names: KCW "LO" (f_λ ≡ 1) is BMWflow `"NLO"`; KCW "NLO" (flowing f_λ) is the `kpzlamnlo*` stage here.

**References.**
- f_D, f_ν rows at generic f_λ and ϖ: BMWflow's NLO kernels with f_λ slots (`flow/models/model_kpz.py`,
  `Integral_pfixed_dD_NLO`; `notes/NLO_flow_equations.tex`), transcribed in `Reference.wl`.
- ϖ = 0: the `KPZLO` and `KPZLOLambda` benchmarks (the f_λ row there is checked against the derivative expansion).
- f_λ row at ϖ ≠ 0: KCW Eq. 39 (zero at p = 0) and the numerical ω integral (Oracle). No published form exists.

## Run
- `math -script benchmarks/KPZNLO/run.wl`: ~45 s, exit 0 = all pass. Suite: `math -script tests/runTests.wl --only KPZNLO`.
- Notebook: `math -script benchmarks/KPZNLO/makeNotebook.wl` writes `KPZNLO.nb`.
- BMWflow (contract v4): `math -script benchmarks/KPZNLO/export.wl` writes `out/KPZ_NLO.py`, `out/KPZLambda_NLO.py` and
  `out/KPZLambda_NLO_points.json` (git-ignored); BMWflow vendors them as `tests/fixtures/`.

## Files
- `Pipeline.wl`: the NLO stages `kpznlo*` (f_λ ≡ 1) and `kpzlamnlo*` (flowing f_λ); needs the KPZLO and KPZLOLambda pipelines.
- `Reference.wl`: the f_D, f_ν NLO loop integrands with f_λ at p, q, Q and ϖ.
- `Checks.wl`: `kpznloChecks[]`.
- `export.wl`: the two modules and the reference points for BMWflow.

## Results (exact, ExactlyZeroQ, symbolic d)
| Check | Result |
|---|---|
| `KPZNLO.contract` | NLO stages; FUNCTIONS, ETAS, P2_DIVIDED as LO |
| `KPZNLO.reference.D`, `.nu` | generic f_λ and ϖ: library = BMWflow's NLO kernels |
| `KPZNLO.flam1` | f_λ = 1: f_D, f_ν rows = the KPZ NLO rows |
| `KPZNLO.slice` | ϖ = 0: loop rows = LO rows; η system and β unchanged |
| `KPZNLO.linear` | ∂s f_X = η_X f_X + p ∂p f_X + (2 − η_ν) ϖ ∂ϖ f_X (KCW Eq. 56) |
| `KPZNLO.flam.ward`, `.flam.flows` | f_λ row = 0 at p = 0 for every ϖ (KCW Eq. 39); it depends on ϖ |
| `KPZNLO.symmetries` | the NLO checks of `BuildLOEquations` pass |
| `KPZNLO.oracle` | even/odd parts at ϖ = 0.3, 1.4, 4 = numerical ω integral, d = 2, 3 |

The suite hook also runs the Python round trip of both modules (including ϖ = 0) and checks that a perturbed row is rejected.

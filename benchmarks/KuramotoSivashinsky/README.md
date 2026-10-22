# Benchmark: Kuramoto–Sivashinsky at order p⁴

**Reference.** L. Gosteva, N. Wschebor, L. Canet, *Scaling regimes of the Kuramoto–Sivashinsky equation from the functional
renormalization group*, arXiv:2607.15784 (2026), Sec. III ("simplest DE ansatz"). The authors' notebook `DE/supmat.nb`
(gitlab.com/goste/kuramoto-sivashinsky-lpa, GPL-3) was used only to check the transcription. It is not copied here.

**Model** (Eq. 16), d = 1:

  Γ_k = ∫ h̄ { ∂t h − (ν − τ ∂x²) ∂x² h − λ/2 (∂x h)² − (D − ζ ∂x²) h̄ },  so Γ⁽¹¹⁾ = iω + p²(ν + τp²) and Γ⁽⁰²⁾ = −2(D + ζp²).

λ and the 1 in front of ∂t are protected by the Galilean and shift symmetries. The regulators are R^ν = D r_ν(q²/k²) and
R^D = D r_D(q²/k²) (Eq. 31), with two independent generic shapes.

**Why it matters.** ∂sτ is the p⁴ coefficient of ∂sΓ⁽¹¹⁾, so it is the first literature check of the p^{2n} projection beyond
n = 1 (design/MOMENTUMPROJECTION.md).

## Run
- `math -script benchmarks/KuramotoSivashinsky/run.wl`: about 25 s, exit 0 = all pass.
- Notebook: `math -script benchmarks/KuramotoSivashinsky/makeNotebook.wl` writes `KuramotoSivashinsky.nb`. Check it headless with
  `math -script examples/runNotebook.wl benchmarks/KuramotoSivashinsky/KuramotoSivashinsky.nb`.
- Suite: `math -script tests/runTests.wl --only KuramotoSivashinsky`.

## Files
- `Pipeline.wl`: the library stages, from `DefineTheory` to `AssembleEquations` (symbolic d).
- `Reference.wl`: Eqs. (34)–(37) as integrands of ∫dq̂/2π over ℝ, in the paper's variables Ĝ = λ²/(kD²), F̂ = ν/D,
  Ĉ = k²ζ/D, T̂ = k²τ/D and η = −∂sD/D.
- `Checks.wl`: `ksLibrary[c, d]` rewrites the library's integrand in those variables. `ksChecks[]` returns the records below.

## Results (all exact, pointwise in q̂)
Each check is tested with `ExactlyZeroQ`: random rationals for every symbol and for every shape atom r, r′, r″, r‴ independently.

| Check | Paper | Result |
|---|---|---|
| `KS.protected` | Sec. III A | the iω and p⁰ monomials of Γ⁽¹¹⁾ flow to 0; λ is declared non-renormalised |
| `KS.equation.D` | Eq. (34) | library = paper |
| `KS.equation.nu` | Eq. (35) | library = paper |
| `KS.equation.zeta` | Eq. (36) | library = paper |
| `KS.equation.tau` | Eq. (37), Γ⁽¹¹⁾ at p⁴ | library = paper |
| `KS.thirdDerivatives` | below Eq. (41) | r‴ appears for general shapes and cancels for r_ν = α r_D |
| `KS.timeReversal.nu` | Eqs. (19), (32) | d = 1, r_ν = α r_D, ν = αD, τ = αζ: ∂sν = α ∂sD |
| `KS.timeReversal.tau` | Eqs. (19), (32) | at the same point, ∂sτ = α ∂sζ (the p⁴ projection against the p² projection) |
| `KS.timeReversal.controlD2` | control | both identities fail at d = 2 |

The suite hook also checks that a 1 % change in one coefficient of Eq. (37) is detected.

## Notes and limits
- **Normalisations.**
  - The paper writes s_ν = ∂sR^ν/ν. The comparison needs s_ν = ∂sR^ν/D, which is what the paper's own
    s_ν = α s_D = −α(η r + 2x r′) implies.
  - D_Λ in Ĝ = λ²D_Λ/(κD²) is a unit, set to 1.
- **The paper's notation.** Eqs. (36)–(37) are written with q̂² in place of x = q̂². The primes are d/dx.
- **Normalisation not tested.** At d = 1 the n = 2 normalisation 8d(d+2) equals 8(d²+2) = 24, so a d-normalisation error would
  pass this benchmark. It stays covered by `Projection.normalisation` and `Projection.explicitVectors` (docs/FUTURE.md,
  "Generic-d p⁴ acceptance").
- **Out of scope (symbolic only).**
  - The sharp-cutoff limits (Eqs. 40–41).
  - The modified Litim regulator (Eq. 42): it needs r″, and `litimSmoothGuard` rejects it.
  - The numerical fixed point (η* = 1/2, T̂* ≈ Ĉ* ≈ −0.003). The exported p⁴ flows need r‴, which FlowRG does not bind yet.

# LOPIPELINE — couplings that are functions of p² (BMW-type LO), exported for BMWflow

Decision record and implementation spec, written for the agent that implements it. Grilled with the user, then attacked
seven times: `design/reviews/LOPipeline-r1-r7.md` (arch-enemy r1–r6, the-enemy r7; all blockers and majors resolved, §8).
Line anchors refer to commit `21719c9` (FRGscalarMSR, branch BMW-LO) and `55108ad` (BMWflow).

## 0. For the implementing agent
- **Read first:** `CLAUDE.md` (map, conventions, gotchas), `docs/MODULES.md`, `design/FRGDESIGN.md` §7,
  `design/CAUSALITY.md` (style of a decision record), `tests/reference/KPZSetup.wl`, `tests/acceptance/KPZAdimensional.wl`,
  and in BMWflow: `notes/LO_flow_equations.tex`, `notes/kpz_notation.tex`, `flow/models/model_base.py`,
  `flow/models/model_kpz.py`. The BMWflow `notes/` folder is git-ignored: it exists only in the local checkout.
- **Where to work:** a git worktree branched from `BMW-LO`. Other sessions may be editing the main checkout.
- **Scope:** the Mathematica side only (steps 0, 2–12 in §5). The BMWflow side (`ModelGenerated`) is done by another agent
  from the contract in §6, written to `docs/handoff/BMWflow_LO.md`. Do not edit BMWflow.
- **Stop and ask the user** if: the step-0 baseline changes; a decision not covered here is needed; a check of §5 cannot be
  made to pass without weakening it; the KPZ benchmark disagrees with `LO_flow_equations.tex`.
- **Rules that hold throughout** (from `CLAUDE.md`): every check is a named `FRGAssert`/`assertCheck` `Module.check` with an
  `expectFail` test; exact zero tests use `ExactlyZeroQ`; no `Simplify` on large expressions; stages are immutable and
  validated; `Dim` stays symbolic; helper names must be distinctive (all modules share ``FRGscalarMSR`Private` ``).
- Step 1 of the plan (this record + its reviews) is done. Start at step 0.

## 1. Goal
A second pipeline next to the derivative expansion (DE). The Γ2 couplings are a reference coupling (k-dependent, for the
scaling analysis) times a dimensionless function of p² (isotropic). The vertices are constants. The flow of each function
comes from the one-loop Wetterich equation for the 1PI 2-point function at zero external frequency. At ν = 0 the loop has a
finite set of poles (from the q line and the Q = p + q line, not shifted), so the internal ω integral is done by residues.
The pipeline checks causality, the declared symmetries and homogeneity, goes dimensionless, and exports a Python module
that BMWflow integrates through a generic adapter.

Benchmark: KPZ LO as written in `BMWflow/notes/LO_flow_equations.tex` + `kpz_notation.tex` and coded by hand in
`BMWflow/flow/models/model_kpz.py` (LO: `Integral_pfixed_dD_LO` :52, `eta_calc` :180, `g_rhs_calc` :236).

Constraints: the DE pipeline keeps working unchanged (step 0 guards it). Choices that make a future **mixed** ansatz (DE
monomials + functions) easier are preferred when they are as reliable as the alternatives.

## 2. Facts the design rests on
- **BMWflow** is Python. No equation-file format: a model is a `ModelBase` subclass (`flow/models/model_base.py:9`),
  registered in `flow/create.py` `MODELS`. One RG step (`flow/evolution.py:190-197`): `f_spline_upd` → `calc_upd` (f at
  q and at Q = |p+q| from cubic splines) → `Integral_upd(g, eta)` → `eta_calc(g)` → `g -= ds*g_rhs_calc` →
  `f -= ds*f_rhs_calc`; base `f_rhs_calc = eta f + p f′ + Integral` (model_base.py:366-370). Quadrature: Gauss–Legendre in
  (q, θ); measure v_{d−1} q^{d−1} sin^{d−2}θ, d = 1 special case θ ∈ {0, π}, v_{d−1} → 1 (model_base.py:184-230). Regulator:
  callables `r(q)`, `r_(q) = dr/d(q²)`, one shape for both sectors. Evolution state: f (n_f, Np, Nw), one scalar g,
  eta (n_f,). The KPZ p=0 extrapolation (model_kpz.py:243-263) is KPZ-only, row 1, and divides by p² itself.
- **FRGscalarMSR, reusable as is:** Theory, Regulator, Propagator (`GammaAtom[a,b,m][s]` is already a function of s = q²,
  Propagator.wl:37-49), Loops (`LoopIntegrand`, Loops.wl:24-29: ∂sR on the q line, G at p+q, ν and p symbolic), Core
  (`ExactlyZeroQ` randomises function atoms), Stages, Projection's `directionalAverage` (Projection.wl:22-35).
- **DE-only today:** Ansatz parser (k-only couplings; `"Pipeline"` option, Ansatz.wl:136, asserted DE at :146), Vertices
  `gamma2Polynomial`, Symmetries (coefficient extraction in p²), Projection, Closure, Assemble, Dimensions (one equation per
  monomial), Adimensional `homVec` (Adimensional.wl:31-45 rejects units inside function arguments).
- **Bug found during review:** the DE frequency dispatch (`integrate`, FrequencyIntegration.wl:49-50) maps every
  `Δ[w, _]` to the same pole factor, so a loop with Δ(ω,q) and Δ(ω,Q) is integrated silently wrong
  (`1/(Δ[w,q²]Δ[-w,Q²])` → π/l_q instead of 2π/(l_q+l_Q)). `docs/FUTURE.md:6-8` and `FRGDESIGN.md:193` wrongly claim that
  shifted poles "plug in without changes".
- **Fingerprints are not reproducible** across runs: `contentHash` (Core.wl:102-104) drops volatile keys only at the top
  level, so nested `Timing`/`Checks`/`Seconds` enter the hash. Do not baseline fingerprints.
- **KPZ needs no second derivatives:** the p→0 limit of the KPZ loop (the η system) contains only f, f′, r, r′ (verified
  by three independent probes, r2, r3, r7).

## 3. Decisions
| # | Decision | Reason |
|---|---|---|
| L1 | Functions only in bilinear (Γ2) terms. Γ3/Γ4 couplings are k-only constants (reference or NonRenormalised). | KPZ LO; momentum-dependent vertices would need Γ3 at three momenta. |
| L2 | Γ2 polynomial in iω; iω coefficients constant and Protected; only ω⁰ functions flow. The parser accepts function × iω, a check rejects it (`Ansatz.functionTimesOmega`). | ν = 0 fixes only ω⁰ coefficients; ∂/∂(iν) at ν = 0 is the later extension. |
| L3 | Syntax `fn[f, X]` = f(−∇²) acting on X (one field, optionally under grad/lap). `"Functions" -> {fD -> Dk, fnu -> nuk}` ties f to its reference and states f(0) = 1; the declared order is the export's `FUNCTIONS` order. | In Fourier space f takes the momentum of X; in a bilinear term that is ±p, so `fn[fnu, lap[h]]` ≡ `lap[fn[fnu, h]]`. For future vertex functions the leg matters, and the syntax already names it. |
| L4 | Coupling kind (constant / function) is per coupling; the pipeline is inferred (`"Pipeline" -> Automatic`); an explicit value must match. A mixed ansatz aborts (`Ansatz.mixedNotYet`). | Mixed ansätze later only lift one check. |
| L5 | η: each reference Running via f(0) = 1; the η system is derived automatically as the exact p→0 limit of the same loop with `directionalAverage`. Scope = KPZ-like models: η rows may contain only f, f′, r, r′, else abort `ExternalMomentum.firstDerivativesOnly`. No integration by parts. | Reuses the DE projection (the seam for mixed ansätze). Integrating f″ by parts brings r″ through ∂sR (r1 E3), and BMWflow has no r″; the user restricted scope instead. |
| L6 | p = 0 grid point: BMWflow-style 'const'/'quad' extrapolation, generalised in the adapter (§6). Exact limit / small-p series → FUTURE.md. | User: stay with how BMWflow works. |
| L7 | Symmetries: structural/causal zero entries checked at every p (ExactlyZeroQ after the ω integral); Protected iω coefficient checked at p = 0 (∂/∂(iν), Taylor; Oracle fallback), carded "Neglected (LO)" at p ≠ 0; Protected low-p entries checked by `dividedOrders` (§5 step 5); NonRenormalised = declaration; relations between functions (d = 1 time reversal) → FUTURE.md. | LO drops the frequency dependence of the functions, so the Galilean Ward identity holds only at p = 0. |
| L8 | Causality: the causal-pole conditions are used only as residue assumptions, not exported (export → FUTURE.md). New check `FrequencyIntegration.regularAtCoincidence`. | User: the current system does not need the export. |
| L9 | Export like DE: linear part (derived; KPZ η_X f + p f′) and loop part, η system linear in η, β of the dimensionless constants (KPZ: `alamk`). BMWflow's g = alamk² is mapped by the adapter. The adapter overrides `f_rhs_calc` with the exported linear part. | Same l/nl split as `betal`/`betanl`; a future ansatz with another scaling part needs no BMWflow change. User chose the adapter mapping. |
| L10 | Numerical Oracle for the LO ω integral. | Independent check of the new two-family residue engine. |
| L11 | Shared modules extended additively; new modules where the method differs; no DE refactor now. | DE stays byte-for-byte; the generic-p flow object is the future seam (DE = its Taylor expansion). |
| L12 | Benchmark symbolic at generic d (this repo); BMWflow acceptance numerical in d = 1 and 2 (other agent). | |
| L13 | Docs now: this record + reviews, MODULES.md, tests, benchmark + notebook, FUTURE.md, CLAUDE.md map. Lessons later, except the one stale line in lesson 01. | User. |

## 4. Pipeline
```
DefineTheory → DefineAnsatz (fn, "Functions"; Pipeline → "LO")
  → DeclareSymmetries (LO: declarations only) → BuildVertices → DefineRegulator → BuildPropagator → BuildLoopIntegrands
  → BuildLOEquations   (new, Kernel/ExternalMomentum.wl; replaces ProjectOnMonomials, IntegrateFrequencies,
                        CheckClosure, AssembleEquations; output "Stage" -> "Equations", "Pipeline" -> "LO")
  → RunOracle (LO branch)
  → DeriveDimensions (LO branch) → AdimensionalEquations (LO branch)
  → ExportPython (new, Kernel/PythonExport.wl) → eqs/<Model>_LO.py
```
KPZ ansatz (names as in `KPZSetup.wl`):
```
DefineAnsatz[kpzTheory,
  "Density" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],
  "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}]
```
In Fourier space: Γ^(1,1) ∋ iω + nuk p² fnu(p²), Γ^(0,2) ∋ −2 Dk fD(p²).

## 5. Implementation steps
Order = `$FRGModules` order. Each step lands with unit tests and an `expectFail` per new check. Run the DE suite after
steps 2–4.

0. **DE baseline** (`tests/reference/DEBaseline.wl` + test `tests/acceptance/DEBaseline.wl`). Record:
   (a) SHA-256 of the exported `.jl` text without its two header lines (date, fingerprint; JuliaExport.wl:108-109) for KPZ
   (`.jl` at d = 2 and 3) and ShallowWater DE4 (setup of `ShallowWaterDE4Adimensional.wl`; `.jl` at d = 2). KS has no Julia
   export (its flows contain r‴, which JuliaExport rejects).
   (b) a hash of the DE `Equations` stage (KPZ, ShallowWater DE4, KS) and `Adimensional` stage (KPZ, ShallowWater DE4) with
   `Timing`, `Checks`, `Fingerprint`, `Warnings`, `Upstream`, `Seconds` removed at every nesting level (reproducible across
   kernels and histories: r5/r6 probes; reference values from r5/r6: KPZ Equations `2306608a…`, KS Equations `3a053cfb…`).
   The test rebuilds what it needs from `tests/reference/*Setup.wl` and the benchmark pipeline (no reliance on other test
   files' globals). First verify that two fresh kernels give the same baseline. The test stays in the suite after the task
   as a DE regression guard. DE stages gain no new keys: Ansatz adds `"Functions"`/kinds keys only when the pipeline is LO.
1. **Design record + review** — done (this file, `design/reviews/LOPipeline-r1-r7.md`).
2. **Ansatz** (`Kernel/Ansatz.wl`): `fn` head, `"Functions"` option; `"Pipeline"` option stays with default Automatic =
   inferred; explicit value must match (`Ansatz.pipeline`; the existing `"MomentumDependent"` expectFail at
   AnsatzTest.wl:64 must still fail). Coupling-kinds key only in LO stages. Checks: `Ansatz.fnBilinear` (L1),
   `Ansatz.fnOneField`, `Ansatz.functionReference`, `Ansatz.functionTimesOmega`, `Ansatz.mixedNotYet`.
3. **Pipeline tag propagation and guards.** Vertices, Regulator, Propagator, Loops, Symmetries, Dimensions
   (Dimensions.wl:189), Adimensional (:264), Oracle (Oracle.wl:122) stamp `"Pipeline"` from upstream instead of the literal
   "DerivativeExpansion". Theory (Theory.wl:38) is pipeline-agnostic and keeps its stamp. Guards with expectFail on an LO
   stage: `JuliaExport.pipeline`, `Projection.pipeline`, `Closure.pipeline`, `Assemble.pipeline` (DE only),
   `PythonExport.pipeline` (LO only). Vertices: `gamma2Polynomial` only for DE; LO check that Γ2 coefficients are
   ref × F(s) or constants. Propagator: `GammaAtom[a,b,m][s]` bound to the functions. Symmetries: LO branch records the
   declarations only (no `Classification` by p² coefficients).
4. **FrequencyIntegration.** Generalise `FrequencyIntegral` (FrequencyIntegration.wl:24-46) to several pole families,
   Δ(ω,q) and Δ(ω,Q) (external ν = 0), deg Δ ≤ 2 per family. The `evenInBeta`/`rational` checks apply per family (the
   residue sum is symmetric in each family's ±β separately). New check `regularAtCoincidence`: on the `Together`-simplified
   sum of residues (single residues carry cancelling (z_q − z_Q) factors), no denominator vanishes when the Q data equal the
   q data. DE dispatch (`integrate`, :49) unchanged except a new abort `FrequencyIntegration.singleMomentum` on
   mixed-momentum Δ (expectFail: the example in §2).
5. **New `Kernel/ExternalMomentum.wl`**, `BuildLOEquations[loops, symmetries]`. From `LoopIntegrand` keep ν and generic p.
   - Ward check at p = 0: ∂/∂(iν), then `directionalAverage`, then ν → 0 (order of Projection.wl:85-86).
   - Flows: ν = 0; `scal` → p, q, c (c = cos θ); for each flowing Γ2 coefficient: residue integral (step 4), divide by the
     ansatz prefactor (KPZ: nuk p² for fnu, −2 Dk for fD), record the power of p² divided out (`P2Divided`).
   - Zero entries (structural/causal) checked at every p with `ExactlyZeroQ`.
   - `ExternalMomentum.dividedOrders`: for a row divided by p^{2n}, `directionalAverage` of every order < n is
     `ExactlyZeroQ`. This is also how the declared Protected ω⁰ low-p entries are checked (as Closure.wl:38 does for DE).
     expectFail: KPZ + a constant Γ4 term `c4 hbar h grad[h].grad[h]` (its tadpole gives p⁰ ∝ c4∫q²C(q) ≠ 0).
   - `ExternalMomentum.isotropicAngles`: after `scal` → p, q, c the integrand depends on the angle only through c (rejects
     `jac`/wedge terms; `FixedDimension` LO → FUTURE.md).
   - η system: p → 0 limit of each row with `directionalAverage`, linear in η (∂sR contains −η r).
     `ExternalMomentum.firstDerivativesOnly` (L5); expectFail: a constant-vertex bubble whose η row needs r″.
   - Output `"Stage" -> "Equations"`, `"Pipeline" -> "LO"`, so `AdimensionalEquations` and `RunOracle`
     (`validateStage[…, "Equations"]`, Adimensional.wl:98, Oracle.wl:53) accept it and branch on `"Pipeline"`.
6. **Oracle** (`Kernel/Oracle.wl`), LO branch: at random rational p, q, c and random positive function values, the residue
   kernel equals direct numerical ω quadrature of the un-integrated loop.
7. **Dimensions / Adimensional**, LO branch.
   - Dimensions: functions are dimensionless; references and constants use the exponent-vector engine. Check
     `Dimensions.functionReferences`: each function's reference is a declared Running reference of `DeriveDimensions`.
   - Adimensional: relabel F_k(p²) = f(p²/k²) = f(ap²), so no unit sits inside a function argument (`homVec` rejects
     F(k² ap²)). The flow is taken at fixed p̂ = p/κ (∂s|p̂ = ∂s|p + p∂p), which gives the +p f′ term. Keep `homVec`,
     `unitInvariance`, `noBeta` (:184) and `unusedEquation` (:168) for the LO branch.
8. **New `Kernel/PythonExport.wl`**, `ExportPython[adim, "Model" -> "KPZ", "File" -> …]`. Own printer (not `toJulia`,
   JuliaExport.wl:27-31: `^` is XOR in Python): Power → `**`, rationals as float literals, `math.pi`; pure numpy,
   broadcasting, no value-dependent branching. Converts internal F′(s) = d/d(p²) to f′/(2q) (f′ = d/d|momentum|); r_ stays
   d/d(q²). Checks: `PythonExport.etaPerFunction` (n_η = n_f; ETAS in FUNCTIONS order), `PythonExport.singleCoupling`
   (n_c = 1), `PythonExport.p2DividedAtMostOne` (for ≥ 2 the θ-quadrature error divided by p₁⁴ breaks the tolerance;
   → FUTURE.md), `PythonExport.noMeasure` (no v_d, v_{d−1}, Γ(d) in the kernels), `PythonExport.modelTable` (the named row
   exists, its coupling matches `COUPLINGS` and its shapes match `SHAPES` as a set; expectFail: KPZ with a renamed
   coupling). The model table of §6 is data in PythonExport.wl. Test with a float d.
9. **Contract** `docs/handoff/BMWflow_LO.md` = §6 verbatim, plus the generated KPZ module as an example.
10. **Benchmark** `benchmarks/KPZLO/` (files as in `KuramotoSivashinsky/`): `Reference.wl` = the .tex kernels A_a, A_c,
    A_d, I_D, I_ν, the η matrix and vector, β (alamk; g via g = alamk²); `Checks.wl` compares with `ExactlyZeroQ` at
    generic d; `run.wl`; notebook; a line in `benchmarks/README.md`; hook `tests/acceptance/KPZLO.wl`.
11. **Docs.** `docs/MODULES.md` sections (ExternalMomentum, PythonExport, LO branches); `docs/CONVENTIONS.md` (fn, LO
    names, measure excluded from Python kernels, f′ = d/d|momentum| in the export); `docs/FUTURE.md` (p = 0 exact limit /
    small-p series, causality export for functions, function relations, function × iω, f″/r″ in η rows, P2Divided ≥ 2,
    n_c > 1, `FixedDimension`/`jac` in LO, mixed ansätze, lessons for the new modules); fix the stale "shifted poles plug
    in without changes" (FUTURE.md:6-8, FRGDESIGN.md:193) and "only DerivativeExpansion" (FRGDESIGN.md:152,
    `lessons/src/01-Theory.md:142`, then `math -script lessons/makeLessons.wl`); FRGDESIGN §7 pointer to this record;
    `CLAUDE.md` map (new modules, design record, benchmark, `--require-python`).
12. **Harness.** `tests/runTests.wl` gets `--require-python` (`requirePythonQ`); `skip` takes a kind ("julia" | "python")
    so each flag turns only its own SKIPs into failures (runTests.wl:12, 21-23).

## 6. Contract (`docs/handoff/BMWflow_LO.md`, for the agent implementing BMWflow's `ModelGenerated`)
Generated file `eqs/<Model>_LO.py`: pure numpy + `math`, no state, no value-dependent branching (elementwise, broadcasts
over BMWflow's (p, w, q, θ) arrays).
- **Constants.** `FUNCTIONS` (order of every f tuple; KPZ `("fD", "fnu")`); `ETAS` (η_X ↔ FUNCTIONS[X]; KPZ
  `("etaD", "etanu")`); `COUPLINGS` (KPZ `("alamk",)`); `SHAPES` (regulator shapes in the order declared in
  `DefineRegulator`; KPZ `("rn", "rd")`); `P2_DIVIDED` (per function, the power of p² divided out of its loop, each ≤ 1;
  KPZ `(0, 1)`); `APPROXIMATION = "LO"`; `MODEL` (key into the model table; KPZ `"KPZ"`).
- **Conventions.** Dimensionless as in `kpz_notation.tex`; p, q in units of κ; f_X(0) = 1; primes on f are d/d|momentum|
  (d/dp in `linear`, d/dq in the η functions); for each shape `r = r(q)`, `r_ = dr/d(q²)`; s = ln(κ/Λ); every exported
  rate is ∂s x (BMWflow steps x → x − δs ∂s x).
- **Coupling state.** Exactly one exported coupling (n_c = 1). The map to BMWflow's scalar `g` is per model, chosen by
  `MODEL`; the adapter applies it (unknown `MODEL` → error).

  | MODEL | g ↔ coupling | β_g | shapes |
  |-------|--------------|-----|--------|
  | KPZ   | g = alamk²   | 2 alamk β_alamk | BMWflow's single `r` binds `rn` and `rd` |
- **Angles.** Every exported function depends on the angle only through c = cos θ.
- **Measure.** Nothing exported contains v_d, v_{d−1}, the measure's 2π, or Γ functions of d. The adapter applies them:
  I_X(p) = v_{d−1}/(2π) ∫_0^{q_max} dq q^{d−1} ∫_0^π dθ sin^{d−2}θ F_X (= ∫dᵈq/(2π)ᵈ F_X); in d = 1 BMWflow's rule
  (v_{d−1} → 1, θ ∈ {0, π}, weights 1). η: ⟨·⟩ = v_d ∫_0^{q_max} dq q^{d−1} (·).
- `linear(p, f, fprime, eta, couplings, d) -> tuple[n_f]`: linear part of ∂s f_X(p) (KPZ: η_X f_X + p f_X′).
- `loop(p, q, c, fq, fQ, rq, rq_, rQ, eta, couplings, d) -> tuple[n_f]`: loop integrand F_X; fq/fQ are tuples of f at q
  and at Q = sqrt(p² + q² + 2pqc) (spline at Q); rq, rq_, rQ are tuples over SHAPES. F_X is already divided by
  p^(2·P2_DIVIDED[X]) (non-finite at p = 0 for those rows: NaN or ±inf). ∂s f_X(p) = linear + I_X.
- **p = 0.** For every X with P2_DIVIDED[X] > 0 the adapter fills the p = 0 row with 'const' (value at p₁) or 'quad'
  (a + b p² through p₁, p₂), model-agnostic. Do not reuse `ModelKPZ`'s code (model_kpz.py:243-263): it is KPZ-only and
  divides by p² itself.
- **η.** `eta_matrix(q, f, fprime, r, r_, couplings, d) -> n_f×n_f nested tuple`, `eta_vector(...) -> tuple[n_f]` (r, r_
  tuples over SHAPES); solve Σ_Y (δ_XY + ⟨M[X][Y]⟩) η_Y = −⟨b[X]⟩ (row X = equation of f_X, column Y multiplies η_Y; M is
  not symmetric).
- `beta(couplings, eta, d) -> tuple[n_c]` (KPZ: alamk(d/2 − 1 + 3η_ν/2 − η_D/2)).
- **Acceptance**, d = 1 and 2, LO, `evo1`, the p/q grids and regulator of `usage_examples/KPZ_1D_twogrids`
  (`wett_alpha_beta`): `ModelGenerated` on the generated KPZ module vs `ModelKPZ` LO, comparing `self.Integral` after
  `Integral_upd` (p = 0 row after extrapolation, both modes), `eta_calc`, `g_rhs_calc`, and an `Evolution` run of N = 2000
  steps (ds = 1e-3). Norm: per function, max_grid |Δ| / max_grid |x| ≤ 1e-12 for single evaluations, ≤ 1e-11 after the
  run; η and g relative ≤ 1e-11. BMWflow's NLO paths untouched. Environment: BMWflow needs scipy and numba (not installed
  for /usr/bin/python3 on this machine; numpy is): set up a venv.

## 7. Verification
- `math -script tests/runTests.wl` exit 0: DE + new unit tests; the step-0 baseline unchanged.
- `math -script benchmarks/KPZLO/run.wl` exit 0; notebook checked with `math -script examples/runNotebook.wl
  benchmarks/KPZLO/KPZLO.nb`.
- Python smoke test (numpy only; SKIP without numpy, failure with `--require-python`): import the exported KPZ module,
  evaluate every function at a few points including a float d, compare with Mathematica `N[]` values to 1e-13.
- The numerical BMWflow acceptance (§6) belongs to the other agent.

## 8. Reviews
Seven rounds, details and probes in `design/reviews/LOPipeline-r1-r7.md`. Blockers: r1 (g vs alamk), r4 and r5 (DE
baseline not reproducible); majors: r1 (p = 0 handling; f″/r″), r2 (contract out of sync; measure stated twice), r3
(no check of low orders before dividing by p²). All resolved above. r6 and r7 found no blockers or majors.

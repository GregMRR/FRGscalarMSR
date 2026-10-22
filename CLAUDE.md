<!-- verified: 6208fe0 2026-10-21 -->
# FRGscalarMSR
Mathematica library: an MSR ansatz for Γ_k of N scalars → checked one-loop Wetterich flow equations (derivative expansion) →
dimensionless β functions → FlowRG (Julia) equations file.

## Commands
- Kernel: `math -script` (Wolfram 15; `wolframscript` is not installed).
- All tests: `math -script tests/runTests.wl` (exit 0 = all pass). One file: `--only DimensionsTest` (substring of the file's
  base name). `--verbose` shows the cards and FAIL blocks; `--require-julia` turns the SKIPs of julia/FlowRG tests into failures.
- Examples: edit `examples/makeNotebooks.wl`, then `math -script examples/makeNotebooks.wl`; check one headless with
  `math -script examples/runNotebook.wl examples/KPZ.nb`.

## Layout
- `FRGscalarMSR.wl` - loader; `$FRGModules` lists `Kernel/*.wl` in load order (each module uses only those above it).
- `Kernel/` - one module per pipeline stage: Theory → Ansatz → Vertices, Regulator → Propagator → Symmetries → Loops → Projection →
  FrequencyIntegration → Closure → Assemble (`Equations` stage) → Dimensions → Adimensional → JuliaExport; Core (asserts, stages,
  `ExactlyZeroQ`), Oracle (numerical cross-check), Stages (save/load).
- `docs/MODULES.md` - per module: objective, assumptions, manipulations, checks, Show commands. `docs/CONVENTIONS.md`,
  `docs/FUTURE.md`.
- `design/` - decision records with their adversarial reviews: `FRGDESIGN.md` (v1), `LEVICIVITA.md` (`jac`), `ADIMENSIONAL.md`
  (dimensions, export), `MOMENTUMPROJECTION.md` (p^{2n} projection, any n), `CAUSALITY.md` (causal-pole conditions exported as
  `causal_poles`); `design/reviews/`.
- `docs/handoff/` - specs for work in other repos (FlowRG: `MUST_FIX`, `D2<shape>`).
- `tests/reference/` - shared pipeline setups (`ShallowWaterSetup.wl`, `KPZSetup.wl`), `JuliaRoundTrip.wl`, `OldProjection.wl`
  (the pre-p^{2n} Laplacian operator, test-only regression reference); `legacy/` holds the old
  notebooks' results used as references (`legacy/README.md` maps the old notation).

## Where to change X
- New module: `Kernel/X.wl` + its name in `$FRGModules` + a section in `docs/MODULES.md` + `tests/unit/XTest.wl`.
- New check: `FRGAssert`/`assertCheck` named `Module.check` (shown on the card) + an `expectFail[…, "Module.check"]` test.
- FlowRG file format: `Kernel/JuliaExport.wl`; the contract is `~/Physics/Sources/FlowRG.jl/docs/equations-interface.md`
  (other repo: changes there go through a `docs/handoff/` spec).
- A design change: record it in `design/` (decision + reason) and update `docs/MODULES.md` in the same change.

## Conventions
- Every check is a named assertion that aborts with a red FAIL block - the user reads cards, not code (grey-box modules).
- Stages are immutable, fingerprinted Associations; constructors call `validateStage` on their inputs - no global state.
- Exact zero tests use `ExactlyZeroQ` (random rationals), never `Simplify[…] == 0`; no `Simplify` on large expressions, no radicals
  (`Together`/`Collect`/`Factor`) - expressions reach LeafCount 10⁵–10⁶.
- `Dim` stays symbolic to export; d is fixed only by `FixedDimension[2]` (required by `jac`).
- ∂s ref = −η ref; dimensionless names `a<>c`; loop momentum q = k aq (`docs/CONVENTIONS.md`).
- User notebooks (`.nb` in `~/sync_cnrs/Projects/…`) and generated `eqs/*.jl` are the user's: hand them instructions, do not edit.
- Never kill `WolframKernel` processes this session did not start (the user's front-end kernels run alongside).

## Gotchas
- Every module's private code shares one context, ``FRGscalarMSR`Private` ``: a helper name defined in one module is visible to,
  and can be overwritten by, every other (e.g. Regulator's `sq`). Pick distinctive names.
- `Context[expr]` holds its argument: use `Context @@ {expr}` on a computed symbol.
- Test files run in alphabetical order in one kernel and share globals; heavy setups are built once behind `If[! ValueQ[…]]`
  (`ShallowWaterSetup.wl` takes ~85 s). The harness defines `test`, `expectFail`, `skip`, `quiet`, `testRoot`.
- Generated DE4plus Julia functions are huge: full compilation takes ~15 min per file. Tests run julia with
  `--compile=min -O0` (`$juliaFlags` in `tests/reference/JuliaRoundTrip.wl`).
- Regulator prefactors may contain reference couplings only (`AdimensionalEquations` aborts otherwise); a model with independent
  diffusion constants in two sectors' regulators hits this (`design/ADIMENSIONAL.md` §4).
- FlowRG path: `~/Physics/Sources/FlowRG.jl`, or `FLOWRG_PATH`.

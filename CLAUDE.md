<!-- verified: 33e790e + branch NLO 2026-11-12; BMWflow handoff v4 2026-11-12 -->
# FRGscalarMSR
Mathematica library: an MSR ansatz for Γ_k of N scalars → checked one-loop Wetterich flow equations (derivative expansion) →
dimensionless β functions → FlowRG (Julia) equations file. Second pipeline (LO): couplings that are functions of p² (`fn[f, X]`) →
flows at generic p → BMWflow Python module (design/LOPIPELINE.md; functions times iω and in vertices, KPZ f_λ: design/LOFLAMBDA.md;
NLO, functions of (ϖ, p) with `BuildLOEquations[…, "ExternalFrequency" -> True]`: design/NLOPIPELINE.md).

## Commands
- Kernel: `math -script` (Wolfram 15; `wolframscript` is not installed).
- All tests: `math -script tests/runTests.wl` (exit 0 = all pass). One file: `--only DimensionsTest` (substring of the file's
  base name). `--verbose` shows the cards and FAIL blocks; `--require-julia` (`--require-python`) turns the SKIPs of julia/FlowRG
  (python3 + numpy) tests into failures.
- Examples: edit `examples/makeNotebooks.wl`, then `math -script examples/makeNotebooks.wl`; check one headless with
  `math -script examples/runNotebook.wl examples/KPZ.nb`.
- Benchmarks: `math -script benchmarks/<Name>/run.wl` (exit 0 = all pass); notebook: edit `benchmarks/<Name>/makeNotebook.wl`,
  run it, check with `examples/runNotebook.wl benchmarks/<Name>/<Name>.nb`.
- Lessons: edit `lessons/src/NN-Name.md`, then `math -script lessons/makeLessons.wl`; check one headless with
  `math -script examples/runNotebook.wl lessons/NN-Name.nb`.
- Docs: `python3 tools/docgen.py reference` (regenerate `docs/reference/`; the pre-commit hook in `.githooks/` runs it),
  `python3 tools/docgen.py check` (stale reference or broken relative link → exit 1; CI runs it). Site: `/project-docs publish`.

## Layout
- `FRGscalarMSR.wl` - loader; `$FRGModules` lists `Kernel/*.wl` in load order (each module uses only those above it).
- `Kernel/` - one module per pipeline stage: Theory → Ansatz → Vertices, Regulator → Propagator → Symmetries → Loops → Projection →
  FrequencyIntegration → Closure → Assemble (`Equations` stage) → Dimensions → Adimensional → JuliaExport; Core (asserts, stages,
  `ExactlyZeroQ`), Oracle (numerical cross-check), Stages (save/load). LO pipeline: Loops → ExternalMomentum (`BuildLOEquations`,
  replaces Projection…Assemble) → Dimensions → Adimensional (LO branches) → PythonExport; `"Pipeline"` is inferred by the Ansatz.
- `README.md` - landing page (install, quickstart, links). `docs/README.md` - docs index + `verified` stamp; `docs/STATUS.md` -
  implemented, checks, limitations, future directions; `docs/guide/` - how-tos (DE parts 1–2, LO, stages, development).
- `docs/reference/` - GENERATED from the `::usage` strings (`docs/docgen.toml` lists the sources in load order): edit the usage
  string, never the page. `tools/docgen.py` - vendored copy of `~/.claude/skills/project-docs/docgen.py`.
- `docs/MODULES.md` - per module: objective, assumptions, manipulations, checks, Show commands. `docs/CONVENTIONS.md`.
- `design/` - decision records with their adversarial reviews: `FRGDESIGN.md` (v1), `LEVICIVITA.md` (`jac`), `ADIMENSIONAL.md`
  (dimensions, export), `MOMENTUMPROJECTION.md` (p^{2n} projection, any n), `CAUSALITY.md` (causal-pole conditions exported as
  `causal_poles`), `LOPIPELINE.md` (functions of p², BMWflow export), `LOFLAMBDA.md` (reference-1 functions at iω and in vertices, ∂/∂(iν)
  flows, contract v3), `NLOPIPELINE.md` (external frequency in the loop, shifted pole family, contract v4); `design/reviews/`.
- `docs/handoff/` - specs for work in other repos (FlowRG: `MUST_FIX`, `D2<shape>`; BMWflow: `BMWflow_LO.md`, the LO contract).
- `tests/reference/` - shared pipeline setups (`ShallowWaterSetup.wl`, `KPZSetup.wl`, `KPZLOSetup.wl`, `KPZNLOSetup.wl`), `JuliaRoundTrip.wl`,
  `PythonRoundTrip.wl`, `DEBaseline.wl` (DE content hashes, checked by `tests/acceptance/DEBaseline.wl`), `OldProjection.wl`
  (the pre-p^{2n} Laplacian operator, test-only regression reference); `legacy/` holds the old
  notebooks' results used as references (`legacy/README.md` maps the old notation).
- `benchmarks/` - literature benchmarks (`KuramotoSivashinsky`, `KPZLO`, `KPZLOLambda`, `KPZNLO`), one folder each (`README.md` is the index): `README.md`, `Reference.wl` (the paper's
  equations), `Pipeline.wl`, `Checks.wl` (`ksChecks[]`-style records, no harness dependency), `run.wl`, `makeNotebook.wl` → `.nb`.
- `lessons/` - one teaching notebook per module (`README.md`: index, source format): sources `src/NN-Name.md` → `makeLessons.wl`
  → `NN-Name.nb` + `Glossary.nb`; `LessonSetup.wl` (`lessonUpTo[stage]` rebuilds KPZ, `showFailure[expr]`).

## Where to change X
- New module: `Kernel/X.wl` + its name in `$FRGModules` + a section in `docs/MODULES.md` + `tests/unit/XTest.wl`.
- New benchmark: `benchmarks/<Name>/` (files as in `KuramotoSivashinsky/`), a line in `benchmarks/README.md`, a suite hook
  `tests/acceptance/<Name>.wl` that only wraps the checks in `test[…]`.
- New check: `FRGAssert`/`assertCheck` named `Module.check` (shown on the card) + an `expectFail[…, "Module.check"]` test.
- BMWflow module format (LO): `Kernel/PythonExport.wl`; the contract (versioned, `CONTRACT_VERSION`) is `docs/handoff/BMWflow_LO.md`,
  the reader BMWflow's `flow/models/model_generated.py` (v4: LO and NLO kinds). A format change updates both repos, bumps the
  version, and re-vendors BMWflow's `tests/fixtures/` (`KPZ_LO.py`, `KPZLambda_LO.py`, and the NLO files written by
  `benchmarks/KPZNLO/export.wl`; `tests/acceptance/BMWflowLO.wl` fails until then).
- FlowRG file format: `Kernel/JuliaExport.wl`; the contract is `~/sources/FlowRG/docs/equations-interface.md`
  (other repo: changes there go through a `docs/handoff/` spec). The two sides point to each other: a change to what
  `JuliaExport.wl` writes must match that page (FlowRG's `CLAUDE.md`, golden rule 6, points back here), and FlowRG answers a
  handoff with a Status section listing what this side must update.
- A design change: record it in `design/` (decision + reason) and update `docs/MODULES.md` and `docs/STATUS.md` in the same change.
- A public symbol: its `::usage` string is its reference entry; new limits or features go to `docs/STATUS.md`.
- A module's behaviour or a quoted symbol: update its lesson `lessons/src/NN-Name.md` and rebuild (the build fails on a missing
  quoted symbol).

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
- The checkouts are in `~/Physics/Sources/` (FRGscalarMSR, FlowRG.jl, BMWflow); `~/sources/` does not exist, but some notebook
  and run headers still name it. FlowRG: set `FLOWRG_PATH`. BMWflow: found at `~/Physics/Sources/BMWflow`, or `BMWFLOW_PATH`; its
  tests need its `.venv` (scipy, numba, pytest; not installed for /usr/bin/python3).
- Stage fingerprints are not reproducible between kernel sessions: compare stage contents (`===` on rows), not fingerprints,
  across runs (design/NLOPIPELINE.md §6).

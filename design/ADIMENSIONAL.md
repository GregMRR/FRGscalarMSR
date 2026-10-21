# ADIMENSIONAL — dimensions, dimensionless flows, FlowRG export

Design record of the three modules after `Equations`: `Dimensions.wl` (`DeriveDimensions`), `Adimensional.wl`
(`AdimensionalEquations`), `JuliaExport.wl` (`ExportJulia`). It covers the decisions and why they were taken, how the design was
reviewed, what was built and checked, and what is still open.
- Module behaviour, checks and Show commands: `docs/MODULES.md` (sections Dimensions, Adimensional, JuliaExport).
- Conventions (∂s ref = −η ref, names, measure): `docs/CONVENTIONS.md`.
- Decision table D1–D18, short form: `design/FRGDESIGN.md` §11.
- Adversarial review of the plan: `design/reviews/Adimensional-r1.md` (2 blockers, 4 majors; all resolved below).
- FlowRG side (separate repo): `docs/handoff/FlowRG-MUST_FIX.md`.

## 1. Why

The step "dimensionful equations → dimensionless β functions → FlowRG `.jl` file" used to be done by hand in
`~/sync_cnrs/Projects/shallow_water/ManualCodeOld/DE4plus_largscale_forcing_LLM/AdimensionalFlowEQGenericRegLargeScale.nb`:
hand-typed `Dimensionality` vectors, `LinEQ` linear parts, `HomClean` homogeneity, a hand-built η system and `ExportDE4`.
`~/sync_cnrs/Projects/shallow_water/AdimensionalLPAsw.nb` (the user's notebook) reloads the LPA `Equations` stage and is meant to
use these modules. The goal: derive every dimension from the ansatz, keep every step inspectable, fail loudly.

## 2. Decisions (grilling with the user, 2026-10-21)

1. **Dimensions are derived, not typed.** Γ_k = ∫dt dᵈx 𝓛 is dimensionless, [x] = 1/k. Each dimension is an exponent vector over
   the basis {k, ref₁…ref_{N+1}}: N references for the field dimensions, one for time. The user chooses the references; exactly
   N+1 free units besides k is asserted (it fails, for instance, when a φ̄ dt φ term carries a coupling).
2. **Amplitudes are declared** (`"Dimensionless" -> {Ak, aDphik}`): in `nk Ak` one monomial cannot split the dimension between two
   couplings. An underdetermined split aborts and names the couplings to declare.
3. **The linear part follows from the dimension**: ∂s ĉ = ĉ(−x + Σ y_j η_j) + (∂s c)/[c]. The user only gives the state
   (`"State"`, the COUPLINGS order). A Leibniz linear solve decides which flowing combination gives which ∂s c; NonRenormalised
   couplings keep the linear part only.
4. **Two kinds of reference** (user's physics argument). A reference that appears alone in an equation (`nk`) is *running*: its η
   comes from that equation. A reference the user wrote times an amplitude (`Dphik aDphik`) is *fixed*: its η will be fixed by hand
   (e.g. a large-scale forcing used as regulator, η_φ chosen to keep the energy injection finite), the amplitude carries the
   generated term, and its equation is the amplitude's flow. The old hand-built convention (η from the amplitude equation,
   amplitude counted as 1) stays available as `"Mode" -> "Running"` (needed for DE4plus, review E3).
5. **Fixed η at run time in FlowRG**, not in Mathematica: placeholder row η = 0 in the η system, `const MUST_FIX = (:phi,)` in the
   file, and FlowRG refuses to load the file unless `fixed_etas` gives the value (user: more flexibility to explore). Considered
   and rejected: a model parameter (interface rule 5 change), a number baked into the file (regenerate per value), η as a state
   variable with β = 0 (zero Jacobian row).
6. **Regulator prefactors: references only.** `nk (1 + Ak) q^2 arnh` would put β_A itself into ∂sR. The exact treatment (one
   unknown vector X = (η…, ∂s Ak…), each β written once, FlowRG `SOLVED`) was judged confusing for now; the user chose to forbid it
   and record the extension in `docs/FUTURE.md`. Symbolic pre-solving of the η system is impossible while the shapes stay abstract:
   its coefficients are aq-integrals that exist only numerically.
7. **Measure** S_d aq^{d−1}/(2π)^d. The old exports lacked the (2π)^{−d} (review E7); Canet's v_d = S_d/(2π)^d settles it.
8. **r″ is exported** as `D2<shape>(aq^2)` (DE4plus g2, D2φ need it, review E5); FlowRG binding is a handoff task. Open question
   left to the user: Litim's r″ contains δ(y − 1).
9. **The old ∂sM_h ∝ η_φ was a bug**; the derivation gives η_h.
10. **Naming**: `a<>c`, dimensionless couplings keep their name, η named by the reference (`etanu`, `:nu`), collisions abort,
    `"Names"` overrides.
11. **Compact export**: the scale-stripped atoms of Assemble become dimensionless atoms, written as Julia local variables.
12. **Three modules**, `Dim` symbolic until export, one file per d; dimensions computable right after `DefineAnsatz`.

## 3. Review findings that changed the design (`design/reviews/Adimensional-r1.md`)

- **E1 (blocker)**: the plan's monomial equation had the measure with the wrong sign; g came out right, H and λ wrong. Fixed
  (+[dt dᵈx]); the unit test checks all couplings against the old notebook's hand values, not only gk.
- **E2 (blocker)**: D4 rejected the user's own LPA regulator `nk (1 + Ak) q^2 arnh`. User decision: change the LPA regulator to
  `nk q^2 arnh` (D18).
- **E3**: in DE4plus every reference appears only times an amplitude, so the strict rule made all η fixed → `"Mode" -> "Running"`
  and `"Equation" -> comb` (nk appears in both `nk ank` and `nk Ak`).
- **E4**: the KPZ reference is Canet Eq. (20), not `LPA_KPZ_thermal_forcing/eqs` (which is shallow-water LPA).
- **E5**: r″ abort would block DE4plus → export `D2`.
- **E6**: the density is not invariant under rescaling, the action is → check restated on the action, with test fields in the raw
  density.
- **E7–E12**: (2π)^{−d} decided in advance, rescaling check moved before unit stripping, summed coefficients handled per term, names
  of derived-dimensionless couplings, S_d defined, `AnsatzOf` for loaded stages, visible SKIP + `--require-julia`, stale paths,
  `betal` never omitted.

## 4. What was built and how it is checked

- `Dimensions.wl`: exact linear solve (normal equations, full column rank asserted first), freedom count by null space, reference
  independence, consistency with the conflicting monomials listed, action-invariance check on polynomial test fields (d = 2, 3),
  names. `AnsatzOf`/`EmbeddedStage` reach the ansatz inside any later stage.
- `Adimensional.wl`: D4 and NonRenormalised-reference checks, R_ab and LHS dimension checks, reference classification, Leibniz
  solve, unit substitution, `homVec` homogeneity (exponents linear in `Dim`, units forbidden in exponents and function arguments),
  random unit rescaling of every loop piece before stripping, compact (AdimAtom) = explicit for betas and η rows.
- `JuliaExport.wl`: FlowRG interface, positional arguments, measure multiplied in, `D`/`D2` shape names, atoms as locals, zero
  entries written `0.0`, declared d = 2 enforced.

Acceptance (all passing):
- `tests/acceptance/KPZAdimensional.wl`: the dimensionless KPZ flow equals Canet's Eq. (20) exactly for symbolic d (fixes the
  (2π)^{−d}); Julia round trip and FlowRG load in d = 1, 2, 3.
- `tests/acceptance/ShallowWaterDE4Adimensional.wl`: DE4plus against the old `eqs/Equations2D.jl`
  (`tests/reference/legacy/DE4plusEquations2D.jl`): `betal` equal; `betanl`, η-system M·η and b equal to old × (2π)^{−2} at
  η_h = η_φ.
- `tests/acceptance/ShallowWaterLPAAdimensional.wl`: the user's LPA ansatz (regulator `nk q^2 arnh`) through SaveStage/LoadStage,
  hand dimensions, η_φ fixed (MUST_FIX), Julia round trip, FlowRG load with `fixed_etas = (phi = 0.1,)`.
- `tests/acceptance/Jacobian2D.wl` (appended): declared d = 2, export in d = 3 aborts, round trip.
- `tests/unit/DimensionsTest.wl`, `tests/unit/AdimensionalTest.wl`: every abort, fixed/running classification.

Found while building:
- The `jac` test model's regulator `mk q^2 jrp` violates D4: [mk] = [nk], so mk cannot be a fourth reference. Any model where two
  fields carry independent diffusion constants as regulator prefactors needs the FUTURE.md extension (decision 6).
- Julia compiles the huge generated functions slowly (DE4plus: ~15 min per file, the old file too). The tests run julia with
  `--compile=min -O0` (seconds). FlowRG's own load of a large file still pays the full compile.

## 5. Open

- FlowRG: `MUST_FIX` guard and `D2<shape>` binding (`docs/handoff/FlowRG-MUST_FIX.md`); Litim's r″ δ term is a user decision.
- The user's notebooks: rerun `ShallowWaterLPA.nb` with `{Hbar, H} -> nk q^2 arnh`, resave `swLPAplusequations.wdx`, then in
  `AdimensionalLPAsw.nb`: `DeriveDimensions[AnsatzOf[equations], …]`, `AdimensionalEquations`, `ExportJulia` (README workflow).
- `docs/FUTURE.md`: non-reference regulator couplings (`SOLVED`), η fixed by a linear relation or a scale-invariant quantity.

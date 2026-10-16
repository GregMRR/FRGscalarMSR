# HANDOFF — FRGscalarMSR, state after milestone 4 (2026-10-16)

This is for the agent who takes over. Read it first, then `FRGDESIGN.md` (the approved design, binding) and `docs/MODULES.md`
(what each module does). Milestones 1–4 of design §8 are done and committed. Milestones 5 and 6 remain.

## 1. Ground rules

- **Workflow.** The user wants a pause after each milestone: implement it, run the full suite, commit, report in a few lines,
  then wait for "go".
  - Commit messages end with the `Co-Authored-By` / `Claude-Session` lines the harness supplies.
  - Push to `origin` (`git@github.com:GregMRR/FRGscalarMSR.git`, SSH; there is no `gh` CLI) only when asked.
- **The user.** A theoretical physicist (FRG, MSR field theory). Explanations must not assume jargon: say what Vieta,
  Routh–Hurwitz and the like mean whenever they appear, in docs and in printouts.
- **Binding computer-algebra rules.** `~/sync_cnrs/Projects/shallow_water/AGENT_PRIMER.md` §7:
  - keep atoms abstract until the end;
  - no radicals;
  - use `Together`/`Collect`/`Factor`, not `Simplify`;
  - verify exactly with random rationals (`ExactlyZeroQ`), never `Simplify[…] == 0`.
- **Kernel.** Only `math -script` (Wolfram 15.0); `wolframscript` is not installed.
  - Never kill `WolframKernel` processes you did not start: the user runs FrontEnd kernels.
  - Never edit `.nb` files or `eqs/*.jl` in the user's projects.
  - Never modify `~/sync_cnrs/Projects/shallow_water/` (the reference project).
- **Global instructions.** The user's `~/.claude/CLAUDE.md` asks for a token-lean style: narrow checks, filtered output, short
  closing reports. Do not create a project `CLAUDE.md` unless asked.

## 2. Running things

```
math -script tests/runTests.wl            # unit + acceptance, ~26 s, exit 0 iff all pass (186/186 now)
math -script tests/runTests.wl --verbose  # also prints every card, table and FAIL block
```

- `runTests.wl` loads the package, sets `$FRGExitOnFail = False`, and runs `tests/unit/*Test.wl`, then
  `tests/acceptance/*.wl`.
- The harness provides:
  - `test[name, actual, expected]`, which compares with `===`;
  - `expectFail[name, expr, "Check.name"]`, which requires an abort through that named check (read from `$FRGLastFailure`);
  - `quiet[expr]`, which hides Print.
- `tests/reference/ShallowWaterSetup.wl` builds the whole shallow-water pipeline once, guarded by `ValueQ`. All acceptance
  files `Get` it. It defines `swaTheory … swaIntegrated`, the atom map `toScriptAtoms`, `refForms` and `regOff`.
- `tests/reference/DE4plusReference.wl` is a mechanical copy of parts of `DE4plus_fast.wl`, with heads renamed: Γ2 →
  refGamma2, Γ3 → refGamma3, Γ → refV, η/γ/k/l/h → refEta/refGam/refK/refL/refH, Δ → refDelta, G → refG. It never runs the
  script.
- `tests/reference/flowsFAST/*.wdx` are the reference equations. Each file = raw RHS / LHS factor, from the `jobs` table of
  `DE4plus_fast.wl` §5.
- Output is often huge. Pipe it through `cut -c1-300 | tail`. When a run "returns unevaluated", look first for a syntax error
  in a Kernel file: load it alone with
  `math -noprompt -run 'Get["Kernel/Core.wl"]; Check[Get["Kernel/X.wl"],"ERR"]; Exit[]'`.

## 3. Architecture (as built)

- `FRGscalarMSR.wl` is the loader. `$FRGModules` lists the modules in load order:
  Core, Theory, Ansatz, Vertices, Regulator, Propagator, Symmetries, Loops, Projection, FrequencyIntegration.
  **Add new modules there.**
- **Contexts.**
  - Public symbols live in `FRGscalarMSR``.
  - Shared internals live in `FRGscalarMSR`Internal``: `assertCheck`, `infoCheck`, `stampStage`, `validateStage`,
    `assertFreshSymbols`, `inputForm`, `displayString`, `reservedNames`.
  - Each file has its own `Begin["`Private`"]`. Private symbols are `FRGscalarMSR`Private`…`, shared across files; be careful
    with names.
- **Stages** are Associations with the keys `Stage`, `Pipeline`, `Upstream` (fingerprints), data, `Checks`, `Warnings` and
  `Timing`.
  - `stampStage` adds `Fingerprint` = SHA256 of `Normal[KeyDrop[stage, volatile], Association]`. `Normal` re-evaluates, so a
    symbol that later gets a value invalidates the stage.
  - Every constructor starts with `validateStage[input, "Name"]`.
  - Cards come from `PrintCard`.
- **Failures.** `FRGAssert[name, cond, details, hint]` (HoldRest) prints a FAIL block, sets `$FRGLastFailure`, then calls
  `Abort[]`, or `Exit[2]` under `math -script`. `assertCheck` returns a PASS record for the card.
- **Data flow** (the user-facing API is design §4):

| stage | constructor | key data |
|---|---|---|
| Theory | `DefineTheory["Scalars"->{H,Psi},"Dimension"->Dim]` | `Fields` (φ1,φ̄1,φ2,φ̄2), `ResponseFields`; bars auto-created as `<name>bar` |
| Ansatz | `DefineAnsatz[theory,"Density"->…,"Couplings"->{…}]` | `Monomials`: coefficient × `FieldFactor[φ, n_t, {SpaceIndex[1]…}]` |
| Vertices | `BuildVertices[ansatz]` | `VertexTable`: sorted field multiset → Γ⁽ⁿ⁾ in private leg variables `wv[j]`, `qv[j]`; accessors `VertexValue`, `Gamma2Matrix`, `Gamma3Matrix` |
| Regulator | `DefineRegulator[ansatz,{a,b}->R,…,"Shapes"->{…}]` | `RMatrix`, `DsRMatrix` (Functions of s = q²), `DsRAbstract` (`dsR[a,b]`), `Scale` (the user's `k`); `dsCoupling[c]` = k dc/dk |
| Propagator | `BuildPropagator[vertices, regulator]` | atoms `GammaAtom[a,b,m][s]` (coefficient of (iω)^m, a ≤ b), `AtomValues`, `G`, `Delta`, `DeltaCoefficients` {c0,c1,c2} in Δ = c0 + c1(iω) + c2(iω)², `ValidityConditions`, `DeltaDegree` |
| Symmetries | `DeclareSymmetries[ansatz, Protected[…], Neglected[…], NonRenormalised[…], "TruncationOrders"->…]` | `Classification`: one record per (pair, m, n), with class Declared, DeclaredNumeric, Protected, Neglected, Causality or Structure, the `GammaK` coefficient and a reason |
| Loops | `BuildLoopIntegrands[vertices, prop, regulator]` | `LoopIntegrand[loops,{a,b},ν,p]` in loop variables `LoopVariables[loops]`; tadpole included when Γ4 ≠ 0 |
| Projection | `ProjectOnMonomials[loops, symmetries]` | `Integrands`: classification records plus `Integrand` (in ω, `LoopMomentumSquared`, atoms, `Δ[±ω,s]`, `DeltaCoefficient`, `dsR`) |
| Integrated | `IntegrateFrequencies[projected]` | `Flows`: `Declared` → `"Value"` (symbolic); others → `"Samples"` (exact values at random points; 3 seeds, 1 for Neglected) |

## 4. Decisions and deviations from FRGDESIGN.md (all reported to the user)

- `FRGAssert` instead of `Assert`, because System`Assert is protected.
- `Kernel/Symmetries.wl` is **new** (not in design §5). `DeclareSymmetries` must run before `ProjectOnMonomials`, and it fails
  early on untagged monomials and on contradictions.
- As the design asks, projection is a separate stage **before** frequency integration, and `CheckClosure` replaces the
  previewed `ProjectFlows`. Both points were flagged to the user.
- **Reserved names** are only `k q p y t x ω ν`. Internal names are private, so the user may name a coupling `nu`, `A`, …
- **Couplings may not sit under an operator.** `dt[nk H]`, `lap[nk H]` and the like all abort.
  - System `Laplacian[H,{x,y}]` evaluates to 0 before the parser sees it, so it cannot be caught; only `Dt` is detected. This is
    documented.
- **Regulator.** The user writes their own `q` and `k`; they are recognised by `SymbolName`. The user's `k` is the scale symbol
  in every output.
- **Pattern matching in tags is unordered.** `{Hbar, _}` also matches `{Hbar, Hbar}`. In KPZ this correctly aborts as a
  contradiction (D is declared); the KPZ shift-symmetry tag must be written `{h, hbar}`.
- Causality pairs (Γ_φφ) are not skipped. Their flows are computed and sampled, and all are exactly 0.
- **Sampled integration.** Protected, DeclaredNumeric, Causality and Neglected flows are integrated **at random rational atom
  points**, never symbolically.
  - The Δ coefficients must be recomputed from the atoms at each point (`samplePoint` in FrequencyIntegration). Randomising
    `DeltaCoefficient` independently breaks the identities that rely on Δ = det B.
  - The same trap applies when testing a symbolic `Value` with `ExactlyZeroQ`: apply `explicitFlow` first (see
    `tests/acceptance/ShallowWaterFlows.wl`).
- **The ω² Neglected values are sampled at arbitrary points.** These points do not respect the validity conditions
  (c1/c2 > 0, c0/c2 > 0). That is fine for zero tests, but not for reporting a "size". **M5 must resample until the validity
  conditions hold** before reporting the size of a Neglected term.

## 5. Next: milestone 5 (Closure, Assemble, Oracle, Stages)

Design §6 items 9–12, §2 Q6/Q10/Q12/Q13/Q16/Q21 and §10 hold the specification. Concretely:

1. **`Kernel/Closure.wl` `CheckClosure[integrated, symmetries]`.**
   - Protected, DeclaredNumeric and Causality flows: all samples must be 0, otherwise abort with the tag's reason.
   - Neglected: print a WARN with the value at a *valid* random point.
   - Declared: equation ∂s(GammaK) = Value / normalisation.
     - The normalisation is the numeric content of GammaK relative to the flowing combination. For example,
       GammaK = −2 Dhk aDhk gives ∂s(Dhk aDhk) = Value/(−2).
     - Q10: the output is ∂s of the coefficient as written in the DSL. Agree with the user how to split the numeric factor
       before coding it: the current plan is "numeric factor out, the rest is the flowing combination".
   - Print the `NonRenormalised` assumptions and the Ward-identity banner.
   - `ShowClosureReport`.
2. **`Kernel/Assemble.wl` `AssembleEquations[flows]`.**
   - Port `compactify`/`collectNum`/`compactFrac` (script §5), with **auto-generated** scale-stripped atoms built from
     `AtomValues`.
   - For each atom value, factor out its coupling/q² scale, the way `toAtoms` writes feta = nk q² xE. The sync assertion becomes
     an exact check of atoms → values.
   - Only at assembly, substitute the derived ∂sR (`DsRRules`) and express the `DeltaCoefficient`s via atoms.
   - The output contract per equation (design §6 item 10, oriented to FlowRG `~/sources/FlowRG/docs/equations-interface.md`):
     - the loop q-integrand;
     - `dsCoupling[…]` entering linearly;
     - shapes r, r′, r″ abstract in y = q²/k²;
     - the validity conditions;
     - `Dim` symbolic.
   - Assert that no intermediate symbols are left.
   - Watch the LeafCounts: g2 is 45k leaves and D2phi 165k in the script.
3. **`Kernel/Oracle.wl` `RunOracle[equations]`.** A generalised script §6c.
   - Numeric inverse of Γ2 + R, `NIntegrate` over ω, finite differences in ν and in an explicit p vector.
   - Run in **d = 2 and d = 3**.
   - It must be independent of Propagator, Projection and FrequencyIntegration; it shares only Γ3 and the loop formula.
4. **`Kernel/Stages.wl`.** `SaveStage`/`LoadStage`, fingerprint-checked, plus the remaining `Show…` dispatchers
   (`ShowEquation`, `ShowClosureReport`).
5. **Tests and acceptance.**
   - The negative tests of design §10 not yet written: an untagged monomial, a protected term with nonzero flow, a stale
     `LoadStage`.
   - The script §6b numbers through `CheckClosure`.
   - The oracle agrees to about 1e-10 on the shallow-water equations.
   - Compacted equations equal `flowsFAST` exactly, and their LeafCounts are comparable.

## 6. Then: milestone 6

- **KPZ.** Compare the dimensionful ν and D flows with a hand derivation (textbook one-loop) and with the dimensionful content
  of `~/sync_cnrs/Projects/shallow_water/LPA_KPZ_thermal_forcing`.
  - This is the only cross-read of the loop normalisation (E3/E4). State honestly what it does and does not test.
  - The Galilean iω check already passes (`tests/unit/SymmetriesTest.wl`).
- **Example notebooks** `examples/ShallowWaterDE4plus.nb` and `examples/KPZ.nb`. Generate them programmatically with
  `Notebook[{Cell[…,"Input"]…}]` + `Export`, then check they run headless.
- **Docs.** Finish README: the known limitation E3/E4, the workflow, `FUTURE.md` cross-links.

## 7. Where things are

- **Design:** `FRGDESIGN.md`; the original is in `~/sources/FRGflowMSR/`.
- **Module documentation:**
  - `docs/MODULES.md` covers every module so far, including the hand derivation of the loop formula.
  - `docs/CONVENTIONS.md`: Fourier convention, legs, DSL, reserved names.
  - `docs/FUTURE.md`: the future seams.
- **Physics and reference:** `~/sync_cnrs/Projects/shallow_water/AGENT_PRIMER.md` and `DE4plus_fast.wl`.
- **Julia target:** `~/sources/FlowRG/docs/equations-interface.md`.
- **Commits:** M1 `6da1937`, M2 `66579cc`, M3 `356862a`, M4 `2a4fe12`.

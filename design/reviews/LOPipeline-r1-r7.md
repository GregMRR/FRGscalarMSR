# Reviews of design/LOPIPELINE.md (plan stage), rounds 1–7

The reviews ran in plan mode and were returned inline; this file records them. Rounds 1–6: the-arch-enemy; round 7:
the-enemy. Attacked commits: FRGscalarMSR `21719c9` (BMW-LO), BMWflow `55108ad`. Section numbers (§, L#, step #) refer to
`design/LOPIPELINE.md`; findings were raised against earlier wordings of the same plan.

## Round 1 — 1 blocker, 2 major, 9 minor, 3 nit
Held: contract measure (S_{d−1}/(2π)^d = v_{d−1}/(2π); d = 1 rule vs model_base.py:200-227); LO kernels A_a, A_c, A_d,
I_D, I_ν rederived from `LoopIntegrand` with exact ω integration at 3 random rational points (ratios exactly 1); η matrix
of `kpz_notation.tex` = p→0 limit of the loop (Together of the difference = 0); float accuracy of F_ν with 1/p² (≤1.2e-14,
d = 1, 16 GL nodes, p ∈ 1e-3…100); RG-time and sign conventions.

| ID | Sev | Finding | Resolution |
|----|-----|---------|------------|
| E1 | blocker | Contract exported `g`, but the pipeline's constant is `alamk` (KPZAdimensional.wl:21,23); BMWflow g = alamk² | Export alamk; contract states g = alamk², adapter maps (user) — L9, §6 table |
| E2 | major | p = 0 rule pointed at KPZ-only code that divides by p² again (model_kpz.py:243-263); which rows carry 1/p² unstated | `P2_DIVIDED`; model-agnostic extrapolation; compare `self.Integral` after `Integral_upd` |
| E3 | major | "IBP leaves only f, f′" holds for KPZ only: removing f″ brings r″ via ∂sR; closedness test fails; BMWflow has no r″ | Scope = KPZ-like, abort `firstDerivativesOnly`, no IBP (user) — L5 |
| E4 | minor | DE `integrate` silently integrates mixed-momentum Δ wrong (π/l_q vs 2π/(l_q+l_Q)) | `FrequencyIntegration.singleMomentum`; fix FUTURE/FRGDESIGN text |
| E5 | minor | Ward check needs ∂/∂(iν) before ν → 0 | step 5 |
| E6 | minor | Functions tied to references in Ansatz, references chosen in DeriveDimensions | `Dimensions.functionReferences` |
| E7 | minor | KPZSetup has two shapes (rn, rd); contract had one | `SHAPES`, per-shape tuples |
| E8 | minor | F′(s) = d/d(p²) vs f′ = d/dq; measure exclusion undocumented | conversion in PythonExport; CONVENTIONS |
| E9 | minor | `"Pipeline"` option already asserted DE; no fingerprint baseline | Automatic + match check; step 0 |
| E10 | minor | n_η = n_f assumed | `etaPerFunction` |
| E11 | minor | `toJulia` keeps `^` (XOR in Python) | own printer |
| E12 | minor | 1e-12 norm undefined | per-function max-norm |
Nits: line ranges; `D` is a System symbol; usage example is NLO/evo2.

## Round 2 — 0 blocker, 2 major, 8 minor, 2 nit
Held: KPZ kernels from `LoopIntegrand` with numerical ω integral (30 digits) vs transcribed model_kpz.py:82-100 to 1e-20;
η rows pointwise vs `kpz_notation.tex` at d ∈ {1, 2, 5/2, 3} (no f″, r″); measure; β and g = alamk²; lagged η; 1e-12
budget (noise injection, N ≤ 20000 Euler steps, drift ≤ 3e-13).

| ID | Sev | Finding | Resolution |
|----|-----|---------|------------|
| E1 | major | Contract draft still had the pre-r1 wording (1/p² target mismatch, no P2_DIVIDED/SHAPES/norm, NLO example) | contract rewritten (§6) |
| E2 | major | Measure constants placed on both sides of the interface; v_{d−1} at d = 1 has Γ(0) | stated once: excluded, adapter applies |
| E3–E10 | minor | line refs; names/order; FRGDESIGN:193; noBeta/unusedEquation; f(p²/k²) relabelling; deg-2 per family; python env + `--require-python`; N unspecified | all folded in (steps 0, 4, 7, 8, 11, 12; §6) |
Nits: `rQ_` unused (∂sR only on the q line); `directionalAverage` already shared.

## Round 3 — 0 blocker, 1 major, 3 minor, 3 nit
Held: η rows symbolically; measure; β/g; tolerances (8.4e-15 d = 1, 5.8e-14 d = 2 on the usage-example grids);
`regularAtCoincidence` attainable.

| ID | Sev | Finding | Resolution |
|----|-----|---------|------------|
| E1 | major | No check that a divided row vanishes below p^{2n}; a shift-symmetry-breaking tadpole would export 1/p² silently | `ExternalMomentum.dividedOrders` (+ Protected low-p) |
| E2 | minor | angle-odd (`jac`/wedge) terms not rejected | `isotropicAngles` |
| E3 | minor | Evolution state is one scalar g; no map for other models | n_c = 1 + model table |
| E4 | minor | `validateStage[…,"Equations"]` in Adimensional/Oracle; Theory stamp | output stage "Equations"/"LO"; Theory unchanged |

## Round 4 — 1 blocker, 0 major, 5 minor, 2 nit
Held: line refs; residue rationality; linear part; η by Cramer vs `eta_calc`; β; measure; r_; dividedOrders expectFail
reachable (Γ4 tadpole of c4 gives a nonzero p⁰ term; (h,h) tadpole vanishes).

| ID | Sev | Finding | Resolution |
|----|-----|---------|------------|
| E1 | blocker | Fingerprint baseline impossible: nested `Timing` enters `contentHash` (probe: same Theory, different Ansatz fingerprint) | content baseline (step 0) |
| E2 | minor | η matrix orientation | row X, column Y |
| E3 | minor | stage name / what ExternalMomentum replaces | `BuildLOEquations`, §4 |
| E4 | minor | no model table / MODEL constant | `MODEL` + table |
| E5 | minor | stamps in Dimensions/Adimensional/Oracle; export guards | step 3 |
| E6 | minor | P2_DIVIDED ≥ 2 breaks tolerance | `p2DividedAtMostOne` |

## Round 5 — 1 blocker, 0 major, 2 minor, 1 nit
| ID | Sev | Finding | Resolution |
|----|-----|---------|------------|
| E1 | blocker | Baseline still not reproducible: `.jl` header has date + fingerprint; KS has no Julia export (r‴); `"Seconds"` nested in Equations | strip `Seconds`; hash `.jl` without header, KPZ/ShallowWater only. Probe: KPZ Equations `2306608ad3ff`, KPZ Adimensional `597b779e7005`, ShallowWater Equations `d3e8af22e17a`, KS Equations `3a053cfbad1f`, identical across kernels |
| E2 | minor | new Ansatz keys would flip the DE hash | keys only for LO |
| E3 | minor | nothing sets MODEL; no row check | `"Model"` option, `modelTable` |

## Round 6 — 0 blocker, 0 major, 4 minor, 2 nit
Held: baseline hashes identical across fresh kernels including one with a perturbed history ($ModuleNumber, Unique,
RandomInteger); all line refs; contract vs BMWflow code.

| ID | Finding | Resolution |
|----|---------|------------|
| E1 | `--require-python` not in the harness | step 12 |
| E2 | baseline targets ambiguous (ShallowWater LPA vs DE4; dims; test file) | step 0 |
| E3 | SHAPES order (KPZSetup declares {rn, rd}) | declared order, set comparison |
| E4 | DE-only stage functions lack pipeline guards | step 3 guards |
| E5 | lesson 01 / FRGDESIGN:152 say DE is the only pipeline | step 11 |

## Round 7 (the-enemy, final) — 0 blocker, 0 major, 3 minor
Held: all round-6 resolutions and line refs; η Cramer orientation; g map; measure; p = 0 rows; signs; KPZ ν-row of η has
only first derivatives (probe expanding model_kpz.py:95-100 to O(p²)).

| ID | Finding | Resolution |
|----|---------|------------|
| E2 | "∂s at fixed p gives p∂p f" is sign-ambiguous | stated at fixed p̂ (step 7) |
| E3 | `regularAtCoincidence` tested on which expression | on the Together-simplified residue sum |
| E4 | garbled "r‴" | reworded |

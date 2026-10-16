# FRGDESIGN — FRGscalarMSR design (hand-off for a fresh agent)

## 0. Read first

- `~/sync_cnrs/Projects/shallow_water/AGENT_PRIMER.md` — physics, conventions, traps, and the computer-algebra rules (§7) that are
  **binding** for this library.
- `~/sync_cnrs/Projects/shallow_water/DE4plus_fast.wl` — the optimised reference script (25 s, 9 equations, exact checks). Most
  numerical/algebraic machinery is reused from it (§9 below).
- `~/sources/FlowRG/docs/equations-interface.md` — the Julia contract the future JuliaExport module must target. (`~/sources/FlowRG`
  is the live path; `LPA_KPZ_thermal_forcing/CLAUDE.md` still says `~/Physics/Sources/FlowRG.jl` — stale, E13.)
- Adversarial review round 1 (the-enemy, 2026-10-16): findings E1–E14 are resolved in this plan; tags "E#" mark where.
  Probe scripts: `/tmp/claude-1002/-home-gtibone-sync-cnrs-Projects-shallow-water/ab3f0a62-a3fc-4e01-86ba-644acfd41efa/scratchpad/probe-*.wl`.
- Never edit `.nb` files or `eqs/*.jl` in the user's projects; the user regenerates them. Never kill `WolframKernel` processes you did
  not start. Kernel: `math -script` (Wolfram 15.0); `wolframscript` is not installed.

## 1. Context

The user (theoretical physicist, FRG / non-equilibrium field theory) derives Wetterich flow equations for stochastic PDEs in the
Martin–Siggia–Rose / Janssen–De Dominicis (MSR) formalism. Today one theory (stochastic shallow water, fields H, Ψ) and one ansatz
("DE4plus") are hardwired in a slow, fragile notebook and in `DE4plus_fast.wl` (hand-written Γ3 matrices, hand-written propagator,
hand-written compaction atoms, four placeholder ∂sR symbols).

**Goal:** a local Mathematica library, used from a notebook, that converts **an ansatz for Γ_k of N MSR scalars** into **a closed,
checked set of flow equations**. Guidelines: modularity, simplicity, self-consistency.

- Each module is a **grey box**: the user need not read its code, but must be able to understand its objective, mathematical
  assumptions and manipulations from `docs/MODULES.md` and from what it prints.
- Every user-facing function is explicitly named and readable/modifiable.
- Every internal step asserts and **fails loudly**.
- The user can inspect intermediate expressions at every stage with readable commands.

**Scope v1:** derivative-expansion pipeline (couplings depend on k only; projections at zero external frequency/momentum).
**Seams to keep open:** momentum-dependent pipeline (couplings depend on ω, p: external variables stay generic, richer pole
structure), ∂sΓ⁽³⁾ projections, symbolic Ward identities, Adimensional module, Julia export to FlowRG.

## 2. Settled design decisions (5 grilling rounds with the user)

| # | Decision |
|---|---|
| Q1 | N physical scalars, each with an auto-generated MSR response field; field order (φ1, φ̄1, φ2, φ̄2, …). Causality built in. |
| Q2 | Ansatz written as an action density (DSL); the library derives the Fourier vertices. **Forced printout** of every Fourier vertex with a convention header (user wants to check the Fourier transform by hand, incl. the δ-stripping convention). |
| Q3 | v1 projects only ∂sΓ⁽²⁾ (incl. the Γ⁽⁴⁾ tadpole when present). Vertex couplings are declared `NonRenormalised` with a reason. Slot kept for ∂sΓ⁽³⁾. |
| Q4 | Γ2 entries polynomial in iω and q²; deg_ω Δ ≤ 2 (asserted, abort beyond). Results rational in Δ's coefficients (Vieta), no radicals; pole-validity conditions emitted. |
| Q5 | Regulator declared separately; ∂sR derived; shape functions stay **abstract** r(y), r′, r″ (y = q²/k²) through to the Julia export, so shapes are swapped at run time to test stability. Amplitude absorbed into r (Q5b, no α symbol). |
| Q6/Q17 | Every (iω)^m (q²)^n monomial of every Γ2 entry up to explicit global `TruncationOrders` is classified: declared / `Protected` (flow must be exactly 0, else abort) / `Neglected` (flow computed, warning with its size). Untagged → abort. Rules (E5): monomials that vanish by structure (e.g. odd ω in a diagonal noise entry, Γ2_aa(ω)=Γ2_aa(−ω)) are detected and reported "vanishes by structure", no tag needed; a declared monomial with a **numeric** coefficient (e.g. the 1 in front of iω) must have zero flow (asserted); a monomial both declared-with-coupling and Protected → abort (contradiction). |
| E2 | Momentum order ≤ 1 (p²) asserted in v1: the script's Laplacian trick is wrong for (∇²_p)ⁿ, n ≥ 2 (gives 24d² instead of 8d(d+2)). Directional-average method for any order → FUTURE.md. Loop momentum q may carry any power. |
| E3/E4 | Loop prefactor (factor, sign, 2π) and Γ⁽⁴⁾ tadpole (−½) are **hand-derived** in MODULES.md, cross-read against textbook one-loop KPZ; no automated independent check (user choice). `flowsFAST` agreement does **not** verify the normalisation (same trace formula) — say so in docs. Cards print "tadpole term hand-derived, untested" whenever Γ⁽⁴⁾ ≠ 0. |
| Q7 / E1 | Fourier convention e^{−i(q·x−ωt)} (∂t → +iω, ∇ → −iq), leg a at (ω,q): reproduces `DE4plus_fast.wl` exactly (§3). Plane-wave cross-check of every vertex + an independent pole-placement check. |
| Q8 | One `.wl` per module, private contexts, immutable stage Associations passed explicitly, no global state. |
| Q9/Q14 | Standalone git repo `~/sources/FRGscalarMSR` (next to FlowRG). Acceptance: shallow-water DE4plus reproduces all 9 `flowsFAST` equations exactly + numerical oracle; one 1-scalar KPZ test. |
| Q10 | Output equation = ∂s(coefficient as written in the DSL, e.g. nk·Ak) = RHS. Interpretation (anomalous dimensions η, bookkeeping couplings = 1) deferred to the Adimensional module. |
| Q12 | Each stage prints a summary card (PASS/FAIL lines, LeafCount, timing); full expressions only through `Show…`. |
| Q13 | Optional `SaveStage`/`LoadStage`, fingerprinted on theory+ansatz+regulator+symmetries; mismatch aborts. |
| Q15 / E7 | DSL = density expression with **package-owned, protected** heads `dt, grad, div, lap` (not System `Dt/Grad/Div`: `Dt[nk H]` silently adds `Dt[nk]` terms). The parser always echoes the table of parsed monomials and asserts no System derivative head appears. |
| Q16 | Final `Equations` stage designed as the input of future Adimensional + JuliaExport → FlowRG `equations-interface.md`. |
| Q18 / E10 | Regulator written in Fourier space per field pair, values are **Γ2-level entries** (what adds to Γ2_ab, incl. the 2 of diagonal pairs), e.g. `{Hbar,H} -> nk q^2 arnh[q^2/k^2]`. ∂sR entries stay abstract symbols until Assemble; the map to the script's placeholders is dsNh = ∂sR_{HH̄}, dsMh = ∂sR_{H̄H̄}, dsNphi = ∂sR_{ΨΨ̄}, dsMphi = ∂sR_{Ψ̄Ψ̄}. |
| Q19 | `Protected[...]`, `Neglected[...]`, `NonRenormalised[...]`, each with a mandatory `"Reason"`, printed in every card. A banner warns that Ward identities are assumed, not derived (future work). |
| Q20 | Verb-first notebook workflow (§4). |
| Q21 | Failed check → red FAIL block (check, module, short offending expression, hint) then `Abort[]`; headless runs exit ≠ 0. |

The user prefers explanations that do not assume jargon (Vieta, Routh–Hurwitz, etc. must be explained in docs and printouts).

## 3. Conventions (`docs/CONVENTIONS.md`, printed as header of every vertex table)

- φ(t,x) = ∫ dω dᵈq/(2π)^{d+1} e^{−i(q·x − ωt)} φ(ω,q) ⇒ **∂t → +iω, ∇ → −iq** (E1: this is what reproduces the script's
  Γ_{HH̄}(ω) = iω + f_η; the opposite sign gives the time-mirror, poles in the wrong half plane).
- Functional derivative normalisation (E8): δφ(ω,q)/δφ(ω′,q′) = (2π)^{d+1} δ(ω−ω′) δ(q−q′), stated explicitly.
- δⁿΓ/δφ_{a1}(ω1,q1)…δφ_{an}(ωn,qn) = (2π)^{d+1} δ(Σω) δ(Σq) · Γ⁽ⁿ⁾_{a1…an}(ω1,q1; …; ω_{n−1},q_{n−1}); all momenta incoming,
  δ stripped, **last leg = −(sum of the others)**, fully symmetrised over identical fields.
- Γ2_ab(ω,q) ≡ Γ⁽²⁾_ab with leg a at (ω,q), leg b at (−ω,−q). Causal (retarded) poles lie in the **upper** half ω-plane; the script's
  `Iint` (closes upward) is reused unchanged.
- Vectors kept as bare symbols with `scal[a,b]` (Orderless, linear over **numbers** only — E14; extend to couplings or keep the wording).
- d is a symbol (`Dim`) throughout. ∂s ≡ k∂k. Momentum projection of p² = ∇²_p/(2d) at p = 0 (script trick, valid only for n = 1;
  n ≥ 2 aborts, E2); frequency projection of (iν)^m = (∂/∂(iν))^m / m! at ν = 0.
- Reserved internal symbols live in private contexts; `DefineTheory` aborts if a user name collides with a library symbol
  (primer trap: a global `s = …` silently broke the old rules; `k`, `l` were both heads and variables).

## 4. Notebook workflow (user-facing API)

```
Get["~/sources/FRGscalarMSR/FRGscalarMSR.wl"]
theory     = DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]       (* Hbar, Psibar auto *)
ansatz     = DefineAnsatz[theory,
  "Density" ->
      Hbar   ( dt[H]  - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]] )
    + Psibar ( dt[Psi] - nk ank lap[Psi] + gk H - g2k lap[H] + lk/2 grad[Psi].grad[Psi] )
    - Dhk aDhk grad[Hbar].grad[Hbar] - Dphik aDphik Psibar^2
    - Dphik D2phik grad[Psibar].grad[Psibar] + Kk grad[Hbar].grad[Psibar],
  "Couplings" -> {nk, Ak, Hk, lk, ...}]                                     (* echoes parsed monomials *)
regulator  = DefineRegulator[ansatz,
  {Hbar, H} -> nk q^2 arnh[q^2/k^2], {Psibar, Psi} -> nk q^2 arnphi[q^2/k^2],
  {Hbar, Hbar} -> -2 Dhk q^2 armh[q^2/k^2], {Psibar, Psibar} -> -2 Dphik armphi[q^2/k^2],
  "Shapes" -> {arnh, arnphi, armh, armphi}]
symmetries = DeclareSymmetries[ansatz,      (* complete tag set for shallow water (E5); precedence Protected > declared > Neglected *)
  Protected[{Hbar, _}, "Momentum" -> 0, "Reason" -> "Hbar -> Hbar + c(t) (mass conservation)"],
  Protected[{Psi, _},  "Momentum" -> 0, "Reason" -> "Psi -> Psi + c(t)"],
  Neglected[{H, Psibar},      "Frequency" -> {1, 2}, "Momentum" -> 0, "Reason" -> "truncation, z = 1 power counting"],
  Neglected[{Psibar, Psibar}, "Frequency" -> 2,      "Momentum" -> 0, "Reason" -> "coloured forcing, dropped"],
  Neglected[_,                "Frequency" -> {1, 2}, "Momentum" -> 1, "Reason" -> "mixed i w q^2 / w^2 q^2 terms, dropped"],
  NonRenormalised[lk, "Reason" -> "Galilean invariance"],
  "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]
(* the iw coefficient 1 in Gamma_{H Hbar}, Gamma_{Psi Psibar} is numeric: its flow is asserted 0;
   odd-w monomials of Gamma_{Hbar Hbar}, Gamma_{Psibar Psibar} are reported "vanish by structure" *)

vertices   = BuildVertices[ansatz]                      (* forced vertex table *)
prop       = BuildPropagator[vertices, regulator]
loops      = BuildLoopIntegrands[vertices, prop, regulator]
projected  = ProjectOnMonomials[loops, symmetries]      (* external ν, p → 0: the pipeline-1 step *)
integrated = IntegrateFrequencies[projected]
flows      = CheckClosure[integrated, symmetries]       (* protected = 0, neglected warned, LHS asserted *)
equations  = AssembleEquations[flows]                   (* compact atom form, FlowRG-oriented contract *)

ShowVertex[vertices, {Hbar, H, Psi}]   ShowGamma2[vertices]   ShowPropagator[prop, {H, Psibar}]
ShowIntegrand[projected, nk*Ak]        ShowEquation[equations, g2k]   ShowClosureReport[flows]
RunOracle[equations]                   SaveStage[equations, "file.wdx"]   LoadStage[...]
```

The quadratic signs above were rederived against primer §3 under the §3 convention (E9: `+Hk lk lap[Psi]` gives h = −Hk lk q² and
real waves ω² ≈ g·Hk·λ·q²; `+Kk grad.grad` gives Γ_{H̄Ψ̄} = +Kk q²). The cubic terms are checked against `AnsatzBMW` in milestone 2. Versus the round-5 preview the user approved,
projection is a separate stage *before* frequency integration, because in pipeline 1 setting ν = p = 0 is what makes the only
poles those of Δ(±ω); `CheckClosure` replaces the previewed `ProjectFlows`. Flag this to the user.

## 5. Repo layout

```
~/sources/FRGscalarMSR/                      (git init; README.md; no CLAUDE.md unless the user asks)
  FRGscalarMSR.wl            loader (Get of Kernel/*.wl in order), version, banner ("Ward identities assumed, not derived")
  Kernel/  Core.wl  Theory.wl  Ansatz.wl  Vertices.wl  Regulator.wl  Propagator.wl  Loops.wl  Projection.wl
           FrequencyIntegration.wl  Closure.wl  Assemble.wl  Oracle.wl  Stages.wl
  docs/    MODULES.md   (per module: objective · input/output stage · math assumptions · manipulations · checks · Show cmds)
           CONVENTIONS.md   FUTURE.md
  examples/ShallowWaterDE4plus.nb   examples/KPZ.nb
  tests/   runTests.wl (math -script, exit code)   unit/<Module>Test.wl
           reference/flowsFAST/*.wdx (copied from the shallow_water project)   ShallowWaterDE4plus.wl   KPZ.wl
```

## 6. Modules (objective · assumptions · manipulations · asserts)

Every stage object: `<|"Stage" -> name, "Pipeline" -> "DerivativeExpansion", "Fingerprint" -> …, "Upstream" -> {…}, data…,
"Checks" -> {…}|>`. Every constructor validates its input's `"Stage"` key and fingerprint.

0. **Core.wl** — `Assert[name, cond, details, hint]` (red block + `Abort[]`), `ExactlyZeroQ` (exact substitution of random
   rationals, several seeds; replaces `Simplify[…]==0`; E12: **`Dim` is randomised too**, and every shape atom r(y), r′(y), r″(y)
   is replaced by its own independent random rational), summary-card printer, fingerprint (`Hash` of canonical inputs), timing,
   LeafCount. Reuse `report`/`abort`/`randomPoint` patterns from the script.
1. **Theory.wl** `DefineTheory` — N scalars (assert N ∈ {1,2} while deg Δ ≤ 2), response fields `<name>bar`, field order, `Dim`.
   Asserts: fresh symbols, no collision with library/reserved symbols.
2. **Ansatz.wl** `DefineAnsatz` — parses the density: expand into monomials coefficient × Π(derivative operators on fields).
   Assumptions: couplings depend on k only; at most quartic in fields. Asserts: every symbol is a declared coupling or field; every
   monomial contains ≥ 1 response field (MSR); no unknown operator, no System `Dt/Grad/Div/Laplacian` (E7); couplings carry no
   t, x dependence; no `dt[coupling]`; no monomial of order p⁴ or higher (E2). Echoes the parsed table
   (coefficient | fields | operators). Seam: `"Pipeline"` option (only `"DerivativeExpansion"` accepted in v1).
3. **Vertices.wl** `BuildVertices` — Fourier transform per monomial; Γ⁽ⁿ⁾ as symmetrised functional derivatives; Γ2 matrix
   (polynomials in iω, q²) and Γ3/Γ4 functions in the leg convention (replaces the script's hand-written `Γ3[i]`). Asserts:
   (a) plane-wave cross-check — substitute Σ plane waves into each monomial independently and compare the coefficient with the
   derived vertex; implementation checks (E11, near-tautological, labelled as such): (b) Γ2_ab(ω,q) = Γ2_ba(−ω,−q); (c) Γ_φφ = 0;
   (d) Γ2 entries polynomial in iω and q². **Convention check (E1/E11)**: for each scalar, the free retarded response G_{φφ̄}
   built from the bare dt + damping terms has its pole in the upper half ω-plane. Forced printout of all nonzero vertices.
4. **Regulator.wl** `DefineRegulator` — R_k matrix from per-pair declarations; ∂sR by chain rule: ∂s(prefactor c) →
   `dsCoupling[c]` (symbol, becomes −η·c downstream), ∂s r(q²/k²) → −2y r′(y). Asserts: quadratic, ω-independent, consistent with
   the pair convention, only φ̄φ / φ̄φ̄ pairs, the declared shapes are the only free functions.
5. **Propagator.wl** `BuildPropagator` — works on abstract entry atoms f_ab(q²) (Γ_k + R per entry; primer rule 1 "abstract until
   last"); G = (Γ2+R)⁻¹ generated automatically (replaces `G2Implicit`/`proprules`); common denominator Δ(ω)Δ(−ω); Δ coefficients
   named (c0, c1, c2). Asserts: G·(Γ2+R) = 1 exactly; deg_ω Δ ≤ 2 (else abort "future extension"). Emits validity conditions:
   all roots of Δ(ω) in the upper half plane (degree 2 ⇒ the generalisation of the primer's A > 0, C > 0).
6. **Loops.wl** `BuildLoopIntegrands` — ∂sΓ2_ab(ν,p) = Tr[∂sR G Γ3_a G Γ3_b G] − ½ Tr[∂sR G Γ4_ab G] (signs/orderings in the §3
   convention), 1/(2π) of the ω-measure included, ∫dᵈq left to Julia. The formula, including the factor 1 (two orderings equal
   after a loop shift), the −½ and the 2π from E8, is **hand-derived step by step in MODULES.md** from ∂sΓ = ½ Tr[∂sR G]
   (E3/E4, user choice: no automated independent check). State plainly that `flowsFAST` and the oracle share this formula, so
   neither verifies it; cross-read against textbook one-loop KPZ. Card warning when Γ4 ≠ 0.
7. **Projection.wl** `ProjectOnMonomials` — for every classified monomial: ν-derivatives and p-Laplacians at 0 (generalises the
   script's `Dimensionalfactors`, `DerivativeCompletion`; asserts momentum order ≤ 1, E2); also projects Γ_k (regulator off)
   for the LHS. Pipeline-1 seam.
8. **FrequencyIntegration.wl** `IntegrateFrequencies` — `CoefficientRules` dispatch on ωⁿ Δ(ω)^−a Δ(−ω)^−b to memoised residues,
   generalised from the script's `Iint`/`ωIntegrate` to arbitrary Δ coefficients, degree 1 or 2. Asserts: convergence at large ω,
   no leftover ω, Δ, roots or radicals; result rational in Δ coefficients.
9. **Closure.wl** `CheckClosure` — LHS asserted = normalisation × coefficient (from Γ_k, regulator off); declared monomial →
   equation ∂s(coef) = RHS/norm; Protected → `ExactlyZeroQ` else abort with the reason; Neglected → value at a random point,
   warning; untagged → abort. `NonRenormalised` printed as an assumption. Must reproduce script §6b (iω, ω² in Γ_{HΨ̄}; ω² in
   Γ_{Ψ̄Ψ̄} are generated).
10. **Assemble.wl** `AssembleEquations` — compactify (script §5) with **auto-generated** scale-stripped atoms derived from the entry
   structure (replaces hand-written `toAtoms`/`atomVals`; the sync assertion stays as an exact check). Output contract per
   equation: loop q-integrand, local part, `dsCoupling[…]` symbols (linear), shapes r, r′, r″ abstract in y, validity conditions,
   `Dim` symbolic. Asserts no intermediate symbols left.
11. **Oracle.wl** `RunOracle` — script §6c generalised: numeric inverse of Γ2+R, NIntegrate over ω, finite differences in ν and an
   explicit p vector, run in **d = 2 and d = 3** (E12). Independent of Propagator/Projection/FrequencyIntegration; shares Γ3 and the
   loop formula (so it does not test E3).
12. **Stages.wl** — `SaveStage`/`LoadStage` with fingerprint check; `Show…` dispatchers.

## 7. Future seams (`docs/FUTURE.md`)

- **Pipeline 2** (couplings as functions of ω, p): Ansatz `"Pipeline"` option; replaces `ProjectOnMonomials` (external ν, p generic);
  the `IntegrateFrequencies` input contract is "rational in ω with listed pole factors", so shifted poles Δ(ω+ν, p+q) plug in.
- **deg Δ > 2**: general symmetric-function reduction; validity = Routh–Hurwitz conditions.
- **p⁴ and higher projections** (E2): write p = |p|e, Taylor-expand in |p|, average over directions with
  ⟨(e·q)²ᵐ⟩ = q²ᵐ (2m−1)!!/[d(d+2)⋯(d+2m−2)]; replaces the Laplacian trick for all orders; validate with the oracle in d = 2, 3.
- **Independent normalisation check** (E3/E4, declined for v1): background-field ½∂̃s Tr ln(Γ2[φ₀]+R) on KPZ + a toy quartic theory.
- **∂sΓ⁽³⁾ projections** (needs Γ5, triangle/box diagrams): new Loops variant; `NonRenormalised` becomes checkable.
- **Ward identities**: symbolic derivation/checks of symmetry constraints (currently user-declared).
- **Adimensional.wl / JuliaExport.wl**: map `Equations` onto FlowRG `equations-interface.md` (COUPLINGS/ETAS/SHAPES, positional args,
  β and η-system linear in η, l/nl split, shapes as `armh(aq^2)`, `Darmh(aq^2)`). FlowRG supports only first shape derivatives
  today; g2/D2phi need r″ — interface extension required.

## 8. Milestones

1. Core, Theory, Ansatz (parser + echo) + unit tests.
2. Vertices + plane-wave check; shallow-water Γ2/Γ3 equal the script's `AnsatzBMW` + `form*` exactly.
3. Regulator, Propagator (G·Γ = 1, Δ coefficients, validity conditions).
4. Loops, Projection, FrequencyIntegration: raw RHS of all 9 equations equal `flowsFAST` exactly, first with ∂sR entries as abstract
   symbols mapped to `dsNh, dsMh, dsNphi, dsMphi` (E10 map), then with the derived ∂sR substituted (random rationals, 30 digits).
5. Closure, Assemble, Oracle.
6. KPZ 1-scalar test, example notebooks, docs, README.

## 9. Reuse from `DE4plus_fast.wl`

`scal` rules (§1), `Iint`/`toXY`/`ωIntegrate` (§4), `Dimensionalfactors`/`DerivativeCompletion` (§3), `compactify`/`collectNum`/
`compactFrac` (§5), LHS assertion (`lhs`, `jobs` normalisations), `randomPoint` + regression (§6a), closure diagnostics (§6b),
oracle (§6c), `report`/`abort`. Primer §7 rules 1–7 are binding: abstract atoms until last, no radicals, Together/Collect/Factor
instead of Simplify, exact random-rational verification.

## 10. Verification

- `math -script ~/sources/FRGscalarMSR/tests/runTests.wl` → unit + acceptance tests, exit 0.
- Shallow water: all 9 equations equal `tests/reference/flowsFAST/*.wdx` exactly at random rationals; `RunOracle` agrees to ~1e-10;
  HEQ ≡ 0; runtime of the same order as the script (~25 s).
- KPZ (E6, restated so v1 can check it): dimensionful ν and D flows match a hand derivation (and the dimensionful content of the
  `LPA_KPZ_thermal_forcing` equations); the iω coefficient of Γ_{hh̄} at p = 0 is a Protected check (Galilean Ward identity) and
  must flow to exactly 0. λ non-renormalisation is **not** checkable in v1 (needs ∂sΓ⁽³⁾) — declared, printed as an assumption.
- Negative tests: a deliberately wrong vertex, an untagged monomial, a protected term with nonzero flow, a stale `LoadStage`, a
  degree-3 Δ, a p⁴ monomial, a System `Dt` in the density, and a declared+Protected contradiction each abort with a red FAIL block.
- Known limitation written in README: the loop normalisation is hand-derived, not independently tested (E3/E4).
- `examples/ShallowWaterDE4plus.nb` runs top to bottom; every stage prints its card; the vertex table is printed.

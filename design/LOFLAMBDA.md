# LOFLAMBDA — a function multiplying iω, and functions in vertices (KPZ f_λ at LO)

Decision record, grilled with the user (2026-11-03), attacked once (`design/reviews/LOFLambda-r1.md`, resolutions in §5).
Extends `design/LOPIPELINE.md`; supersedes its L1 and L2, and L7 for an iω entry that carries a function (§2 below).
Notation: ν is the viscosity coupling (nuk) and also the external frequency; I_X(p) is the loop part of ∂s f_X (the `Loop` row).
Line anchors refer to commit `6df4a91` (branch BMW-LO).

## 1. Goal
The Kloss–Canet–Wschebor LO ansatz for KPZ with a third function f_λ(p²):

    Γ_k = ∫ h̄ f_λ(−∇²) { ∂t h − λ/2 (∇h)² } − ν h̄ f_ν(−∇²) ∇²h − D h̄ f_D(−∇²) h̄

so Γ⁽¹¹⁾ = iω f_λ(p²) + ν p² f_ν(p²), Γ⁽⁰²⁾ = −2 D f_D(p²), and the cubic vertex is λ f_λ(k_h̄²)(k₁·k₂) with f_λ at the momentum of
the h̄ leg (Γ⁽¹¹⁾ here is Γ_{h̄h}(ω, p) with (ω, p) on the h̄ leg, the library's `VertexValue` convention) (p, q or Q = p + q inside the loop). Galilean invariance ties the f_λ of the cubic term to the f_λ of ∂t h. The functions
stay frequency independent (LO). The flow of f_λ(p²) is the iν coefficient of ∂sΓ⁽¹¹⁾(ν, p): ∂/∂(iν) taken **before** the ω
integral, at ν = 0, at generic p.

Acceptance (user): (1) in the analytic equations, f_λ ≡ 1 gives exactly the KPZLO equations of f_D and f_ν; (2) the exported
module integrated by BMWflow in d = 1, 2, 3; the exponents with and without f_λ are compared by the user (visual).

## 2. Decisions
| # | Decision | Reason |
|---|---|---|
| F1 | Syntax: `fn[flam, hbar] (dt[h] - lamk/2 grad[h].grad[h])`; `"Functions" -> {fD -> Dk, fnu -> nuk, flam -> 1}`. A function may sit in any monomial (supersedes L1), also times iω (supersedes L2). It must sit in exactly one bilinear monomial, with coefficient number × reference (that entry gives its flow); in a vertex monomial its coefficient is any constant. Checks: `Ansatz.functionInGamma2` (expectFail: fn only in a vertex), `Ansatz.oneFunctionPerCoefficient` (each (Γ2 entry, ω order) carries at most one function; KPZ: Γ(h̄,h) has f_λ at ω¹ and f_ν at ω⁰), `Ansatz.omegaFunctionReference` (a function times iω has reference 1; a running Z_k f_Z iω is out of scope, docs/FUTURE.md). | Fixes the "one function ↔ one flowing Γ2 coefficient" map; Vertices already puts f at the named leg's momentum (Vertices.wl:28). |
| F2 | Reference `1`: f(0) = 1 is not an η equation but a symmetry statement; no η, no η row, no reference coupling, `Dimensions.functionReferences` accepts it. The declared Protected iω tag at p = 0 (`ExternalMomentum.ward`) proves I_λ(0) = 0 exactly, i.e. f_λ(0) stays 1. A function with reference 1 and no matching Protected iω tag (or a ω⁰ entry) aborts: `ExternalMomentum.unitReference`. | Galilean Ward identity: the iω coefficient at p = 0 is not renormalised. |
| F3 | Flow of an iω function: `D[LoopIntegrand[…, nu, pV], nu] /. nu -> 0`, `deltaCompletion`, divide by I and by the entry's normalisation, then `loIntegrate` (the two-family residues already take pole order a_j ≥ 2, FrequencyIntegration.wl:75-91). P2_DIVIDED = 0. ν enters only G(Q) (the cubic vertex has no ω), so the double poles are in the Q family only. | Reuses the residue engine; ∂ν before integrating is what the user asked. |
| F4 | F(0) and F′(0): a vertex at the external leg gives F(p²); in the η rows (p → 0) it gives F(0) → 1 (normalisation of every function). `ExternalMomentum.etaAtZero`: no F′(0) or higher survives in an η row (KPZ: I_ν = f_λ(p) J(p), J = O(p²), so only f_λ(0)). | BMWflow's η functions get f, f′ at q only. |
| F5 | Contract v3: `loop(p, q, c, Q2, fp, fq, fQ, rq, rq_, rQ, eta, couplings, d)` (fp = f at p); `ETAS` entry `None` for reference-1 functions; `eta_matrix`/`eta_vector` only over functions with an η (n_η × n_η, in ETAS order without the Nones); `linear` of a reference-1 function = p f′; atoms f + r only for ω⁰ functions (`pyAtomPairs`). `CONTRACT_VERSION = 3`. `KPZ_LO.py` re-exported (same kernels, new signature). | f at p is needed for vertex functions; the η system has one row per η. |
| F6 | BMWflow side: a handoff (`docs/handoff/BMWflow_LO.md` v3), not edited here. | CLAUDE.md: other-repo changes go through `docs/handoff/`. |
| F7 | Benchmark `benchmarks/KPZLOLambda/`; `KPZLO` unchanged and is the reference of acceptance 1: `KPZLOLambda.flam1.*` substitutes f_λ → 1, f_λ′ → 0 in the Equations and Adimensional stages and compares fD, fν integrands, η matrix/vector, linear parts and β with KPZLO (`ExactlyZeroQ`). Python (numpy only): the new module at f_λ ≡ 1 equals `KPZ_LO.py` to 1e-13. | Acceptance 1, also through the printer. |
| F8 | Acceptance 2: `benchmarks/KPZLOLambda/runBMWflow.py` (BMWflow `.venv`, after BMWflow implements v3): both modules, d = 1, 2, 3, grids/regulator of `KPZ_1D_twogrids`; prints z = 2 − η_ν, χ = (2 − d + η_D − η_ν)/2 at the fixed point, PNGs in git-ignored `out/`. The user compares visually. | User. |
| F9 | Oracle LO: rows with frequency order 1 compare with the exact ∂ν (symbolic ν in the numeric 2×2 matrices) of the NIntegrate ω integral, against the integrand × normalisation (before the division by the normalisation), same tolerance 10⁻⁸ relative. | Checks the ω integration of F3 (double poles); the normalisation is checked by the DE route (§5 E1). |
| F10 | Independent references for the f_λ model (review E1, E2): the f_D, f_ν rows at generic f_λ equal the Kloss et al. kernels with f_λ slots (BMWflow `model_kpz.py`, `benchmarks/KPZLOLambda/Reference.wl`); the f_λ, f_ν, f_D rows at f_λ = 1 + x p², f_ν = f_D = 1 have p², p⁰ Taylor coefficients equal to the DE flows of x, ν, D of the ansatz (h̄ − x ∇²h̄)(∂t h − λ/2 (∇h)²) − ν h̄ ∇²h − D h̄² (DE projection, single-family residues, Assemble's normalisation: a different route); a sign flip of the f_λ row is rejected. | The f_λ row has no published form to compare with. |

## 3. Implementation (in `$FRGModules` order; every new check = named `assertCheck` + `expectFail`)
1. Ansatz (`loMonomialChecks`, Ansatz.wl:257-275; `"Functions"` parsing :167-188): F1, references `_Symbol | 1`.
2. Vertices (`loGamma2Check`, Vertices.wl:113): reference 1 already passes (number × 1 × F).
3. ExternalMomentum: `flowRec` reads the coefficient of the function in `VertexValue[v, pr, {{nu, pV}}]` per ω order (0 or 1),
   records `"Frequency"`; F3 for order 1; `"HasEta"` (reference ≠ 1); F2, F4 checks; EtaIntegrand only for HasEta rows.
4. Oracle: F9.
5. Dimensions: F2. Adimensional LO: η system over HasEta rows; linear = 2 ap² f′ for reference-1 functions; loops of
   reference-1 rows not divided by a reference.
6. PythonExport: F5.
7. Benchmark, handoff, driver, docs (MODULES, CONVENTIONS, FUTURE, LOPIPELINE note, CLAUDE.md map; BMWflow path
   `~/sources/BMWflow`).

## 4. Physics caveats
- d = 1: KPZ's time-reversal symmetry (f_D = f_ν relation, exact χ = 1/2, z = 3/2) holds for f_λ ≡ 1; ∫ h″ f_λ(−∂²)(h′)² ≠ 0 for
  f_λ ≠ 1 (review E2 probe), so the f_λ model breaks it and its d = 1 exponents are not protected. Acceptance 2 shows this;
  freezing f_λ = 1 in d = 1 is the KPZLO model.
- Causal poles: Δ(ω) = c₀ + c₁ iω with c₁ = f_λ(q²): the residues assume f_λ > 0 at every q and Q (stated on the ExternalMomentum
  `causalPoles` card; not checked at run time, like the other causal-pole conditions of LO, docs/FUTURE.md).
- The Ward check is structural for KPZ (the vertex with the external h leg is ∝ p, so the whole ∂sΓ_{h̄h} integrand vanishes at
  p = 0); it still guarantees f_λ(0) = 1 for any ansatz that passes it.

## 5. Review resolutions (`design/reviews/LOFLambda-r1.md`)
- E1 [blocker] fixed: F10's DE route checks the sign and normalisation of the f_λ row (`KPZLOLambda.de.flam`, sign flip rejected).
- E2 [major] fixed: F10's generic-f_λ reference (`KPZLOLambda.reference.D`, `.nu`); f_λ at p beyond p = 0 is tested by it.
- E3 [major] fixed: §4 (d = 1).
- E4 fixed: §4 (causal poles). E5 fixed: anchor. E6 fixed: F9. E7 fixed: `Ansatz.oneFunctionPerCoefficient`.
- E8 fixed: notation line above; `ExternalMomentum.wardNeglectedLO` is now emitted only for Protected iω entries without a
  function (supersedes L7 there). E9 fixed: §1.

## 6. Implementation notes
- `ExternalMomentum.etaAtZero` has no reachable expectFail: a vertex function at the external leg that survives p → 0 also
  breaks the Ward check, which runs first (tried: c3 f_λ h̄h² and c3 (f_λ − 1) h̄h²). Kept as a card check.
- Python round trip tolerance: 1e-13, or 10× Mathematica's own machine-precision error at that value (some points are
  cancelling sums: 3e-13 in machine precision for the f_D row at one point). `fp` is drawn last so the v2 points are unchanged.
- `ExternalMomentum` stage rows gain `"Frequency"` and `"HasEta"`; `EtaIntegrand` is `None` for reference-1 rows. The
  Adimensional LO stage gains `"HasEta"`; `"Etas"`/`"EtaNames"` list only functions with an η.

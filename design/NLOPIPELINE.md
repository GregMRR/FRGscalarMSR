# NLOPIPELINE — functions of (ϖ, p): the external frequency kept in the loop (KCW NLO)

Decision record, 2026-11-12. Attacked twice (`design/reviews/NLOPipeline-r1.md`, the-arch-enemy; `NLOPipeline-r2.md`,
the-enemy); resolutions in §5. Extends `design/LOPIPELINE.md` and `design/LOFLAMBDA.md`. Brief:
`~/sync_cnrs/Projects/Burgers/Burgers_NLO/NLO_IMPLEMENTATION.md`.
Notation: ν_k is the viscosity coupling; ϖ is the external frequency (the symbol `nu`/`loNu` in `LoopIntegrand`, `Global`w`
in stages, `aw` dimensionless). Line anchors refer to commit `33e790e` (branch main).

## 1. Goal
The NLO of Kloss, Canet, Wschebor, PRE 86, 051124 (2012), arXiv:1209.4650 ("KCW", Eqs. 43, 44, 49): the left-hand-side
functions depend on (ϖ, p); every function on the right-hand side (vertices, propagators) is at zero frequency (KCW Eq. 43); ϖ
enters the loop only through G(ω + ϖ, Q). For the KPZ/Burgers ansatz with f_D, f_ν, f_λ:

    Γ_k = ∫ h̄ f_λ(−∇²){∂t h − λ/2 (∇h)²} − ν h̄ f_ν(−∇²) ∇²h − D h̄ f_D(−∇²) h̄

and Γ⁽¹¹⁾ = iϖ f_λ + ν_k p² f_ν, Γ⁽⁰²⁾ = −2 D_k f_D (KCW Eq. 38). The LO pipeline is the ϖ = 0 slice of this one, and that
slice is closed: η, g and f_X(p, 0) of the NLO equal those of LO.

Before this record, `BuildLOEquations` removed ϖ before the ω integral (`ExternalMomentum.wl:90-91`), `loIntegrate` aborted
on Δ(±(w + ϖ), sQ) (`FrequencyIntegration.wl:93-106`), `Adimensional` had no ϖ ∂ϖ term and the BMWflow contract (v3) had no ϖ
argument.

## 2. Decisions
| # | Decision | Reason |
|---|---|---|
| N1 | Residues with a shift parameter per family: `FrequencyIntegralShift[j]`, family spec `{a, b, deg, "Shifted"}`. Factors Δ_j(w + s) = lead Π(w + s − z), Δ_j(−w − s) = lead (−1)^deg Π(w + s + z); residues of the X factors at w = z − s. A real s moves no root across the real axis. | ϖ stays one symbol; unshifted memoised results unchanged; degree 2 needs no change to `loEliminateBeta` (both roots carry the same shift). |
| N2 | NLO is an option of the existing constructor: `BuildLOEquations[loops, symmetries, "ExternalFrequency" -> True]` (default `False`). The LO code path is unchanged. Only an NLO stage gets the keys `"Approximation" -> "NLO"` and `"ExternalFrequency" -> Global`w`; an LO stage has no new key, so its fingerprint (`Core.wl:102-105` hashes all keys) is unchanged. Downstream reads `Lookup[stage, "Approximation", "LO"]`. `"Pipeline"` stays `"LO"`. | All LO checks (zero entries, ward, divided orders, η rows) are reused at ϖ = 0; no duplicate constructor; the LO fingerprints stay a regression test. |
| N3 | Rows: f_D = F_{h̄h̄}(ϖ, p)/(−2D_k); f_ν = E[F_{h̄h}]/(ν_k p²); f_λ = O[F_{h̄h}]/(iϖ), with E[F] = [F(ϖ) + F(−ϖ)]/2, O = [F(ϖ) − F(−ϖ)]/2, built by `nu -> -nu`. In general: a (iϖ)⁰ function's row is E[F]/(norm p^(2n)), a (iϖ)¹ function's row is O[F]/(iϖ norm p^(2n)). An entry with a function at (iϖ)⁰ and none at (iϖ)¹ (KPZ with f_λ ≡ 1) drops O[F] (INFO `ExternalMomentum.wardNeglectedNLO`). | KCW Eq. 38, 49. For real ϖ and real functions Γ(−ϖ) = Γ(ϖ)*, so E is the real part and O is i × the imaginary part. |
| N4 | The (iϖ)¹ row: ϖ cancelled exactly per term of `loIntegrate`, before the `cv` substitution. For a term t = c(ϖ) R(ϖ) (c from the numerator, R the residue result), O[t]/(iϖ) = Together[(E[c] O[R] + O[c] E[R])/(iϖ)]; E[R] and O[R]/(iϖ) are memoised per exponent pattern, E[c], O[c] are computed per term. One exported expression, finite at ϖ = 0. Check `ExternalMomentum.nloLimit`: its ϖ → 0 value equals the LO ∂/∂(iϖ) row. | BMWflow evaluates the loop on a ϖ grid whose first point is 0 (`grid_log0`). Multi-field numerators carry (w + ϖ), so c depends on ϖ in general. |
| N5 | External-frequency symbol `Global`w` (dimensionful), `aw` (dimensionless); w = aw × (frequency unit of `DeriveDimensions`). Linear part: −(∂s ln of the frequency unit) aw ∂aw f, derived from the dimension vectors, not typed (= (2 − η_ν) aw ∂aw f for Burgers/KPZ, KCW Eq. 56). | Same derivation as the p ∂p term. |
| N6 | Contract v4 = v3 + the NLO kind. `ExportPython` writes `CONTRACT_VERSION = 4` for both kinds; LO signatures unchanged. NLO: `APPROXIMATION = "NLO"`, `loop(p, w, q, c, Q2, fp, fq, fQ, rq, rq_, rQ, eta, couplings, d)` (fp, fq, fQ at ϖ = 0), `linear(p, w, f, fprime, fprime_w, eta, couplings, d)` (fprime_w = ∂f/∂aw), `eta_matrix`, `eta_vector`, `beta` as v3. BMWflow accepts v3 (LO) and v4 (LO, NLO). | Existing v3 modules (e.g. `Burgers_LOplus/eqs/Burgers_LOplus.py`, user-owned) keep working without regeneration. |
| N7 | New checks at ϖ ≠ 0, each a named `assertCheck` with an `expectFail` test: `FrequencyIntegration.shiftToZero`, `ExternalMomentum.nloLimit`, `ExternalMomentum.nloNoiseEven` (O of a (iϖ)⁰-only entry such as F_{h̄h̄} vanishes), `ExternalMomentum.nloWardAllFrequencies` (an (iϖ)¹ row with a Protected iω tag at p = 0 vanishes at p = 0 for every ϖ, KCW Eq. 39), `ExternalMomentum.nloDividedOrders` (the orders below p^(2n) of a divided (iϖ)⁰ row vanish at every ϖ). Route for the p = 0 checks: the integrated NLO row (ϖ symbolic) at p = 0 (then Q² = q², finite for ϖ ≠ 0), `ExactlyZeroQ`; not `loProjectedAtZero`, whose `project` sets `nu -> 0` (`Projection.wl:86`). The `expectFail` of an all-ϖ check breaks the identity only at ϖ ≠ 0; the `expectFail` of `nloLimit` and `shiftToZero` (compared at ϖ = 0) uses a ϖ-independent error. | `nloDividedOrders` is required for the 1/p² division of the f_ν row at ϖ ≠ 0 (time-gauged shift symmetry); the brief did not list it. |
| N8 | The η system, β and the ϖ = 0 slice stay as LO. BMWflow fills p = 0 of a divided row per ϖ column, from that column (`_fill_p0`). | The brief (§2) says the p = 0 fill uses ϖ = 0 data only: wrong for ϖ ≠ 0 columns. |
| N9 | Oracle: an NLO stage is compared with the NIntegrate ω integral of the trace formula at ϖ ≠ 0, E and O parts before the normalisation, relative 1e-8; the ϖ = 0 comparisons stay. | Independent of the residues. |

## 3. Implementation (in `$FRGModules` order; every new check = named `assertCheck` + `expectFail`)
0. Before any change: record the fingerprints of the KPZLO and KPZLOLambda Equations and Adimensional stages; after the
   change they must be identical (N2).
1. FrequencyIntegration: N1 in `FrequencyIntegralFamilies`; `loCheckRegularAtCoincidence` substitutes
   `FrequencyIntegralC[2, m] -> FrequencyIntegralC[1, m]` and `FrequencyIntegralShift[2] -> 0` together;
   `loIntegrate[e, w, {s1, s2}, deg, nu]` maps Δ(w + ν, s2), Δ(−w − ν, s2) to the shifted family (Mathematica evaluates
   `-(nu + w)` to `-nu - w`: one pattern). Single-family terms at Q: substitute w → w − ν in the whole polynomial numerator,
   re-expand in powers of w, then `FrequencyIntegral[n', …]` per power (the shift-invariant integral); `FrequencyIntegral`
   unchanged. `FrequencyIntegration.ready` still aborts on any other argument.
2. ExternalMomentum: option `"ExternalFrequency"`; N3, N4, N7.
3. Oracle: N9 (`oracleLO`; `Oracle.wl:82` sets ν = 0 today).
4. Adimensional (LO branch), gated on `"Approximation" === "NLO"`: w → aw × frequency unit; functions f(ap², aw); N5;
   check `Adimensional.nloSlice` (aw → 0 rows equal the LO stage's rows).
5. PythonExport: N6; `docs/handoff/BMWflow_LO.md` v4.
6. Tests, `benchmarks/KPZNLO/` (KPZ with f_λ ≡ 1: f_D, f_ν rows equal BMWflow's NLO kernels A_a, A_c, A_d of
   `notes/NLO_flow_equations.tex`; KPZ with flowing f_λ: N7; exports `KPZ_NLO.py`, `KPZLambda_NLO.py`), Python round trip
   (1e-13, includes w = 0), docs.
7. BMWflow (`flow/models/model_generated.py`): N6 version table; `loop_shape` includes the ϖ axis in NLO;
   `Integral_pfixed_dD_NLO` (fp at the ϖ = 0 column); `f_rhs_calc` passes w and ∂f/∂aw from `f_spl_w`; tests with
   f_λ ≡ 1 against `ModelKPZ` NLO (relative 1e-10) and the ϖ = 0 column against the LO module (relative 1e-10).
   Order: from step 5 until BMWflow reads v4 and the fixtures are re-vendored, `tests/acceptance/BMWflowLO.wl` fails; run
   the full suite with `--require-python` after step 7.

## 4. Physics caveats
- d = 1: KCW fix f_λ ≡ 1 to keep time reversal (KCW Eq. 42); a flowing f_λ breaks it (KCW Appendix C, `LOFLAMBDA.md` §4).
- Causal poles: the residues assume f_λ(q, 0) > 0 at every momentum; nothing checks it at run time (docs/STATUS.md).
- Cost: an NLO BMWflow run is about N_ϖ times an LO run.
- KCW Table I: χ moves from their LO to their NLO by 0.04 in d = 2 and 0.006 in d = 3; KCW consider NLO reliable up to
  d ≈ 3.5. KCW LO = BMWflow `"NLO"` with f_λ ≡ 1 (`Burgers_NLO/LO_vs_NLO_KPZ.pdf`).

## 5. Review resolutions
Round 1 (the-arch-enemy, `reviews/NLOPipeline-r1.md`: 0 blocker, 3 major, 6 minor):
- E1 [major] fixed: N2 adds keys only to NLO stages; §3 step 0 records the LO baseline before any change.
- E2 [major] fixed: N7 route (integrated NLO row at p = 0, not `project`) and ϖ-only breakers.
- E3 [major] fixed: §3 step 7 (`loop_shape`).
- E4 fixed: §3 order note. E5 fixed: N4 (E[c]O[R] + O[c]E[R]). E6 fixed: N3 (`wardNeglectedNLO`).
- E7 fixed: the d = 1 literature target applies to f_λ ≡ 1 only (handoff `Burgers_NLO/RUN_NLO.md`).
- E8 fixed: ϖ = 0 column tolerance 1e-10 (different algebraic forms). E9 fixed: notation line.
- Nits fixed: `flow/create.py` docs; Oracle point count; one pattern for `-(nu + w)`.
Round 2 (the-enemy, `reviews/NLOPipeline-r2.md`: 0 blocker, 1 major, 3 minor):
- E1 [major] fixed: N7, ϖ-independent breakers for `nloLimit` and `shiftToZero`.
- E2 fixed: §3 step 1 (numerator shift before `FrequencyIntegral`). E3 fixed: §3 step 4 (gate).
- E4 fixed: the BMWflow test compares both modules with Mathematica values at random points.

## 6. Implementation notes
- **Fingerprints are not reproducible between kernel sessions** (two runs of unchanged `main` give different hashes), so §3
  step 0 as written could not work. Replacement: the LO Equations rows and the Adimensional parts (Loop, Linear, EtaMatrix,
  EtaVector, BetaL, Functions, Etas, HasEta, P2Divided) and the check names of KPZLO and KPZLOLambda, built by `main` and by this
  branch, compared with `===`: identical.
- `Together` keeps conjugate denominator factors (L ± iϖ f_λ f_λ)ⁿ apart, so real rows carried `I`; `nloRealDenominator`
  multiplies each pair into (L² + ϖ² f_λ² f_λ²)ⁿ, and `ExternalMomentum.nloReal` asserts that no `I` remains.
- The N7 checks, `nloP2Divided`, `nloReal` and `Adimensional.nloSlice` are helpers: their expectFail tests feed synthetic rows
  (e.g. w² q², zero at w = 0). `nloNoiseEven` also has a real-ansatz expectFail (KPZ without the Protected iω tag).
- `ExternalMomentum.nloZeroEntries` was added (entries outside the ansatz vanish at every ϖ; causality of Γ_{hh}).
- Cost (KPZ with f_λ): NLO Equations 6 s, Adimensional 9 s, export 5 s; LeafCount of the f_ν row 2.6·10⁵ (LO: 2.3·10⁴); module
  37 kB (LO 15 kB).
- Independent results: the f_D, f_ν rows at generic f_λ and ϖ equal BMWflow's hand-written NLO kernels exactly; the generated
  NLO module at f_λ ≡ 1 equals `ModelKPZ` NLO to 1.5·10⁻¹³ (Integral and f_rhs, d = 1, 2).
- Reference points: at w = 3.375 the f_D kernel value 1.7·10⁻⁷ is a cancelling sum; Python and Mathematica differ by 2·10⁻¹¹
  relative there, so the BMWflow point test uses 10⁻¹⁰.
- The Oracle needs no separate unit test for NLO: `KPZNLO.oracle` runs it (24 comparisons, d = 2, 3).

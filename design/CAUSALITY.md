# CAUSALITY — causal-pole conditions, dimensionless and minimal, in the FlowRG file

Decision record. Grilled with the user, then attacked once: `design/reviews/Causality-r1.md` (1 blocker, 3 majors, 5 minors,
all resolved below, §5). Line anchors refer to commit `1b9fae1`.

## 1. Goal
`BuildPropagator` derives the conditions under which every root of Δ(ω) lies in the upper half ω-plane (causal poles); every
frequency integral assumes them. Until now they were only printed as an ASSUMPTION (`Propagator.validity`,
`Assemble.validity`) and dropped by `AdimensionalEquations`. They are now carried through Adimensional in a minimal
dimensionless form and written into the equations file as an extra function, so that whoever runs the flow can evaluate them.

- **Export only.** FlowRG is not changed. An extra definition is allowed by rule 6 of
  `FlowRG.jl/docs/equations-interface.md`. A FlowRG-side check is in `docs/FUTURE.md`.

## 2. Decisions
| # | Decision | Reason |
|---|---|---|
| C1 | Scope: exactly the propagator's validity conditions (degree 2: c1/c2 > 0, c0/c2 > 0; degree 1: c0/c1 > 0). | No new physics; the criterion is re-derived in r1 (u = iω, Im ω > 0 ⇔ Re u < 0). |
| C2 | `AdimensionalEquations` computes them (stage key `"CausalPoles"`); `ExportJulia` only writes. | Adimensional owns the unit substitution, q → k aq and the atom compaction. |
| C3 | Assumption: every reference coupling > 0, k > 0, aq > 0. ASSUMPTION line `Adimensional.referencesPositive`. | Needed to strip their powers; already implicit (fractional powers of references, ∂s ref = −η ref keeps the sign). |
| C4 | Minimal = clear denominators, strip manifestly positive factors (positive numbers, powers of references, k, aq), strip common positive content. Not: shape positivity, implication between conditions. | Exact and independent of the regulator: shapes are bound in Julia, so `1 + rn(aq²) > 0` stays. |
| C5 | Julia: `causal_poles(couplings…, aq) = (e1, …)`, each entry must be > 0 at every aq. | Positional like the other functions; `aq` always last even if unused. |
| C6 | The contract is the **sign only**. | After stripping, the scale of each entry is arbitrary (r1 E9): a "margin" means nothing. |
| C7 | Names `CausalPoles` / `causal_poles` / `Adimensional.causalPoles*` / `JuliaExport.causalPoles`. | "Causality" already names the Γ_φφ = 0 class (r1 E8). |
| C8 | A condition that reduces to a positive constant is dropped in Adimensional; to 0 or a negative constant: abort `Adimensional.causalPoles`. None left: `causal_poles` is omitted, INFO "trivially satisfied". | Same parts in the stage, `JuliaParts` and the Julia round trip (r1 E6). An acausal Δ at every state is an error. |
| C9 | Follows `"Form"`: Compact (dimensionless atoms as locals) or Explicit. | As `betanl`. |

## 3. Algorithm (per condition `e > 0`, compact form of `equations["ValidityConditions"]["Atoms"]`)
1. **Clear denominators.** `{n, d} = NumeratorDenominator[Together[e]]`; if d is a number, `p = Sign[d] n`, else `p = n d`
   (same sign as n/d).
2. **Units.** Substitute the dimensionless atoms (`ScaleAtom → [atom] AdimAtom`), then `toUnits` (non-references → name × form,
   q → k aq). `dimensionVector` splits p = (basis monomial, positive by C3) × rest; rest is free of k and the references.
3. **Content.** On the expanded rest (a small polynomial in atoms, couplings, aq; never the explicit form): the positive numeric
   content `GCD` of the term coefficients, and `aq^m` with m the smallest aq exponent over the terms, counted outside function
   arguments (`rn[aq^2]` does not count). Minimal = rest / (content); the stripped factor is basis monomial × content.
4. **Explicit.** Minimal compact with the `AdimAtom` values substituted. Constant explicit → C8.

**Check `Adimensional.causalPolesMinimal`**: p − stripped × minimal = 0 exactly (`ExactlyZeroQ`), and stripped is a product of a
positive number and powers of the basis units and aq. The stripping is done by the private helper `causalPolesStrip`; the unit
tests replace it with `Block` (as `tests/unit/ProjectionTest.wl:60`) by a wrong factor and by a negative content (r1 E3).

## 4. Results
- KPZ: Compact `{AdimAtom[h, hbar, 0, 0]}`, Explicit `{1 + rn[aq^2]}`.
- Shallow water (full, `ShallowWaterSetup.wl`): `Ak + ank + arnh + arnphi`, `agk aHk alk + aq² (ag2k aHk alk + (Ak+arnh)(ank+arnphi))`.
- Shallow water LPA: `1 + Ak + aHk alk amk + arnh + arnphi`, `agk aHk alk + aq² (Ak+arnh)(1+arnphi)`.

## 5. Review r1 resolutions
- E1 (blocker), E2: expected values above taken from a dimension count and the probe, each in the file whose ansatz produces it.
- E3: the check is tested by mutating `causalPolesStrip` (§3).
- E4: one wording everywhere: "> 0 at every aq" (the cards' "at every q"). How to check it (nodes, endpoints, sampling between
  nodes) is not decided here: the function is only exported (`docs/FUTURE.md`).
- E5: numeric content is the positive GCD; the aq power is taken separately, outside function arguments.
- E6: C8. E7: §4 names both forms. E8: C7. E9: C6.

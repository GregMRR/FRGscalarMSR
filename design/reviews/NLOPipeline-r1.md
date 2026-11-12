# the-arch-enemy — NLO plan (now design/NLOPIPELINE.md), round 1
Attacked commit: 33e790e (plan file, before the record) · Date: 2026-11-12 · Class: plan
Verdict: 0 blocker, 3 major, 6 minor, 3 nit. Report returned inline (plan mode); copied here. Resolutions: NLOPIPELINE.md §5.

## Held
- File:line claims (ExternalMomentum.wl:90-91, FrequencyIntegration.wl:93-106, Adimensional.wl:357, model_generated.py:48,175-181).
- N1 shifted-family residues: factorisation checked by hand for deg 1, 2; Mathematica probe on w X1² Y1 X2 Y2 (shifted) agrees
  with NIntegrate to 15 digits at ν = 3/5 and ν = −11/7. regularAtCoincidence with C2→C1 and shift→0 together is valid.
- N3: O/(iν) at ν → 0 equals the LO ∂/∂(iν) row; the h̄h̄ bubble is even in ν at fixed (q, c).
- N4 for KPZ; N5 (∂s ln(ν_k k²) = 2 − η_ν, sign as model_kpz.py:179); N8 (`_fill_p0` per column).

## Findings
- E1 [major] The LO fingerprint test had no baseline, and the new stage key changed the LO fingerprint (Core.wl:102-105).
- E2 [major] The N7 "all ν" checks had no route: `loProjectedAtZero` → `project` sets nu → 0 (Projection.wl:86).
- E3 [major] BMWflow `loop_shape` is (Np, 1, q, θ): broadcast fails for N_w > 1.
- E4 [minor] BMWflowLO.wl fails between the v4 export and the BMWflow re-vendor.
- E5 [minor] Per-pattern memo of O assumes c free of ν; O[cR] = E[c]O[R] + O[c]E[R].
- E6 [minor] O[F_{h̄h}] dropped silently when f_λ ≡ 1.
- E7 [minor] d = 1 χ = 1/2 is not a target for a flowing f_λ.
- E8 [minor, suspected] 1e-12 between two algebraic forms is too tight.
- E9 [minor] ν used for viscosity and frequency.
- Nits: flow/create.py "LO only" docs; Oracle point count; redundant `-(nu + w)` pattern.

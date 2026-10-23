# the-enemy — design MOMENTUMFUNCTIONS (plan read-claude-md-i-greedy-platypus.md), round 1
Attacked commit: 21719c9 (plan file outside the repo)  ·  Date: 2026-10-23
Verdict: 2 blocker, 5 major, 6 minor, 1 nit

(Plan mode prevented writing design/reviews/MOMENTUMFUNCTIONS/r1.md; this file holds the report verbatim. Save it there, with
the probes of the appendix as r1-probe-*.wl.)

Confirm green: `math -script tests/runTests.wl` → 486/486 passed, 2 SKIPPED (7 min 56 s); tree clean afterwards.

## Attack surface
- D1 / Ansatz.functionQuadraticOnly (plan:11,40) — nonsense-ansatz hunt — held (function in Γ3/Γ4 named)
- D2 frequency polynomial (plan:12) — consistency w/ FrequencyIntegration — held
- D3 mixed ansatz, p→0 for constant slots (plan:13) — consistency w/ D10 — finding E12
- D4 JuliaExport abort (plan:15) — held
- D5 DSL, "F(p²) unambiguous" (plan:16-17) — write-order/placement ambiguity: Times is Orderless, F(−∇²) self-adjoint on a quadratic form → held; deferred Couplings-list decision → E6; applied symbol `Fnu[0]` → E8; slot multiplicity → E5
- D6 prefactor rule (plan:18-20) — degeneracy hunt — E4
- D7 normalisation, 2ŷ∂ŷF̂ (plan:21-23, 67) — independent rederivation (probe-linear) — sign held; F̂·loop(0) subtraction held by hand; ŷ=0 evaluation → E1; option misuse → E9
- D8 one function per slot (plan:24-25) — ambiguity hunt — E4, E2
- D9 measure (plan:26-27) — identity S_{d−1}∫(1−c²)^{(d−3)/2}dc = S_d verified (probe-linear) → held; K_d, d=1, mixed measures → E10
- D10 constraints (plan:28-32) — semantics vs Symmetries.wl:79-80, KS Checks.wl:18-19 — E7, E12
- D11 Oracle (plan:33) — shares vertices/loop formula (existing caveat) — held
- D12 benchmark (plan:34-36) — reference integrity — E13
- Vertices cascade (plan:44-46) — check integrity — E11
- Symmetries/Projection cascade (plan:50-56) — division by p^{2n₀} — E1, E2, E11
- FrequencyIntegration cascade (plan:57-60) — "s₁ ≠ s₂ for p ≠ 0" — E3 (probe-twofamily)
- Closure/Adimensional cascade (plan:61-69) — E1, E7, E11
- Regulator "already aborts" (plan:47) — vs Regulator.wl:69-71 — E6

## Findings

### E1 [blocker] The functional integrand divided by p^{2n₀} is pointwise O(1/ap) as ap→0; ŷ = 0 (D7) and the "exact" p→0 regression are undefined as specified
- Where: plan:22-23 (D7 "η from Fnu's equation at ŷ = 0"), plan:26-27 (D9 stored integrand), plan:53-55 (Projection: divide by p^{2n₀}, no direction average), plan:81-82 (regression); symmetrisation appears only for constraint checks (plan:29-30, 62).
- Evidence (r1-probe-pointwise.wl): KPZ (hbar,h) bubble from the library's own `LoopIntegrand`, d = 1, aq = 7/10, c = +1, ν = 0, ω integrated numerically, with a momentum-dependent ν: GammaAtom[h,hbar,0](s) = s/(1+s) + s e^{−s}:
  p = 1e-1, 1e-2, 1e-3, 1e-4 → raw I/p² = 0.194, 3.20, 33.5, 336.5 (∝ 1/p); symmetrised (I(q)+I(−p−q))/2/p² = −0.2530, −0.2391, −0.2379, −0.2378 (finite).
  With constant ν (F ≡ 1) the raw ratio is finite (−0.0942 … −0.1178): constant-coupling KPZ tests cannot see this.
  Reason: the O(p) term is p·V(q) with V odd; it dies only under the angular integral or under q ↔ −p−q.
- Consequence: F̂ evaluated at ŷ = 0 gives 0/0 or ±∞ depending on c. The Normalise η row is undefined. The "ap → 0 equals today's KPZ exactly" regression cannot be done pointwise. The downstream integrator cancels ±1/ap terms between c and −c, so precision is lost at small ap.
- Hint: store the q ↔ −p−q symmetrised integrand. Assert, as a named check, that it is O(p^{2n₀}) pointwise. Define ŷ = 0 as a series limit, and do the c-integral (⟨c²⟩ = 1/d) in the regression.
- Resolution: _(author)_

### E2 [blocker] Dividing by p^{2n₀} assumes the lower Taylor orders of the flow vanish. Nothing classifies or checks them, and existing tags lose their meaning
- Where: plan:24-25 (D8: no other monomial allowed in the slot), plan:50-51 ("truncation applies to function-free slots only"), plan:54 (divide by p^{2n₀}); Symmetries.wl:79-80 (`Symmetries.contradiction`: a Protected slot that carries a coupling aborts); tests/reference/KPZSetup.wl (`Protected[{h,hbar},"Frequency"->0,"Momentum"->0]`, the shift symmetry).
- Evidence: in KPZ_LO the slot ((hbar,h), iω⁰) is functional (Fnu, n₀ = 1). The p⁰ part of its flow is then (a) not a classified monomial, (b) not allowed as an ansatz term (D8), and (c) not expressible as a tag, because D10 defines Protected on a functional slot as "the flow vanishes ∀p", which is false for KPZ. The `"Momentum" -> 0` spec has no meaning on a functional slot. Take a model with no shift symmetry (any non-derivative vertex such as `gk hbar h^2`): its generated p⁰ term m² goes into F̂ as m²/p², silently, and no named check fires.
- Consequence: a β_F̂ that is divergent or meaningless at small ŷ, accepted silently.
- Hint: classify orders p⁰…p^{2n₀−2} of every functional slot with today's p→0 machinery (each must be Protected or numeric-zero) and assert the O(p^{2n₀}) behaviour. Define what "Momentum" means for tags on functional slots.
- Resolution: _(author)_

### E3 [major] "Generic s₁ ≠ s₂ assumed (p ≠ 0)" is false: s₁ = s₂ on the surface c = −ap/(2aq), and the per-family residues blow up there
- Where: plan:57-60.
- Evidence (r1-probe-twofamily.wl): Δᵢ = aᵢ + iω and f = 1/(Δ₁(ω)²Δ₁(−ω)Δ₂(ω)Δ₂(−ω)). Family 1 gives −π(5a₁² − a₂²)/(2a₁²(a₁−a₂)²(a₁+a₂)²) and family 2 gives π/((a₁−a₂)²a₂(a₁+a₂)): double poles at a₁ = a₂. Their sum π(2a₁+a₂)/(2a₁²a₂(a₁+a₂)²) is finite, and NIntegrate agrees (0.036651914). Solve[|p+q|² = q²] gives c = −ap/(2aq), which lies in [−1, 1] whenever aq ≥ ap/2, i.e. inside the integration domain for every ap. In degree 2 the spurious factor is Res(Δ(·,s₁), Δ(·,s₂)).
- Consequence: unless the two families are combined and cancelled, the stored integrand is 0/0 on a hypersurface of every (aq, c) integral. NIntegrate unit tests at generic points will not catch it.
- Hint: combine the families before storing. Add a named check that the result is finite at coincident atoms (s₂ → s₁), and a NIntegrate unit test at s₁ = s₂.
- Resolution: _(author)_

### E4 [major] The form of a functional slot's coefficient is not pinned down: several ansätze are degenerate or ambiguous
- Where: plan:18-20 (D6), plan:24-25 (D8).
- Evidence:
  - (a) `-nk hbar lap[h] - nk Fnu hbar lap[h]` gives the slot nk(1+F)p². D8's "any other p^{2n} monomial" does not say whether a constant at the same power counts. The parametrisation is degenerate (F → F + c).
  - (b) `Fnu hbar lap[h] + Fnu hbar lap[lap[h]]` gives F(p²)(−p² + p⁴): n₀ is undefined.
  - (c) `nk Ak Fnu`, `gk Fnu`: D6's "exactly one reference" says nothing about further couplings, and splitting the coefficient between them is undetermined.
  - (d) `Fnu^2`, `1/Fnu`: D6 does not cover powers of a function.
  - D6 is asserted only in Dimensions (plan:20), which the dimensionful chain Symmetries → … → Equations does not require. These ansätze therefore reach the Equations stage silently.
- Consequence: underdetermined or ill-defined flows.
- Hint: one named check in Ansatz/Symmetries: each functional slot is a single term, number × (at most one non-function coupling) × F¹ × (p²)^{n₀}.
- Resolution: _(author)_

### E5 [major] How many slots a function may occupy is undefined: two slots (duplicates Related) and zero slots (silent) are both possible
- Where: plan:16-17, 24.
- Evidence:
  - `nk Fnu hbar lap[h] - Dk Fnu hbar^2` puts one function in two slots. This gives two equations for one function, and it is the same physics that `Related` would express.
  - Zero slots: `Fnu div[hbar grad[h]]` and `Fnu hbar dt[hbar]` vanish in Γ2. Probe (r1-probe-ansatz.wl) on today's code: `nuk div[hbar grad[h]]` is accepted, its coefficients are {−nuk, −nuk, 1, −Dk}, and Γ2 = {{0, I w}, {−I w, −2 Dk}}, with no nuk.
- Consequence: overdetermined equations, or a function that never flows, with no named check.
- Hint: `Symmetries.functionSlot`: every declared function sits in exactly one non-vanishing slot.
- Resolution: _(author)_

### E6 [major] The deferred decision "functions in Couplings too or separately" decides whether several checks fire
- Where: plan:40-42; plan:47 ("Regulator.knownSymbols already aborts").
- Evidence: Regulator.wl:69-71 allows `Join[cps, {k, q}]`. If functions are in "Couplings", `{hbar,h} -> nk Fnu q^2 rn` passes. The same holds for `NonRenormalised[Fnu]` (Symmetries.wl:66) and for a function given as a reference, as `"Dimensionless"` or as `"Equation"` in Dimensions. If functions are kept separate, the `assertFreshSymbols` and couplingsNotFields checks (Ansatz.wl:150-152) are missing for them.
- Hint: decide now, and list the named checks this requires.
- Resolution: _(author)_

### E7 [major] The D10 tag semantics are underspecified or self-contradictory
- Where: plan:28-32.
- Evidence:
  - Protected or Neglected on a "functional slot" contradicts D8 (such a slot holds a function, which is a coupling; Symmetries.wl:79-80 aborts that case).
  - `Related[FD -> Fnu]`: "flows must agree" has no stated normalisation. Closure is dimensionful, and ∂s(Dk FD) ≠ ∂s(nk Fnu); the identity, if any, is ∂s(Dk FD)/Dk = ∂s(nk Fnu)/nk.
  - The comparable identity in the repo, KS time reversal (Checks.wl:18-19, 31-33), holds only at ν = αD with r_ν = α r_D. The tag has no slot for these conditions.
  - With `"Dimension" -> 1`, c must be ±1. ExactlyZeroQ at a random rational c tests non-physical points, so a true d = 1 identity aborts.
  - The fixed-value tag has no syntax, and its interaction with Normalise is undefined.
- Hint: define each tag on functional slots precisely, including the comparison normalisation, the conditions it holds under, and c² = 1 at d = 1.
- Resolution: _(author)_

### E8 [minor] An applied function symbol (`Fnu[0]`, `Fnu[1]`) is accepted silently
- Where: plan:16; Ansatz.wl:158 (knownSymbols checks symbols and heads only).
- Evidence (r1-probe-ansatz.wl): today `nuk[1] lap[h]` is ACCEPTED with coefficient −nuk[1]. `D[f[0], f]` = 0, so `couplingPowers` (Dimensions.wl) treats `Fnu[0]` as dimensionless.
- Hint: `Ansatz.functionBare`: function symbols appear only bare.

### E9 [minor] Nonsensical uses of `"Normalise"` and of references have no named check
- Where: plan:21-23.
- Evidence: the plan names no abort for any of these: Normalise to a function that does not multiply that reference; Normalise to a dimensionful function; Normalise together with `"Mode" -> "Fixed"`; a function declared as a reference (D7 says "references stay constant"); Normalise plus a fixed-value tag on the same function.
- Hint: add named checks in Dimensions/Adimensional.

### E10 [minor] The D9 measure is loosely declared
- Where: plan:26-27, 68.
- Evidence: K_d is undefined; the identity holds only with K_d = S_{d−1}/(2π)^d. At Dim → 1 the symbolic form is S₀·∫(1−c²)^{−1} = 0·∞, so the d = 1 sum over c = ±1 must be stored as a separate case. A mixed ansatz carries a 1D measure (constant rows) and a 2D measure (functional rows), and the Normalise η row mixes them. The Adimensional contract (MODULES:693-696) is not updated.
- Hint: define the measure per equation and per η row.

### E11 [minor] The cascade omits existing checks that the new atoms break, and one check is circular
- Where: plan:44-56, 67.
- Evidence:
  - `Vertices.gamma2Polynomial` (Vertices.wl:135) and Propagator/Symmetries `Coefficient` assume polynomial entries; `PolynomialQ[x f[x], x]` = False.
  - `Adimensional.linearNoMomentum` requires BetaL free of aq, but BetaL now holds 2ŷ∂ŷF̂ in ap.
  - `Symmetries.truncationCoversAnsatz` (Symmetries.wl:93-96) demands Momentum ≥ n₀ for functional slots, although "truncation applies to function-free slots only".
  - Contradictory `"Pipeline"` values are unspecified.
  - The "independent" planeWave route puts F(q²) in by hand on both sides. Today's planeWave also uses independent leg momenta (Vertices.wl:64-67), where F acting on hbar and F acting on h differ.
- Hint: list each adaptation in the cascade.

### E12 [minor] D10's "WARN with size at p ≠ 0" for numeric slots contradicts D3/D8
- Where: plan:31-32 vs plan:13, 25.
- Evidence: D3/D8 say function-free slots keep the p → 0 classification, yet D10 computes untagged numeric function-free slots at p ≠ 0. Non-numeric constant slots get no such report.
- Hint: state which slots run at generic p.

### E13 [minor] The D12 reference comparison is not yet guaranteed to come from the derivation
- Where: plan:34-36.
- Evidence: "Pointwise in (q,p,c)" holds only under the same loop routing. Integrands that differ by q ↔ −p−q have equal integrals but differ pointwise. Also:
  - no equation numbers are given;
  - I believe there is an erratum to PRE 84, 061128, but I could not check it offline;
  - the plan does not confirm that the paper's LO f's are independent of frequency (D2).
- Hint: compare the symmetrised integrands on both sides, as a check, not as a rewrite of the reference. Cite the arXiv version, the equations and the errata.

### E14 [nit] Notation mixes dimensionful and dimensionless variables
- Where: plan:54 ("scal[p,p] → ap²… (dimensionful p², p·q = p q c)"); ŷ is undefined (y is reserved).

## Nonsensical / ambiguous ansatz cases
| case | named check in the design? |
|---|---|
| function in a cubic or quartic monomial | yes: Ansatz.functionQuadraticOnly |
| function under an operator (`lap[Fnu h]`) | yes: Ansatz.functionNotInOperator |
| undeclared function | yes: Ansatz.functionDeclared / knownSymbols |
| two functions in one slot | yes: Symmetries.functionalSlot |
| function in the regulator | only if functions are kept out of Couplings (E6) |
| quadratic `jac` with a function | yes (existing Ansatz.jacQuadratic) |
| `(nk + Dk) Fnu` | caught late by Dimensions.consistent (existing) |
| constant + function at the same p^{2n₀} | ambiguous (E4a) |
| one function at two powers in a slot | none: n₀ undefined (E4b) |
| function × extra coupling (`nk Ak Fnu`, `gk Fnu`) | late or none (E4c) |
| `Fnu^2`, `1/Fnu` | none (E4d) |
| one function in two slots | none (E5) |
| function whose Γ2 vanishes (`Fnu div[hbar grad h]`, `Fnu hbar dt[hbar]`) | none: silent (E5, probe) |
| applied symbol `Fnu[0]` | none: silent (E8, probe) |
| functional slot whose lower p-orders do not vanish | none: silent 1/p² (E2) |
| Protected/Neglected/Related on a functional slot | contradictory or undefined (E7) |
| NonRenormalised[Fnu], function as reference/Dimensionless | none (E6, E9) |
| Normalise misuse (wrong reference, dimensionful, Fixed, double) | none (E9) |
| `"Pipeline"` contradicting the function list | unspecified (E11) |
| function absent from the density | unspecified (WARN as for couplings?) |
| writing order / factor placement in a quadratic monomial | not ambiguous (held) |

## Appendix: probes (run with `math -script`)
- r1-probe-pointwise.wl: load the package, `quiet[x_] := Block[{Print}, x]`, Get KPZSetup.wl. Take `li = LoopIntegrand[kpzLoops,{hbar,h},0,pp]`, apply DsRRules, DeltaRules and the atoms {GammaAtom[h,hbar,0] → s/(1+s) + s e^{−s}, GammaAtom[hbar,hbar,0] → −2/(1+2s) − 2e^{−s}}, then scal[a,b] → a b, rn, rd → e^{−y}, couplings → 1, dsCoupling → 0. NIntegrate over ω at q = 7/10 for p = 10^{−1…−4}, raw and symmetrised (q → −p−q).
- r1-probe-twofamily.wl: residues of 1/((a1+Iw)²(a1−Iw)(a2+Iw)(a2−Iw)) at I a1 and I a2, Factor, the sum, the limit, and Solve[ap²+aq²+2 ap aq c == aq², c].
- r1-probe-ansatz.wl: DefineAnsatz with `nuk[1] lap[h]`, and with `nuk div[hbar grad[h]]` followed by Gamma2Matrix.
- r1-probe-linear.wl: S_{d−1}∫(1−c²)^{(d−3)/2}dc / S_d = 1, and ∂s F̂|_ŷ = (∂s F|_p)/k^x − x F̂ + 2ŷF̂′ (both 0/1 as expected).

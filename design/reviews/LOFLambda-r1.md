# the-enemy — design/LOFLAMBDA.md, round 1
Attacked commit: 6df4a91 (record untracked; working tree has uncommitted WIP in Kernel/, not attacked)  ·  Date: 2026-11-03
Verdict: 1 blocker, 2 major, 5 minor, 1 nit

## Attack surface
- Confirm green: the target is a decision record. It has no machine checks yet, so there is nothing to run. Code anchors were read at 6df4a91 (`git archive`). Probes were run on that tree.
- §1 l.9–11: Γ⁽¹¹⁾ = iω f_λ + ν p² f_ν; Γ⁽⁰²⁾ = −2D f_D; cubic vertex λ f_λ(k_h̄²) k₁·k₂. Checked by hand: −λ/2 (∇h)², with ∇ → −ik, gives +λ/2 k₁·k₂; the 2 symmetrisations give λ k₁·k₂. — held (ν overloaded, see E8)
- §1 Galilean tie (f_λ on the whole covariant bracket) — re-derived by hand. Under the Galilean shift, D_t h picks up a constant, and f_λ(0)·const is harmless. — held
- §1 l.13 / F3: ∂ν before the ω integral; ν enters only G(Q). Checked against Loops.wl:26-31 (ν only in `gm` and in Γ3, and the KPZ Γ3 has no ω). — held
- F3: `deltaCompletion` (Projection.wl:18-19) handles Δ(−w−ν) through the chain rule. Pole order a_j ≥ 2 is supported (FrequencyIntegration.wl:83-86). — held; the anchor is wrong, see E5
- F3 sign/normalisation of the iν coefficient, and Risks l.46-47 ("acceptance 1 the f_λ → 1 limit"). — finding E1
- F2: "the Protected iω tag at p = 0 proves I_λ(0) = 0". Probe E1: true, but trivially (see E1). — held as a statement; the check is vacuous for F3
- F4: I_ν = f_λ(p) J(p) with J = O(p²), so only f_λ(0) survives. Re-derived through the Laplacian product rule: f′(0) multiplies J(0) ≡ 0. This matches BMWflow model_kpz.py:94,101 (`f_lambda_p**2`, `f_lambda_p`). — held
- F4/L5: η rows need only f_λ′(q). The external h vertex is O(p), so G(Q) and f_λ(Q²) are needed only to first order. This matches eta_calc `f_lambda_derq` (model_kpz.py:193,204). — held
- F5: linear = p f′ (no η) for reference-1 functions. ∂s|p̂ F(κp̂) = ∂s|p F + p∂pF; this agrees with §3.5 "2 ap² f′". — held
- F5 `fp` argument; ETAS `None`; n_η × n_η system. — held
- F7 acceptance 1 (f_λ ≡ 1 versus KPZLO). — finding E2 (and E1)
- F8 acceptance 2: d = 1, 2, 3, z = 2 − η_ν, χ = (2 − d + η_D − η_ν)/2. The formulas are standard. — finding E3 (d = 1)
- F9 Oracle central difference. — finding E6
- F1 `Ansatz.oneFunctionPerEntry` wording. — finding E7
- F1/§3.2 Vertices.wl:113 `loGamma2Check` passes a reference of 1. Read: `c/((f/.funcs) f[ss])` = ±I, numeric. — held
- Vertices.wl:28 "f at the named leg's momentum" (`factorSymbol` with `fnf[scal[qv[l],qv[l]]]`). — held
- Ansatz.wl:167-188, :257-275 anchors. — held
- Causal-pole assumption with c₁ = f_λ(q). — finding E4
- LOPIPELINE L7 / `ExternalMomentum.wardNeglectedLO` (ExternalMomentum.wl:108-109). — finding E8
- Undefined I_λ, J; overloaded ν. — finding E8
- §1 "Γ⁽¹¹⁾ = iω f_λ" without saying which leg carries ω. — finding E9

## Findings
### E1 [blocker] No planned check can catch a sign or normalisation error in the f_λ flow; Risks claims one does
- Where: design/LOFLAMBDA.md:46-47 (Risks: "the Oracle (F9) checks the loop side, acceptance 1 the f_λ → 1 limit"), :24 (F3), :23 (F2), :28 (F7).
- Evidence: Mutate F3 to ∂s f_λ = −(∂ν loop)/(I·norm), or take `norm` with the wrong sign. Every planned check survives:
  - `ExternalMomentum.ward` (F2) is vacuous. The whole integrand of ∂sΓ_{h,h̄}(ν, p) is identically zero at p = 0, for every ν, ω and q. This is because the vertex carrying the external h leg is λ f(k_h̄²) p·k. Running `math -script design/reviews/LOFLambda-r1-E1.wl <6df4a91 tree>` prints `{h, hbar}: ... at p = 0 identically zero? True   at generic p zero? False` (same for `{hbar, h}`). The f_λ vertex has the same p·k structure.
  - F7 (acceptance 1) compares only the f_D and f_ν rows, the η system and β at f_λ ≡ 1. The f_λ row is never compared, and the I_D and I_ν rows do not contain I_λ.
  - F9 differentiates the same `LoopIntegrand` numerically. It checks the ω integration only. If it compares after the division by I·norm, it shares the normalisation with the code under test.
- Consequence: the one new flow equation of the record, the one that differs from KPZLO, has no independent check. A sign error would let f_λ run the wrong way, and that changes the exponents of acceptance 2. The user would only see it "visually", with nothing to compare against.
- Hint: add an independent reference for the f_λ row, either a hand-derived or KCW 2012 kernel at generic f_λ or a closed-form limit, and state that F9 compares before normalisation.
- Resolution: fixed: design/LOFLAMBDA.md F10 — DE route (benchmarks/KPZLOLambda, KPZLOLambda.de.flam; sign flip rejected).

### E2 [major] Acceptance 1 tests only f_λ ≡ 1, although a generic-f_λ reference for the f_D, f_ν and η rows is in reach
- Where: design/LOFLAMBDA.md:16, :28 (F7).
- Evidence: Setting f_λ → 1 and f_λ′ → 0 erases where f_λ sits: at p, q or Q, and with what power. Yet BMWflow `flow/models/model_kpz.py:78-101` (I_D, I_ν with `f_lambdaq`, `f_lambdaQ`, `f_lambda_p`, "kloss2012_omega_integration.nb") and `:193-204` (η with `f`, `qdf`) keep the Kloss-2012 f_λ placeholders. A spot check shows they are usable. `math -script design/reviews/LOFLambda-r1-E3.wl` prints `Int dw/2Pi C C / (4 kq kQ) - A3a = 0`: their A3a is the correct ω integral with Γ⁽¹¹⁾ = iω f_λ + l.
- Consequence: a misplaced f_λ, for example f_λ(q) where f_λ(Q) belongs, passes F7.
- Hint: add an F7 variant that compares with ModelKPZ's LO formulas at symbolic f_λ(p), f_λ(q), f_λ(Q), f_λ′(q). Those placeholders were never exercised with f_λ ≠ 1, so check them too.
- Resolution: fixed: F10 — generic-f_λ reference from model_kpz.py (KPZLOLambda.reference.D, .nu).

### E3 [major] In d = 1, f_λ ≠ 1 breaks the KPZ time-reversal symmetry; the record does not discuss it, but acceptance 2 runs d = 1
- Where: design/LOFLAMBDA.md:16-17, :29 (F8: d = 1, 2, 3).
- Evidence: In d = 1 the symmetry h̄ → h̄ + (ν/D)∂²h (t → −t) leaves the λ term invariant only because ∫ h″ (h′)² = 0. With the ansatz the term becomes ∫ h″ f_λ(−∂²)(h′)². The ∂t term stays a total derivative, because f_λ is a symmetric operator. `math -script design/reviews/LOFLambda-r1-E2.wl` gives `flam = 1: 0` and `flam = 1 + a k^2: (26*a*Pi)/5` (= −a ∫(h″)³), which is nonzero.
- Consequence: the exact d = 1 values χ = 1/2, z = 3/2, which the ansatz keeps at f_λ ≡ 1, need not survive. A d = 1 deviation in the visual comparison could be blamed on the code or on the physics. My recollection, not checked, is that KCW 2012 fix f_λ ≡ 1 in d = 1.
- Hint: decide and record what is done in d = 1 (f_λ ≡ 1, or an expected deviation), and add the exact d = 1 values as a criterion.
- Resolution: fixed: §4 of the record (d = 1 caveat), README and runBMWflow.py docstring.

### E4 [minor] f_λ(q) > 0 becomes an undeclared residue assumption on a flowing function
- Where: design/LOFLAMBDA.md:24 (F3), :26 (F5, handoff). The code is Propagator.wl:91 (deg 1: c₀/c₁ > 0) and ExternalMomentum.wl:142-145 (`causalPoles`, not exported, LOPIPELINE L8).
- Evidence: with c₁ = f_λ(q²), the residue choice needs f_λ(q), f_λ(Q) > 0 at every q. In KPZLO, c₁ = 1.
- Consequence: if f_λ crosses zero during the BMWflow run, the kernels are silently wrong. Nothing on either side detects it.
- Hint: state f_λ > 0 in the record and in the v3 handoff, with a cheap runtime check suggested for BMWflow.
- Resolution: fixed: §4 (causal poles).

### E5 [minor] Wrong code anchor for the pole order a_j ≥ 2
- Where: design/LOFLAMBDA.md:24 ("FrequencyIntegration.wl:120-134").
- Evidence: `git show 6df4a91:Kernel/FrequencyIntegration.wl | sed -n 120,134p` shows `samplePoint`/`sampleValue`. The residues at order a_j are at :83-86 (`FrequencyIntegralFamilies` :75-91), and the single family is at :19-24.
- Consequence: misleads the implementer and later reviewers.
- Hint: change the anchor to :75-91.
- Resolution: fixed: anchor FrequencyIntegration.wl:75-91.

### E6 [minor] F9 states no tolerance and does not say which side of the normalisation it compares
- Where: design/LOFLAMBDA.md:30.
- Evidence: a central difference with Δν = 10⁻⁶ has a truncation error of O(Δν²·∂ν³) ≈ 10⁻¹². The record does not say whether the comparison is relative or absolute, or whether it is made before or after /(I·norm) (see E1).
- Hint: give the tolerance, and compare the raw ∂ν loop.
- Resolution: fixed: F9 (exact ∂ν, before normalisation, 1e-8).

### E7 [minor] `Ansatz.oneFunctionPerEntry` as worded would reject the target ansatz
- Where: design/LOFLAMBDA.md:22 (F1): "a Γ2 entry carries at most one function over its ω orders".
- Evidence: the KPZ Γ⁽¹¹⁾ entry carries f_λ (ω¹) and f_ν (ω⁰). §3.3 ("per ω order") implies the intended meaning is "at most one per ω order".
- Hint: write "at most one function per ω order of an entry".
- Resolution: fixed: Ansatz.oneFunctionPerCoefficient.

### E8 [minor] L7 is partly superseded without saying so; loose symbols
- Where: design/LOFLAMBDA.md:3 (supersedes only L1, L2); LOPIPELINE.md L7 ("Protected iω … Neglected (LO) at p ≠ 0"); ExternalMomentum.wl:108-109 (INFO card `wardNeglectedLO`); design/LOFLAMBDA.md:9-13, :23, :25.
- Evidence: with f_λ, the iω coefficient at p ≠ 0 is the f_λ flow, not neglected, so the INFO card becomes false. §3.3 does not plan to update it. Other loose symbols: ν means the viscosity (l.9, 11) and also the external frequency (l.13, F3); I_λ and J are never defined.
- Hint: add L7 to "supersedes", update the card, and rename the external frequency (or the viscosity) in the record.
- Resolution: fixed: notation line; wardNeglectedLO only without a function.

### E9 [nit] "Γ⁽¹¹⁾ = iω f_λ" does not say which leg carries ω
- Where: design/LOFLAMBDA.md:11.
- Evidence: in the VertexValue convention, Γ_{h,h̄}(ω,p) = +iω f_λ and Γ_{h̄,h}(ω,p) = −iω f_λ.
- Hint: write Γ_{h h̄}(ω, p) with h at (ω, p).
- Resolution: fixed: §1.

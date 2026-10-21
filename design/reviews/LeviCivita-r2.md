# the-enemy — /home/greg/.claude/plans/i-want-to-use-greedy-clover.md, round 2
Attacked commit: 1fbb12e  ·  Date: 2026-10-21
Verdict: 0 blocker, 0 major, 4 minor, 2 nit

Note: plan mode is active, so this report is written here and not at `reviews/i-want-to-use-greedy-clover/r2.md`. The probes are
inline `math -noprompt -run` commands (reproduced in the findings). The target is a design plan with no machine checks of its own, so
"confirm green" does not apply. Round 1 is in `i-want-to-use-greedy-clover-agent-ac0522a3db1dae127.md`.

## Attack surface
- R1-E1 resolution (`Ansatz.singleEpsilon`, distinct `EpsIndex[j,_]` labels) — emulated plan canonical on `jac[..]^2` (Power shares one label, Ansatz.wl:108) — finding E1
- R1-E2 resolution (κ Hbar grad[Psi].grad[H]) — hand: G block-diagonal, Γ3_{Hbar Psi H} = λ·wedge + κ·scal, bubble cross term λκ has e = 1; κ² gives nonzero n = 0 flows — held
- R1-E3 resolution (declared `"Dimension" -> 2`) — Theory.wl:28 today vs plan rule, probe — finding E2; grep of every `Dimension` use (Theory, Ansatz:122/129, Regulator:44, Symmetries:105, Projection:45, Oracle:26, Assemble:128/133) — held, except E6
- R1-E4, E6(a)(b)(c), E7, E8 resolutions — held (E6(c) wording: E5)
- F2 + `dimFactors` with dim = 2 vs explicit 2D Laplacian — library `FRGscalarMSR`Private`project` on (p²q² − (p·q)²) f: gives q² f(0,q²,0)/2, equal to explicit-component Laplacian/(2·2); with symbolic d the result differs ((d−1) factor), confirming d = 2 must be fixed — held
- F4 pointwise in q (Oracle at fixed qVec valid) — held (round 1 P1 rerun by hand)
- Every 3-vertex wedge ∝ wedge[p,q] (momentum conservation) ⇒ λ never enters n = 0 flows, λ² only n = 1 — held (consistent with the E2 resolution's "n = 0 flow nonzero" coming from κ)
- canonical completeness for one ε (only symmetry = slot swap; two ε indices never land on one factor) — held
- Quadratic jac monomials — finding E3
- Acceptance "shallow water with Dimension -> 2 = existing with Dim -> 2" — finding E4
- Oracle `"Dimensions" -> Automatic` vs byte-identical card — finding E5
- Literal-2 hazard list and stale body text — finding E6
- Vertices/planeWave 2D, unit hand value −c wedge[q_H,q_Psi] — held (round 1)

## Findings
### E1 [minor] (re R1-E1) `singleEpsilon` counts distinct labels, so `jac[a,b]^2` passes it; the "canonical is now safe" claim rests only on momentumOrder
- Where: plan §"Resolutions" E1 ("counted by distinct EpsIndex[j,_] labels … This makes the anonymous EpsIndex[s] in canonical safe"); Ansatz.wl:108 (`S[Power[x_, n_]] := Expand[S[x]^n]` squares one expansion, so both ε share one j).
- Evidence: emulating the plan's toFactors + canonical on `Hbar jac[H+Psi, H-Psi]^2` (true value 4 jac[H,Psi]² ≠ 0): "distinct eps labels per raw term (max): 1", "surviving monomials: 3 of 6": three groups are declared zero by the collapsed swap. A search over 144 pairs (a,b) of simple field polynomials found no case where all monomials cancel, so momentumOrder (Ansatz.wl:207) still aborts. The abort shows the wrong message, and safety is empirical, not proven.
- Consequence: an unproved safety claim. If some input cancels fully, an unsupported p⁴ term is dropped with only an INFO line, which is exactly the R1-E1 failure.
- Hint: count ε slots (`Count[…, EpsIndex[_, 1]] > 1`), not distinct labels, and add `Hbar jac[H,Psi]^2` to the unit test.
- Resolution: _

### E2 [minor] Accepting a literal 2 lets a stale `Dim = 2` silently become a declared d = 2
- Where: plan §Theory ("accepts … `"Dimension" -> 2` … So d is never fixed silently"); Theory.wl:28-29.
- Evidence: `Get[".../FRGscalarMSR.wl"]; ToExpression["Dim = 2; CheckAbort[DefineTheory[\"Scalars\" -> {H}, \"Dimension\" -> Dim], \"ABORTED\"]"]` today prints `FAIL check: Theory.dimension … offending: 2` and aborts. Under the plan, the evaluated option value 2 passes. The ε-free pipeline then runs in d = 2, and the Oracle runs in d = 2 only.
- Consequence: a user who meant d to stay symbolic gets d = 2 equations. The only signal is "dimension 2" on the cards. This removes the `.values` protection that TheoryTest.wl:18-21 enforces for the other symbols.
- Hint: take the fixed dimension through a separate option (e.g. `"FixedDimension" -> 2`), or give DefineTheory a Hold attribute and inspect the unevaluated value.
- Resolution: _

### E3 [minor] A quadratic `jac` monomial survives the parser but is identically zero in Fourier space, so its coupling silently does nothing
- Where: plan F3 ("Γ2 is unaffected, because wedge[q, −q] = 0"); §Ansatz canonical (L ≠ L′ for distinct fields); Ansatz.wl:213 (unusedCouplings looks only at monomial coefficients).
- Evidence: for `c jac[Hbar, H]`, canonical gives {Hbar{1}, H{2}} vs {Hbar{2}, H{1}}, which differ, so the monomial is kept with coefficient c. Its vertex is (−i)² c wedge[q, −q] = 0 (the term is the total derivative ∂_i(ε_ij Hbar ∂_jH)). The plane-wave check also gives 0, so Vertices.planeWave passes, and Symmetries finds no Γ2 coefficient containing c. No check fires, and the ansatz table shows a term that has no effect.
- Consequence: a user's declared quadratic coupling gets no equation and no warning.
- Hint: in Ansatz, abort or WARN on any ε monomial with FieldCount = 2 (a total derivative in d = 2).
- Resolution: _

### E4 [minor] The acceptance criterion "exactly the existing equations with Dim -> 2" is undefined against Assemble's form selection
- Where: plan §Resolutions "Declared-dimension tests" (shallow water); Assemble.wl:55-65 (smallest-LeafCount candidate, `TimeConstrained[Simplify[#], 10, #]`).
- Evidence: Assemble picks the output form by LeafCount and by a time-limited Simplify. An input with the number 2 in place of `Dim` therefore need not give a structurally identical (`===`) integrand to `(old /. Dim -> 2)`.
- Consequence: read as `===`, the test fails or is flaky. Read loosely, it is untestable as written.
- Hint: state the comparison as `ExactlyZeroQ[new − (old /. Dim -> 2)]` per equation, with the same keys and classes.
- Resolution: _

### E5 [nit] `"Dimensions" -> Automatic` changes the Oracle card text on the ε-free path, and its abort rule is ambiguous for a symbolic d
- Where: plan E6(c); Oracle.wl:14 (usage), :37, :96 (`"d in " <> ToString[OptionValue["Dimensions"]]`).
- Evidence: if Automatic is not resolved before :96, the INFO line reads "d in Automatic" instead of "d in {2, 3}". This breaks "same card text". "An explicit list containing anything other than the declared value aborts" read literally rejects `{2,3}` when the declared value is the symbol Dim.
- Hint: resolve Automatic first and print the resolved list, and restrict `Oracle.declaredDimension` to integer declarations.

### E6 [nit] The literal-2 hazard list is misstated, and the plan body still contradicts its resolutions
- Where: plan §Theory "Literal-2 hazards"; §Ansatz Checks ("caught by momentumOrder anyway"); §Oracle (`Oracle.epsilonDimension`, abort only on 3).
- Evidence: with only d = 2 allowed, the Oracle rule `eqs["Dimension"] -> d` is `2 -> 2`, a no-op. Ansatz.wl:122/129 and Regulator.wl:44 take Intersection/Complement with the integer 2, which is harmless. The guards cost nothing, but the stated reason is wrong. The body text still states superseded designs.
- Hint: fold the resolutions into the body sections before implementing.

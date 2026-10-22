# Future seams

These are extensions the v1 design keeps room for. None of them is implemented yet.

- **Pipeline 2: couplings that depend on ω and p.**
  - Selected with the Ansatz `"Pipeline"` option. It replaces `ProjectOnMonomials`, and the external ν, p stay generic.
  - The input contract of `IntegrateFrequencies` is "rational in ω with listed pole factors", so shifted poles Δ(ω+ν, p+q) plug
    in without changes.
- **Δ of degree greater than 2 in ω** (N ≥ 3 scalars, or richer frequency structure).
  - Needs a general reduction by symmetric functions of the roots.
  - The validity condition becomes the Routh–Hurwitz criterion: an algebraic test on the coefficients of Δ that tells whether
    all its roots lie in one half plane, without computing the roots.
- **Causal poles checked in FlowRG.** `causal_poles(couplings…, aq)` is exported (design/CAUSALITY.md) but FlowRG does not call
  it. FlowRG could evaluate it on its aq nodes and at the ends of the domain (aq → 0, where the shallow-water condition reduces to
  agk aHk alk > 0) and stop the flow or flag a fixed point on a violation; sampling is not a proof between nodes. Through a
  `docs/handoff/` spec.
- **p⁴ in FlowRG.** The symbolic chain projects any (p²)ⁿ (design/MOMENTUMPROJECTION.md), but `JuliaExport` exports only r, r′,
  r″: p⁴ flows need r‴, r⁗ bound by FlowRG (an interface extension through a `docs/handoff/` spec), and a smooth (C⁴) shape.
- **More than one `jac` per monomial.** The ε reduction (F2) handles one pair of ε indices; two would need the product rule for
  two ε pairs and the matching canonical form.
- **Generic-d p⁴ acceptance.** `benchmarks/KuramotoSivashinsky` checks the p⁴ flows end to end against the literature
  (Gosteva, Wschebor, Canet 2026), but only at d = 1, where the n = 2 normalisation 8d(d+2) equals 8(d²+2). A literature p⁴ flow
  at symbolic d would close the loop. Candidates: the Lai–Das Sarma equation near d = 4 − ε, whose O(ε) result z = 4 − ε/3 does
  not depend on the regulator; the KPZ derivative expansion of Sánchez-Villalobos, Canet, Wschebor (in preparation); or the
  shallow-water DE4plus flows.
- **Independent check of the Γ⁽⁴⁾ tadpole normalisation (E4).**
  - The bubble is checked against Canet's KPZ flow (`tests/acceptance/KPZ.wl`). The tadpole needs a theory with a quartic
    vertex.
  - Options: the background-field form ½ ∂̃s Tr ln(Γ2[φ₀] + R) on a toy quartic theory, or a published one-loop flow of a
    model with Γ⁽⁴⁾ ≠ 0.
- **∂sΓ⁽³⁾ projections.**
  - These need Γ5 and the triangle and box diagrams, which means a new variant of Loops.
  - With them, `NonRenormalised` becomes checkable.
- **3D Levi-Civita (curl-type) operators.**
  - `jac` covers ε_ij in d = 2 only (design/LEVICIVITA.md). ε_ijk ∂_j A_k needs vector-valued operands and the d = 3 identity
    ε_ijk ε_lmn = det[δ], plus its own declared dimension and oracle.
- **Ward identities.** Derive and check the symmetry constraints symbolically. Today the user declares them.
- **Non-reference couplings in the regulator** (Adimensional, plan D4). Today a prefactor such as `nk (1 + Ak) q^2 arnh` aborts,
  because ∂sR would contain β_A itself. The exact treatment: one vector of unknowns X = (η_1 … η_{N+1}, ∂s c for every
  regulator coupling c); each β written once and linear in X; the rows ∂s c = β_c(X) join the η system; FlowRG gets a
  `const SOLVED = (:A,)` declaring couplings whose β is the solved unknown.
- **η fixed by a relation** (Adimensional). Today a fixed η is a number given at run time (FlowRG `fixed_etas`). Extensions:
  a linear relation among the η's (`"Fixed" -> <|Dphik -> 2 etanu - etah + 1|>`), or a quantity kept scale invariant
  (`"ScaleInvariant" -> Dphik^a nk^b k^c`), from which the library derives the relation through `Dimensions`.
- **Second shape derivatives in FlowRG.** The export writes r″ as `D2<shape>(aq^2)` (DE4plus g2, D2φ); FlowRG must bind them
  (`docs/handoff/FlowRG-MUST_FIX.md`, task 2). For Litim, r″ contains δ(y − 1), which a pointwise `D2` cannot represent.

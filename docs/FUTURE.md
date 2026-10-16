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
- **Projections at p⁴ and higher (E2).**
  - Write p = |p| e, Taylor-expand in |p|, and average over directions with
    ⟨(e·q)^{2m}⟩ = q^{2m} (2m−1)!! / [d(d+2)⋯(d+2m−2)].
  - This replaces the Laplacian trick at every order. Validate it with the oracle in d = 2 and d = 3.
- **Independent check of the loop normalisation (E3/E4; declined for v1).**
  - Use the background-field form ½ ∂̃s Tr ln(Γ2[φ₀] + R) on KPZ and on a toy quartic theory.
- **∂sΓ⁽³⁾ projections.**
  - These need Γ5 and the triangle and box diagrams, which means a new variant of Loops.
  - With them, `NonRenormalised` becomes checkable.
- **Ward identities.** Derive and check the symmetry constraints symbolically. Today the user declares them.
- **Adimensional.wl and JuliaExport.wl.**
  - Map the `Equations` stage onto FlowRG's `docs/equations-interface.md`: COUPLINGS/ETAS/SHAPES, positional arguments, β and an
    η-system linear in η, the l/nl split, and shapes as `armh(aq^2)`, `Darmh(aq^2)`.
  - FlowRG supports only first derivatives of the shapes today. g2 and D2phi need r″, so the interface must be extended.

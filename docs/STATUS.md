# Status

What the library does today, how it is checked, its known limits, and the planned extensions. Each design record in
[design/](../design/) holds the decision and the reasons behind a feature.

## Implemented

| Feature | Since | Design record |
|---|---|---|
| Derivative-expansion pipeline, N = 1 or 2 scalars: ansatz → vertices → propagator → symmetry classes → loops → projection → residues → closure → equations | v1 | [FRGDESIGN.md](../design/FRGDESIGN.md) |
| 2D Levi-Civita term `jac[A, B]` with `FixedDimension[2]` | v1 | [LEVICIVITA.md](../design/LEVICIVITA.md) |
| Momentum projection on any (p²)ⁿ | after v1 | [MOMENTUMPROJECTION.md](../design/MOMENTUMPROJECTION.md) |
| Dimensions derived from the ansatz, dimensionless β functions and η system, FlowRG (Julia) export | after v1 | [ADIMENSIONAL.md](../design/ADIMENSIONAL.md) |
| Causal-pole conditions, minimal and dimensionless, exported as `causal_poles` | after v1 | [CAUSALITY.md](../design/CAUSALITY.md) |
| LO pipeline: couplings that are functions of p², flows at generic p, BMWflow (Python) export | after v1 | [LOPIPELINE.md](../design/LOPIPELINE.md) |
| LO functions at iω and in vertices (KPZ's f_λ), ∂/∂(iν) flows, BMWflow contract v3 | latest | [LOFLAMBDA.md](../design/LOFLAMBDA.md) |

Usage: [flow equations](guide/flow-equations.md), [dimensionless flows and export](guide/dimensionless-and-export.md),
[LO pipeline](guide/lo-pipeline.md).

## Checks against independent results

- **Shallow water**: the 9 DE4plus equations reproduce `DE4plus_fast.wl`'s `flowsFAST` exactly, and the numerical oracle confirms
  them in d = 2 and 3.
- **LPA**: the LPA ansatz (no K, g2, D2phi) reproduces the old notebooks' LPA equations exactly (`flowsDE4plus_agreed` reduced to
  LPA, and the original `flowsLPAold`).
- **DE4 vs the old notebooks**: exact agreement with `flowsDE4plus_agreed` (g, Dϕ, N, Dh, A, K) and `flowsNEW` (g, Dϕ, N, Dh, A,
  K, H, g2). The old agreed g2 and D2ϕ disagree and are considered wrong (`tests/reference/legacy/README.md`).
- **KPZ**: the one-loop flow reproduces L. Canet's closed-form flow (arXiv:2509.01472, Eq. 20) exactly, for any d, before and
  after adimensionalisation.
- **2D `jac`**: the equations agree with the numerical oracle in d = 2 (`tests/acceptance/Jacobian2D.wl`).
- **FlowRG export**: generated files are evaluated in Julia and loaded by FlowRG; DE4plus agrees with the old hand-built
  `Equations2D.jl` up to its two known errors.
- **Literature benchmarks** ([benchmarks/](../benchmarks/README.md)): Kuramoto–Sivashinsky DE flows in d = 1 including p⁴ (9/9
  exact), KPZ LO (8/8 exact), KPZ LO with f_λ (14/14 exact).

## Known limitations

- **Ward identities** are assumed (declared by the user as `Protected`), not derived. `NonRenormalised` vertex couplings are
  not checkable, because that needs ∂sΓ⁽³⁾.
- **Loop normalisation (E3/E4).** The trace formula is hand-derived in [MODULES.md](MODULES.md) (Loops).
  - Its bubble part (prefactor, sign, 2π) is checked independently against Canet's KPZ flow (`tests/acceptance/KPZ.wl`).
  - The −½ Γ⁽⁴⁾ tadpole is **untested**, and the cards say so whenever Γ⁽⁴⁾ ≠ 0.
  - Agreement with `flowsFAST` and with the oracle does not test the normalisation, because they share the formula.
- **Truncation.** Δ (the determinant whose roots are the propagator poles) must be of degree ≤ 2 in ω, so at most 2 scalars.
  The symbolic chain projects any (p²)ⁿ, but the FlowRG export binds only r, r′, r″ (p⁴ flows need r‴, r⁗).
- **Validity conditions.** The causal-pole conditions must hold at every loop momentum. The library does not check them: they
  are exported, minimal and dimensionless, as `causal_poles` in the equations file, and FlowRG does not call it yet.
  Stripping their positive factors assumes the reference couplings are > 0.
- **Regulator prefactors** may contain reference couplings only (`nk q^2 arnh`, not `nk (1 + Ak) q^2 arnh`).
- **Second shape derivatives.** Flows of momentum-derivative couplings (DE4plus g2, D2φ) need r″ (`D2<shape>`). FlowRG binds
  it for its smooth regulators (Wetterich, Exponential, ForcingExponential) but not for Litim, whose r″ contains δ(y − 1)
  ([handoff/FlowRG-MUST_FIX.md](handoff/FlowRG-MUST_FIX.md), Status).

## Not implemented and future directions

The v1 design keeps room for these extensions. None of them is implemented yet.

- **Pipeline 2: couplings that depend on ω and p.** The LO pipeline (functions of p² at ν = 0) is implemented. Shifted poles do
  **not** plug into `IntegrateFrequencies` without changes: its dispatch maps every Δ(ω, ·) to the same pole factor (it now
  aborts with `FrequencyIntegration.singleMomentum`); a loop with Δ at q and at Q = p + q needs `FrequencyIntegralFamilies`.
  Δ(ω + ν, ·) with ν ≠ 0 (functions of ω) needs a further family shift.
- **LO pipeline extensions** ([LOPIPELINE.md](../design/LOPIPELINE.md)):
  - p = 0 rows: the exact limit or a small-p series instead of BMWflow's 'const'/'quad' extrapolation (L6).
  - Causal-pole conditions for functions exported (today a residue assumption on the card, L8).
  - Relations between functions (d = 1 time reversal, f_ν = f_D) declared and checked (L7).
  - Function × iω with a running reference (Z_k f_Z(p²) iω: its η from ∂/∂(iν) at p = 0), or × ω² (`Ansatz.omegaFunctionReference`);
    reference-1 functions at iω (KPZ's f_λ) are implemented ([LOFLAMBDA.md](../design/LOFLAMBDA.md)).
  - Causal-pole condition f_λ(q) > 0 for a function at iω checked at run time ([LOFLAMBDA.md](../design/LOFLAMBDA.md) §4).
  - f″ / r″ in η rows (models beyond KPZ-like, `ExternalMomentum.firstDerivativesOnly`): needs r″ in BMWflow or an integration by
    parts that does not bring r″ through ∂sR (L5).
  - Rows divided by p⁴ or more (`PythonExport.p2DividedAtMostOne`): BMWflow's θ quadrature error divided by p₁⁴ breaks its tolerance.
  - More than one exported coupling (`PythonExport.singleCoupling`): BMWflow's state has one scalar g.
  - `FixedDimension`/`jac` in LO (`ExternalMomentum.isotropicAngles`): the angle dependence beyond cos θ.
  - Mixed ansätze, DE monomials next to functions (`Ansatz.mixedNotYet`): the generic-p flow is the seam (DE = its Taylor expansion).
  - Speed of the two-family residues for N = 2: a two-field KPZ-like LO model (Δ of degree 2, four families of powers) does not
    finish in 2 minutes; a partial-fraction or numeric-root route would be needed.
  - Lessons for ExternalMomentum and PythonExport (L13).
- **Δ of degree greater than 2 in ω** (N ≥ 3 scalars, or richer frequency structure).
  - Needs a general reduction by symmetric functions of the roots.
  - The validity condition becomes the Routh–Hurwitz criterion: an algebraic test on the coefficients of Δ that tells whether
    all its roots lie in one half plane, without computing the roots.
- **Causal poles checked in FlowRG.** `causal_poles(couplings…, aq)` is exported ([CAUSALITY.md](../design/CAUSALITY.md)) but
  FlowRG does not call it. FlowRG could evaluate it on its aq nodes and at the ends of the domain (aq → 0, where the
  shallow-water condition reduces to agk aHk alk > 0) and stop the flow or flag a fixed point on a violation; sampling is not a
  proof between nodes. Through a `docs/handoff/` spec.
- **p⁴ in FlowRG.** The symbolic chain projects any (p²)ⁿ ([MOMENTUMPROJECTION.md](../design/MOMENTUMPROJECTION.md)), but
  `JuliaExport` exports only r, r′, r″: p⁴ flows need r‴, r⁗ bound by FlowRG (an interface extension through a `docs/handoff/`
  spec), and a smooth (C⁴) shape.
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
  - `jac` covers ε_ij in d = 2 only ([LEVICIVITA.md](../design/LEVICIVITA.md)). ε_ijk ∂_j A_k needs vector-valued operands and
    the d = 3 identity ε_ijk ε_lmn = det[δ], plus its own declared dimension and oracle.
- **Ward identities.** Derive and check the symmetry constraints symbolically. Today the user declares them.
- **Non-reference couplings in the regulator** (Adimensional, plan D4). Today a prefactor such as `nk (1 + Ak) q^2 arnh` aborts,
  because ∂sR would contain β_A itself. The exact treatment: one vector of unknowns X = (η_1 … η_{N+1}, ∂s c for every
  regulator coupling c); each β written once and linear in X; the rows ∂s c = β_c(X) join the η system; FlowRG gets a
  `const SOLVED = (:A,)` declaring couplings whose β is the solved unknown.
- **η fixed by a relation** (Adimensional). Today a fixed η is a number given at run time (FlowRG `fixed_etas`). Extensions:
  a linear relation among the η's (`"Fixed" -> <|Dphik -> 2 etanu - etah + 1|>`), or a quantity kept scale invariant
  (`"ScaleInvariant" -> Dphik^a nk^b k^c`), from which the library derives the relation through `Dimensions`.

# Conventions

This header is printed above every vertex table (from milestone 2 on). Every module follows it.

## Fields

- There are N physical scalars φ_a, each with an MSR response field φ̄_a. `DefineTheory` creates φ̄_a automatically as the symbol
  `<name>bar`.
- The field order is (φ1, φ̄1, φ2, φ̄2, …). It indexes every matrix (Γ2, R, G).
- Causality is built in: Γ2_{φφ} = 0 for the physical fields.

## Fourier transform

- φ(t,x) = ∫ dω dᵈq/(2π)^{d+1} e^{−i(q·x − ωt)} φ(ω,q). Hence **∂t → +iω** and **∇ → −iq**.
  - This sign choice reproduces Γ_{HH̄}(ω) = iω + f_η of `DE4plus_fast.wl`.
  - The opposite sign gives the time-mirrored theory, with its poles in the wrong half plane.
- Functional derivatives are normalised by δφ(ω,q)/δφ(ω′,q′) = (2π)^{d+1} δ(ω−ω′) δ(q−q′).
- Vertices:
  δⁿΓ/δφ_{a1}(ω1,q1)…δφ_{an}(ωn,qn) = (2π)^{d+1} δ(Σω) δ(Σq) · Γ⁽ⁿ⁾_{a1…an}(ω1,q1; …; ω_{n−1},q_{n−1}).
  - All momenta are incoming.
  - The delta function is stripped.
  - The last leg carries −(the sum of the others).
  - Vertices are fully symmetrised over identical fields.
- Γ2_ab(ω,q) ≡ Γ⁽²⁾_ab with leg a at (ω,q) and leg b at (−ω,−q). Causal (retarded) poles lie in the **upper** half ω-plane.

## Vectors, dimension, scale

- Vectors are bare symbols. Their scalar product is `scal[a, b]`, which is Orderless and linear over sums and over explicit
  **numbers** only. It is not linear over couplings: `scal[nk a, b]` stays as it is.
- The space dimension d is a user symbol (e.g. `Dim`) throughout. The one exception is `"Dimension" -> FixedDimension[2]`,
  which `jac` requires. d is never fixed silently: a literal 2, or a `Dim` that got the value 2, aborts.
- In d = 2, `wedge[a, b]` = ε_ij a_i b_j, with **ε₁₂ = +1**. It is antisymmetric (kept in canonical order, not Orderless), and
  linear over sums and over explicit numbers only, like `scal`.
- ∂s ≡ k∂k, where k is the RG scale.
- Shape functions are functions of y = q²/k².

## Projections (derivative-expansion pipeline)

- The p² coefficient is ∇²_p/(2d) at p = 0. This is valid **only for p²**: monomials of order p⁴ abort. The general method is in
  FUTURE.md.
- The (iν)^m coefficient is (∂/∂(iν))^m / m! at ν = 0.

## Ansatz density (DSL)

- Γ_k = ∫dt dᵈx density. The density is written with the package-owned, protected, inert operators `dt`, `grad`, `div`, `lap`,
  plus `.` between two gradients, and `jac[a, b]` = ε_ij ∂_i a ∂_j b (d = 2 only, a and b scalars, at most one per monomial).
  With ∇ → −iq, `jac[A, B]` with A at leg l and B at leg m gives −`wedge[q_l, q_m]`.
  - Do not use System `Dt`, `D`, `Grad`, `Div` or `Laplacian`. `Dt[nk H]` silently adds `H Dt[nk]`, and the System operators
    applied to explicit coordinates evaluate to 0 before the library sees them.
- Couplings depend on k only. They stand outside the operators: write `nk dt[H]`, not `dt[nk H]`.
- Each parsed monomial is coefficient × Π `FieldFactor[φ, n_t, {indices}]`. `FieldFactor[φ, n_t, {i, …}]` is ∂t^{n_t} ∂_i … φ.
  A repeated `SpaceIndex[1]` is summed over. The two indices of ε are `EpsIndex[1]` (on the factor of a) and `EpsIndex[2]` (on
  the factor of b).

## Reserved names

Fields and couplings may not be named `k q p y t x ω ν`. These names have a fixed meaning in library input and output: RG
scale, loop momentum, external momentum, shape argument, time, position, frequencies. They may not be System symbols or
FRGscalarMSR symbols either, and they must have no definitions.

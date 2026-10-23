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

## Dimensions and dimensionless couplings (Dimensions, Adimensional, JuliaExport)

- Γ_k = ∫dt dᵈx 𝓛 is dimensionless; [x] = 1/k. A dimension is written k^x ref₁^{y₁} ⋯ ref_{N+1}^{y_{N+1}} in the basis of the RG
  scale k and N+1 user-chosen reference couplings (N scalars).
- Anomalous dimensions: **∂s ref_j = −η_j ref_j**. A coupling of dimension k^x Π ref_j^{y_j} has ĉ = c/(k^x Π ref_j^{y_j}) and
  ∂s ĉ = ĉ (−x + Σ_j y_j η_j) + (∂s c)/[c]. The frequency exponent is z = x_ω − Σ_j y_{ω,j} η_j for [ω] = k^{x_ω} Π ref_j^{y_{ω,j}}
  (e.g. [ω] = k² nk gives z = 2 − η_ν).
- Names: ĉ = `a<>c` for a dimensionful c (gk → agk); a dimensionless coupling keeps its name; the η of a reference named "nu" is
  `etanu` (FlowRG name `:nu`).
- Loop momentum: q = k aq. ∫dᵈq/(2π)ᵈ f(q²) = k^d ∫daq LoopMeasure(d, aq) f(k² aq²), LoopMeasure = S_d aq^{d−1}/(2π)^d,
  S_d = 2π^{d/2}/Γ(d/2). The exported loop parts include this measure; the old hand-built exports did not have the (2π)^{−d}.
- Reference couplings are **> 0** (and k, aq > 0). The causal-pole conditions are exported with their powers stripped
  (`causal_poles`, design/CAUSALITY.md). "Causality" alone names Γ2_{φφ} = 0 (the `Causality` symmetry class); "causal poles"
  names the roots of Δ(ω) in the upper half plane.

## Projections (derivative-expansion pipeline)

- The (p²)ⁿ coefficient, any n, is the direction average of the s^{2n} Taylor coefficient along p = s e:
  scal[p,p] → s², scal[p,q] → s c, ⟨c^{2k}⟩ = q^{2k}(2k−1)!!/[d(d+2)⋯(d+2k−2)]. It equals (∇²_p)ⁿ at p = 0 / N_n,
  N_n = 4ⁿ n! (d/2)_n (N_1 = 2d). design/MOMENTUMPROJECTION.md.
- p^{2n} flows contain shape derivatives up to order 2n: an explicit shape substituted later must be C^{2n} (Litim is not C³;
  the KPZ tests' Litim helpers abort with `Litim.smoothShape`).
- The (iν)^m coefficient is (∂/∂(iν))^m / m! at ν = 0.

## Ansatz density (DSL)

- Γ_k = ∫dt dᵈx density. The density is written with the package-owned, protected, inert operators `dt`, `grad`, `div`, `lap`,
  plus `.` between two gradients, and `jac[a, b]` = ε_ij ∂_i a ∂_j b (d = 2 only, a and b scalars, at most one per monomial).
  With ∇ → −iq, `jac[A, B]` with A at leg l and B at leg m gives −`wedge[q_l, q_m]`.
  - Do not use System `Dt`, `D`, `Grad`, `Div` or `Laplacian`. `Dt[nk H]` silently adds `H Dt[nk]`, and the System operators
    applied to explicit coordinates evaluate to 0 before the library sees them.
- Couplings depend on k only. They stand outside the operators: write `nk dt[H]`, not `dt[nk H]`.
- Momentum functions (`"MomentumFunctions" -> {Fnu, ...}`, design/MOMENTUMFUNCTIONS.md) are written bare in quadratic terms, outside the
  operators, at most one per slot, as number × (a reference, or nothing) × F: `- nk Fnu hbar lap[h]` is −nk Fnu(p²) h̄ ∇²h.
- Each parsed monomial is coefficient × Π `FieldFactor[φ, n_t, {indices}]`. `FieldFactor[φ, n_t, {i, …}]` is ∂t^{n_t} ∂_i … φ.
  A repeated `SpaceIndex[1]` is summed over. The two indices of ε are `EpsIndex[1]` (on the factor of a) and `EpsIndex[2]` (on
  the factor of b).

## Reserved names

Fields and couplings may not be named `k q p y t x c ap ω ν`. These names have a fixed meaning in library input and output: RG
scale, loop momentum, external momentum, shape argument, time, position, the cosine c = cos(p, q) and the dimensionless external
momentum ap of the functional equations, frequencies. They may not be System symbols or
FRGscalarMSR symbols either, and they must have no definitions.

## Flows at generic p (momentum functions)

- Functional equations are written at generic external momentum p with the loop momentum q routed at the midpoint: the two
  propagator lines carry q − p/2 and q + p/2; the integrand is symmetrised under q → −q, so it is even in p at fixed q.
  Variables: q = |q|, p = |p|, c = cos(p, q) (dimensionless: aq, ap, c); a function is F[x] of the momentum squared x.
- ∫dᵈq/(2π)ᵈ = ∫dq ∫_{−1}^{1}dc S_{d−1} q^{d−1}(1−c²)^{(d−3)/2}/(2π)ᵈ, S_{d−1} = 2π^{(d−1)/2}/Γ((d−1)/2) (`FunctionalMeasure`);
  in d = 1 it is ∫dq/(2π) Σ_{c=±1}.
- A slot Σ_{n<n₀} c_n p^{2n} + p^{2n₀} ref F(p²) gives ∂s(ref F(p²)) = [loop(p) − its orders below n₀]/p^{2n₀}.
- F̂(ŷ) = F(k² ŷ)/[F], ŷ = ap²; ∂s at fixed ŷ adds 2ŷ F̂′(ŷ) to the linear part. F̂(0) is free unless `"Equation" -> ref F` normalises it to 1.

# LEVICIVITA — `jac` operator (2D Levi-Civita term): implementation hand-off

This document is a self-contained implementation plan for an agent with fresh context.
- The library's own design is in `design/FRGDESIGN.md`.
- This design survived two adversarial reviews: `design/reviews/LeviCivita-r1.md` and `LeviCivita-r2.md`. The final review
  found 0 blockers and 0 majors, and all of its findings are folded in below.
- Every line anchor refers to commit `1fbb12e`.

## 0. Before you start
- Read `README.md`, `docs/MODULES.md` (the sections of the modules you touch), `docs/CONVENTIONS.md` and
  `design/FRGDESIGN.md` §decisions.
- Test command: `math -script tests/runTests.wl`. Use `--verbose` to see the cards. Run it before any change, so you know the
  baseline is green.
- Library conventions you must follow:
  - every check is a named `FRGAssert` / `assertCheck`, and a failure aborts;
  - stages are immutable and fingerprinted;
  - exact zero tests use `ExactlyZeroQ`, never `Simplify[…] == 0`;
  - unit tests use `test[...]` / `expectFail[name, expr, "Check.name"]`.

## 1. Goal and constraints
- Users can write ε_ij ∂_iA ∂_jB in the density as `jac[A, B]`, with ε₁₂ = +1. The motivating term is cubic:
  c φ̄ ε_ij ∂_iA ∂_jB.
- `jac` exists only in d = 2. The user must **declare** it with `DefineTheory[…, "Dimension" -> FixedDimension[2]]`. d is
  never fixed silently.
- **No regression.** For a symbolic dimension and no `jac`, every stage is byte-identical to today: same keys, fingerprints,
  checks and card text. Gate every new code path on ε (or on the declared 2), and add stage keys only when ε is present.

## 2. Mathematics (all verified numerically in review r1)
- **F1.** wedge[a,b] ≡ ε_ij a_i b_j is bilinear and antisymmetric, and wedge[a,a] = 0.
- **F2 (d = 2 only).** wedge[a,b] wedge[c,d] = scal[a,c] scal[b,d] − scal[a,d] scal[b,c]. With a symbolic d the p²
  projection would be off by a factor of d − 1, which is why d must be 2.
- **F3.** Convention ∇ → −iq (CONVENTIONS.md): ε_ij ∂_iA ∂_jB at legs (l, m) → −wedge[q_l, q_m]. Γ2 is unaffected, because
  wedge[q, −q] = 0.
- **F4.** After F2, each term of the loop integrand is wedge[p,q]^e f(p², p·q, q², ν) with e ∈ {0, 1}. The e = 1 terms
  vanish at p = 0, and so does their p-Laplacian (2 ε_ij q_j q_i ∂_{p·q} f = 0). Dropping them is exact for every projection
  (m, n ≤ 1).
- **F5.** A quadratic jac term is a total derivative (∂_i(ε_ij A ∂_jB)), and its vertex vanishes.

## 3. Implementation, in this order (run the narrow tests after each step)

### Step 1: Theory, declared dimension (`Kernel/Theory.wl`)
- New public inert head `FixedDimension`, with a usage message.
- `DefineTheory` (Theory.wl:28-30) accepts either of two forms:
  - `"Dimension" -> Dim`, a symbol: unchanged, including the abort when the symbol has a value. A leftover `Dim = 2` must
    still abort;
  - `"Dimension" -> FixedDimension[2]`: the stage stores `"Dimension" -> 2`.
- `FixedDimension[n]` with n ≠ 2, or any other non-symbol value, aborts with `Theory.dimension`. Hint: "a symbol, or
  FixedDimension[2] (needed by jac)".
- `Theory.dimensionDistinct` and `assertFreshSymbols` run only for a symbol.
- Card line (Theory.wl:36): print "dimension 2 (declared)".
- Review r2 checked every downstream use of `"Dimension"` and found them all harmless with the integer 2. The Oracle rule
  `eqs["Dimension"] -> d` at Oracle.wl:26 becomes `2 -> 2`, which does nothing.
  - Still to fix: the Assemble card text "… symbolic" at Assemble.wl:128 must print "d = 2" when the dimension is 2.

### Step 2: Core (`Kernel/Core.wl`, next to `scal` at :133-139)
- New public symbol `wedge`, with a usage message. Rules:
  - linear over sums and explicit numbers in both slots;
  - `wedge[0,_] = wedge[_,0] = 0`;
  - `wedge[-a,b] = wedge[a,-b] = -wedge[a,b]`;
  - `wedge[a,a] = 0`;
  - canonical order `wedge[b,a] := -wedge[a,b]` when `! OrderedQ[{b,a}]`. Do **not** make it Orderless.
- Internal `reduceWedge2D[e]`, declared in the `FRGscalarMSR`Internal`` block: applies F2 to products and powers until each
  term has at most one wedge.
- `displayString` (:43): add `wedge[a, b] → wedge(a,b)`. Do not use `a^b`, which reads as a power.

### Step 3: Ansatz (`Kernel/Ansatz.wl`)
- **New heads.** New protected inert heads `jac` and `EpsIndex` (usage messages; `SetAttributes` at :23).
  - Add `jac` to `opHeads` (:25) and to the `badSyms` exclusion list (:129).
- **Operator checks.**
  - `Ansatz.operatorArity` (:136) becomes per head: `jac` takes 2 arguments, the others take 1.
  - `Ansatz.operatorOnField`: every argument must contain a field.
  - `Ansatz.couplingsOutsideOperators` (:137): test every argument, not only `First`, so that `jac[H, nk Psi]` fails.
- **`Ansatz.epsilonDimension`.** New check: if the density contains `jac`, then `theory["Dimension"] === 2`. Otherwise abort
  with the hint: "jac (ε_ij) exists only in d = 2: declare \"Dimension\" -> FixedDimension[2]".
- **Rank rule.** `rankOf[jac[a_, b_], fld_]`: both arguments must have rank 0 (`needRank`), and the result is 0.
- **Parsing (`toFactors`, near :67):**
  `S[jac[a_, b_]] := With[{j = new[]}, Expand[derive[EpsIndex[j, 1]][S[a]] derive[EpsIndex[j, 2]][S[b]]]]`.
  - `derive[i_]` already appends any index that is not "t".
- **`Ansatz.singleEpsilon`.** New check, run on every term after `splitTerm` (:84) and **before** the `GroupBy` (:150).
  - Count the occurrences of `EpsIndex[_, 1]` in the term's flattened factor list. `splitTerm` already expands powers into
    repeated factors.
  - The count must be ≤ 1. Otherwise abort with the hint "at most one jac per monomial (p^4 is unsupported)".
  - Count slots, not distinct labels: `jac[..]^2` repeats one label (review r2-E1).
- **`canonical` (:81).** It now returns `{sign, list}`:
  - Rename `SpaceIndex` as today, excluding `EpsIndex`. Then map `EpsIndex[_, s] → EpsIndex[s]`.
  - Let L be the sorted list and L′ the sorted list with slots 1 ↔ 2 swapped.
    - If `L′ === L`, the monomial is zero.
    - Otherwise take the smaller of L and L′ in canonical order, with sign −1 when L′ is chosen.
  - Change the `GroupBy` (:150) to sum sign × coefficient.
  - ε-free monomials get sign +1 and exactly today's list.
- **`Ansatz.jacQuadratic`.** New check: a monomial with ε and `FieldCount == 2` aborts with the hint: "a quadratic jac term
  is a total derivative and does not contribute: remove it".
- **`Ansatz.indexPairs` (:164).** Every `SpaceIndex` appears twice; `EpsIndex[1]` and `EpsIndex[2]` appear once each or are
  both absent.
- **Stage key.** Add `"Epsilon" -> True` to the stage only when ε is present. Read it everywhere with
  `Lookup[ansatz, "Epsilon", False]`.
- **Display.** `fieldString` / `MonomialString` (:90-96) print the ε pair as `jac[f_slot1, f_slot2]`, with time derivatives
  wrapped as `dt[...]`.

### Step 4: Vertices (`Kernel/Vertices.wl`)
- `factorSymbol` / `monomialAt` (:27-30):
  - `SpaceIndex` groups → `scal`, as today;
  - the `EpsIndex[1]` and `EpsIndex[2]` momenta → `wedge[q_slot1, q_slot2]`;
  - the `(−I)^Length[ix]` prefactor already counts both indices.
- Plane-wave check (`planeWaveCoefficient` :61, `symbolAsPlaneWave` :71, card :117):
  - component count `nc` = 2 when `theory["Dimension"] === 2`, otherwise 3 (as today);
  - with ε: `jac -> (D[#1, x1] D[#2, x2] − D[#1, x2] D[#2, x1] &)` on the raw density;
  - `wedge[a,b]` maps to a₁b₂ − a₂b₁;
  - the card prints `nc`.
- The Γ2 checks are unchanged.

### Step 5: Projection (`Kernel/Projection.wl`)
- `dim` is the declared value (:45). If ε is present, assert `dim === 2`.
- In `project` (:26), as the first operation when ε is present: `e = reduceWedge2D[Expand[e]]`.
- New check `Projection.wedgeOdd`:
  - every remaining wedge is `±wedge[p, qL]`, to power ≤ 1;
  - drop that part (F4);
  - as evidence, evaluate the dropped part's n = 0 and n = 1 projections with explicit 2-vectors at a random rational point,
    and require exactly 0.
- `Projection.complete` (:35) also excludes `wedge`.

### Step 6: Assemble and Oracle
- **Assemble.** The card (:128) prints "d = 2" for a declared 2. The `"Dimension"` key carries the declared value.
- **Oracle (`Kernel/Oracle.wl`).**
  - `"Dimensions" -> Automatic` (:37). Resolve it **before** any use, including the card text at :96: {2, 3} for a symbolic
    dimension (today's card text unchanged), {2} for a declared 2.
  - An explicit list must be a subset of {2, 3} when d is symbolic, and equal to {2} when d is declared. Otherwise abort with
    `Oracle.declaredDimension`.
  - `dot[]` (:64) also maps `wedge[a,b]` to a₁b₂ − a₂b₁. Apply it after substituting the numeric `pVec`, to avoid
    `Part::partd`.

### Step 7: Docs
- `docs/CONVENTIONS.md`: `jac`, ε₁₂ = +1, `FixedDimension[2]`.
- `docs/MODULES.md`: Theory, Core, Ansatz, Vertices, Projection, Oracle. Cover the new checks and F2/F4/F5.
- `design/FRGDESIGN.md`: the exception to E12/Q16 (symbolic d, Oracle in d = 2 and 3) for a declared d = 2.
- `docs/FUTURE.md`: the 3D ε_ijk (curl-type) operator.

## 4. Tests
- **`tests/unit/TheoryTest.wl`:**
  - `FixedDimension[2]` is accepted and the stage stores 2;
  - `FixedDimension[3]` → `Theory.dimension`;
  - a leftover `Dim = 2` still aborts.
- **`tests/unit/CoreTest.wl`:** wedge antisymmetry, linearity in both slots, the zero rules; `reduceWedge2D` against
  explicit random 2-vectors.
- **`tests/unit/AnsatzTest.wl`:**
  - `jac` with a symbolic `Dim` → `Ansatz.epsilonDimension`;
  - `Hbar (jac[A,B] + jac[B,A])` → all monomials cancel;
  - `Hbar jac[A,A]` → cancelled;
  - `jac[dt[A], B]` is parsed;
  - `jac[grad[A], B]` → `Ansatz.vectorRank`;
  - `jac[A]` → `Ansatz.operatorArity`;
  - `jac[H, nk Psi]` → `Ansatz.couplingsOutsideOperators`;
  - `Hbar jac[A,B] jac[B,A]` and `Hbar jac[H+Psi, H-Psi]^2` → `Ansatz.singleEpsilon`;
  - `c jac[Hbar, H]` → `Ansatz.jacQuadratic`;
  - `MonomialString` shows `jac[…]`;
  - the fingerprint of an ε-free, symbolic-Dim ansatz is unchanged. Record it on the baseline first.
- **`tests/unit/VerticesTest.wl`:** for `c Hbar jac[H, Psi]`, Γ3_{Hbar H Psi} = −c wedge[q_H, q_Psi] (hand value, F3). The 2D
  plane-wave check passes.
- **`tests/unit/OracleTest.wl` (new):** `{2, 3}` with a declared 2 → `Oracle.declaredDimension`.
- **`tests/acceptance/Jacobian2D.wl` (new):**
  - Density:
    `Hbar (dt[H] + lk jac[Psi, H] + kap grad[Psi].grad[H] − nk lap[H]) − dk Hbar^2 + Psibar (dt[Psi] − mk lap[Psi]) − ek Psibar^2`.
    Do not use the names `D`, `E` or `ν`: they are System or reserved symbols.
  - Setup modelled on `tests/reference/KPZSetup.wl`: regulator on H̄H, H̄H̄, Ψ̄Ψ and Ψ̄Ψ̄, `TruncationOrders`, and
    Protected/Neglected tags for every monomial not in the ansatz.
  - The `kap` term is essential: it creates single-wedge cross terms (a non-trivial F4) and nonzero n = 0 flows.
  - Run the full pipeline, then `RunOracle` in d = 2 within 10⁻⁸.
- **Declared-2 regression (acceptance):** shallow water with `FixedDimension[2]`. Each equation must satisfy
  `ExactlyZeroQ[new − (old /. Dim -> 2)]`. Do not use `===`, because Assemble picks its output form by LeafCount.
- **Done when** `math -script tests/runTests.wl` exits 0 with all old and new tests, and the docs are updated.

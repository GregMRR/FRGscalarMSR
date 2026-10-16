# Legacy equations (old notebooks)

These files are used by `tests/acceptance/ShallowWaterLPA.wl` and `tests/acceptance/ShallowWaterDE4Legacy.wl`. They were
copied read-only from `~/sync_cnrs/Projects/shallow_water/` on 2026-10-16, with `ϕ` in file names renamed to `phi`; every copy
is byte-identical to its source.

## Folders

### `flowsDE4plus_agreed/`

- **Source:** `flowsDE4plus_agreed/`. It is identical to `LPA_KPZ_thermal_forcing/flowsDE4plus_agreed/`, the input of
  `AdimensionalFlowEQGenericReg.nb`.
- **Content:** DE4 flows in the old `Subscript` notation, already normalised (∂s of the coupling), with aDhk = ank = aDphik = 1
  implicit. `AdimensionalFlowEQGenericReg.nb` applies `DE4toLPA` = {K, D2ϕ, g2 → 0} to them for the LPA.

### `flowsLPAold/`

- **Source:** `Legacy-backup/flowsLPAold/`.
- **Content:** the oldest LPA set, with λ = 1, A = 0 and no regulator on the H̄H̄ noise (no rm_h, no dsMh). The depth is written
  H_0.
- **Normalisation:** raw right-hand sides. The flows are 1 × gEQ, −1/2 × DϕEQ, 1/(2 Dim) × NEQ and −1/(4 Dim) × DhEQ, the
  factors of `shallow_water_DE4plus.nb` §14.

### `flowsNEW/`

- **Source:** `flowsNEW/`, the exports of `shallow_water_DE4plus.nb`.
- **Content:** new notation, same normalisation as flowsFAST. `DphiD2phiEQ` is **not copied**: the source file is corrupt (it
  contains `Part[...]`; see AGENT_PRIMER §6).

## Symbol map

`OldNotation.wl` holds the map, taken from `OldtoNew` and `DimensionlessRegulators` of `shallow_water_DE4plus.nb`:

| old | new |
|---|---|
| N_k | nk |
| λ_k | lk |
| K_k | Kk |
| A_k | Ak |
| D_hk | Dhk |
| D_ϕk | Dphik |
| D_2ϕk | D2phik |
| g_k | gk |
| g_2k | g2k |
| H_k, H_0 | Hk |
| rn_h(q²) | arnh(q²/k²); likewise rn_ϕ, rm_h, rm_ϕ, with 1/k² on first derivatives |
| dsNϕ, dsMϕ | dsNphi, dsMphi |

## Results (2026-10-16)

**Exact agreement.**
- `flowsDE4plus_agreed`: g, Dϕ, N, Dh, A and K, in full DE4. In LPA they also agree, with {K, D2ϕ, g2 → 0}.
- `flowsLPAold`: g, Dϕ, N and Dh, with the factors above.
- `flowsNEW`: g, Dϕ, N, Dh, A, K, H (≡ 0) and g2.

**Not usable as references.**
- `flowsDE4plus_agreed/HEQ.wdx` holds only the placeholder symbol `HkEQ`.
- `flowsDE4plus_agreed/g2EQ.wdx` and `DphiD2phiEQ.wdx` contain unsubstituted second derivatives f″(q²) of fη, fγ, fDh and fDϕ.
  - Even after completing them with the DE4 forms fη = nk s (Ak + arnh), fγ = nk s (1 + arnphi), fDh = 2 Dhk s (1 + armh) and
    fDϕ = 2 Dphik (1 + armphi + D2phik s), they **disagree** with this library.
  - The ratio varies from point to point, every ∂sR component differs, and the difference does not vanish at Dim = 1.
  - This library's g2 and D2ϕ agree with `flowsNEW` (g2), with `flowsFAST` (both), and with the numerical oracle in
    d = 2, 3 (both, to 2·10⁻⁹).
  - The old agreed g2 and D2ϕ are therefore wrong, unless they were derived for a different ansatz.

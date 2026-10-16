# Legacy LPA equations (old notebooks), for `tests/acceptance/ShallowWaterLPA.wl`

Copied read-only from `~/sync_cnrs/Projects/shallow_water/` on 2026-10-16; `ϕ` in file names renamed to `phi`.

| folder | source | content |
|---|---|---|
| `flowsDE4plus_agreed/` | `LPA_KPZ_thermal_forcing/flowsDE4plus_agreed/` (gEQ, DϕEQ, NEQ, DhEQ, ANEQ) | DE4 flows in the old `Subscript` notation, already normalised (∂s of the coupling), aDhk = ank = aDphik = 1 implicit. Input of `AdimensionalFlowEQGenericReg.nb`, which applies `DE4toLPA` = {K, D2ϕ, g2 → 0} and exports the Julia LPA equations. |
| `flowsLPAold/` | `Legacy-backup/flowsLPAold/` (gEQ, DϕEQ, NEQ, DhEQ) | The oldest LPA set: λ = 1, A = 0, no regulator on the H̄H̄ noise (no rm_h, dsMh); depth written H_0. Raw RHS: the flows are 1 × gEQ, −1/2 × DϕEQ, 1/(2 Dim) × NEQ, −1/(4 Dim) × DhEQ (factors of `shallow_water_DE4plus.nb` §14). |

Symbol map (from `OldtoNew` and `DimensionlessRegulators` of `shallow_water_DE4plus.nb`): N_k → nk, λ_k → lk, K_k → Kk,
A_k → Ak, D_hk → Dhk, D_ϕk → Dphik, D_2ϕk → D2phik, g_k → gk, g_2k → g2k, H_k, H_0 → Hk, rn_h(q²) → arnh(q²/k²) (and
rn_ϕ, rm_h, rm_ϕ alike; first derivatives get 1/k²), dsNϕ → dsNphi, dsMϕ → dsMphi.

(* ::Package:: *)

(* old notebook notation -> FRGscalarMSR / DE4plus_fast.wl notation, shared by the legacy acceptance tests.
   OldtoNew and DimensionlessRegulators copied from shallow_water_DE4plus.nb, plus H_0 = H_k and dsN\[Phi] = dsNphi, dsM\[Phi] = dsMphi. *)

oldToNew = {Subscript[N, k] -> nk, Subscript[\[Lambda], k] -> lk, Subscript[K, k] -> Kk, Subscript[A, k] -> Ak, Subscript[D, hk] -> Dhk,
   Subscript[D, \[Phi]k] -> Dphik, Subscript[D, 2 \[Phi]k] -> D2phik, Subscript[g, k] -> gk, Subscript[g, 2 k] -> g2k, Subscript[H, k] -> Hk, Subscript[H, 0] -> Hk,
   Subscript[rn, h][q^2] -> arnh[q^2/k^2], Subscript[rn, \[Phi]][q^2] -> arnphi[q^2/k^2], Subscript[rm, h][q^2] -> armh[q^2/k^2], Subscript[rm, \[Phi]][q^2] -> armphi[q^2/k^2],
   Derivative[1][Subscript[rn, h]][q^2] -> Derivative[1][arnh][q^2/k^2]/k^2, Derivative[1][Subscript[rn, \[Phi]]][q^2] -> Derivative[1][arnphi][q^2/k^2]/k^2,
   Derivative[1][Subscript[rm, h]][q^2] -> Derivative[1][armh][q^2/k^2]/k^2, Derivative[1][Subscript[rm, \[Phi]]][q^2] -> Derivative[1][armphi][q^2/k^2]/k^2,
   dsN\[Phi] -> dsNphi, dsM\[Phi] -> dsMphi};
newNotationSymbols = {nk, lk, Kk, Ak, Dhk, Dphik, D2phik, gk, g2k, Hk, k, q, Dim, dsNh, dsMh, dsNphi, dsMphi, arnh, arnphi, armh, armphi, aDhk, ank, aDphik};
legacyRef[set_, f_] := Import[FileNameJoin[{testRoot, "reference", "legacy", set, f <> ".wdx"}]] /. oldToNew;
(* explicit equation with d_s R abstract under the script's names (E10) *)
legacyE10 = {dsR[H, Hbar] -> dsNh, dsR[Hbar, Hbar] -> dsMh, dsR[Psi, Psibar] -> dsNphi, dsR[Psibar, Psibar] -> dsMphi};
legacyForm[eqs_, c_] := EquationOf[eqs, c]["Compact"] /. Normal[eqs["Atoms"]] /. legacyE10;

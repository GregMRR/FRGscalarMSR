(* r2-probe-related: the D12 example Related[FD -> al Fnu, "When" -> {nk -> al Dk, rn -> al rD}] against the time-reversal
   (FDT) condition used by the KS benchmark (Checks.wl:20,32: nu = al D, tau = al zeta, r_nu = al r_D), written on Gamma2 + R:
     Re(Gamma11 + R11)/p^2 == T * ( -(Gamma02 + R02)/2 )   for every p, with one constant T.
   LO ansatz: Gamma11 = I w + nk Fnu(p^2) p^2, Gamma02 = -2 Dk FD(p^2). Two ways of writing the regulator. *)
cond[r11_, r02_] := Module[{lhs, rhs},
   lhs = nk Fnu[y] + r11; rhs = T (Dk FD[y] - r02/2);
   lhs - rhs];
plan = {FD -> (al Fnu[#] &), nk -> al Dk, rn -> (al rD[#] &)};
(* KS-style regulator: R11 = Dk p^2 rn, R02 = -2 Dk rD  (benchmarks/KuramotoSivashinsky/Pipeline.wl:21) *)
ks = cond[Dk rn[y], -2 Dk rD[y]] /. FD[y] -> FD[y] /. plan;
(* KPZ-style regulator: R11 = nk p^2 rn, R02 = -2 Dk rd  (tests/reference/KPZSetup.wl) *)
kpz = cond[nk rn[y], -2 Dk rD[y]] /. plan;
solveT[e_] := Module[{cs = CoefficientList[Expand[e /. {Fnu[y] -> u, rD[y] -> v}], {u, v}] // Flatten // DeleteCases[0]}, Solve[Thread[cs == 0], {T}]];
Print["KS regulator:  conditions ", Union[Flatten[CoefficientList[Expand[ks /. {Fnu[y] -> u, rD[y] -> v}], {u, v}]]], " -> T: ", solveT[ks]];
Print["KPZ regulator: conditions ", Union[Flatten[CoefficientList[Expand[kpz /. {Fnu[y] -> u, rD[y] -> v}], {u, v}]]], " -> T: ", solveT[kpz]];
Print["with al = 1:   ", solveT[ks /. al -> 1], solveT[kpz /. al -> 1]];
(* consistent triples: KS regulator -> {FD -> Fnu, nk -> al Dk, rn -> al rD} ; KPZ regulator -> {FD -> Fnu, nk -> al Dk, rn -> rD} *)
okKS = cond[Dk rn[y], -2 Dk rD[y]] /. {FD -> (Fnu[#] &), nk -> al Dk, rn -> (al rD[#] &)};
okKPZ = cond[nk rn[y], -2 Dk rD[y]] /. {FD -> (Fnu[#] &), nk -> al Dk, rn -> (rD[#] &)};
Print["FD -> Fnu, KS regulator, rn -> al rD: T ", solveT[okKS], "   KPZ regulator, rn -> rD: T ", solveT[okKPZ]];

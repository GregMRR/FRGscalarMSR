(* r2-probe-flam: today's behaviour for a coupling in front of hbar dt[h] (the slot of f^lambda in LO/NLO KPZ), as a stand-in
   for a MomentumFunction there: D6 gives it no reference prefactor, D5 forbids declaring a function "Dimensionless". *)
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
th = Block[{Print}, DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
an = Block[{Print}, DefineAnsatz[th, "Density" -> hbar (Flam dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {Flam, nuk, lamk, Dk}]];
r = CheckAbort[DeriveDimensions[an, "References" -> {nuk -> "nu", Dk -> "D"}], "ABORTED"];
Print["result: ", If[r === "ABORTED", r, "accepted"]];

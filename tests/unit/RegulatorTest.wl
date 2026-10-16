(* ::Package:: *)

(* Regulator: Gamma2-level entries, d_s R by the chain rule, every assertion *)

rTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
rAnsatz = quiet[DefineAnsatz[rTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
reg = quiet[DefineRegulator[rAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}]];

test["R matrix (symmetric, pair order irrelevant)", RegulatorMatrix[reg, s0], {{0, nuk s0 rn[s0/k^2]}, {nuk s0 rn[s0/k^2], -2 Dk rd[s0/k^2]}}];
test["d_s R of nu q^2 r(q^2/k^2)", Expand[DsRegulatorMatrix[reg, s0][[1, 2]] - (dsCoupling[nuk] s0 rn[s0/k^2] - 2 nuk s0 (s0/k^2) rn'[s0/k^2])], 0];
test["d_s R of -2 D r(q^2/k^2)", Expand[DsRegulatorMatrix[reg, s0][[2, 2]] - (-2 dsCoupling[Dk] rd[s0/k^2] + 4 Dk (s0/k^2) rd'[s0/k^2])], 0];
test["abstract dsR matrix", DsRAbstractMatrix[reg], {{0, dsR[h, hbar]}, {dsR[h, hbar], dsR[hbar, hbar]}}];
test["dsCoupling symbols", reg["DsCouplings"], {dsCoupling[Dk], dsCoupling[nuk]}];
test["scale symbol is the user's k", reg["Scale"], k];

badReg[args___] := DefineRegulator[rAnsatz, args];
expectFail["physical-physical pair", badReg[{h, h} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}], "Regulator.responsePair"];
expectFail["pair given twice", badReg[{hbar, h} -> nuk q^2 rn[q^2/k^2], {h, hbar} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}], "Regulator.pairOnce"];
expectFail["odd power of q", badReg[{hbar, h} -> nuk q rn[q^2/k^2], "Shapes" -> {rn}], "Regulator.evenInQ"];
expectFail["shape with another argument", badReg[{hbar, h} -> nuk q^2 rn[q^2/k], "Shapes" -> {rn}], "Regulator.shapeArgument"];
expectFail["shape derivative in the input", badReg[{hbar, h} -> nuk q^2 rn'[q^2/k^2], "Shapes" -> {rn}], "Regulator.shapeArgument"];
expectFail["undeclared function", badReg[{hbar, h} -> nuk q^2 Exp[-q^2/k^2] zf[q^2/k^2], "Shapes" -> {rn}], "Regulator.onlyDeclaredShapes"];
expectFail["unknown symbol (e.g. a frequency)", badReg[{hbar, h} -> nuk q^2 rn[q^2/k^2] + wz, "Shapes" -> {rn}], "Regulator.knownSymbols"];
expectFail["shapes required", badReg[{hbar, h} -> nuk q^2], "Regulator.shapeList"];
expectFail["shape named like a coupling", badReg[{hbar, h} -> nuk q^2 Dk[q^2/k^2], "Shapes" -> {Dk}], "Regulator.shapesDistinct"];
expectFail["unknown field in a pair", badReg[{hbar, zz} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}], "Regulator.pairFields"];
expectFail["no scale k", badReg[{hbar, h} -> nuk q^2, "Shapes" -> {rn}], "Regulator.momentumSymbol"];
expectFail["wrong input stage", DefineRegulator[rTheory, {hbar, h} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}], "Stage.input"];

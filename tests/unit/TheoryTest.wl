(* ::Package:: *)

(* Theory: response fields, field order, collisions *)

Module[{th = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]]},
  test["two scalars: field order", th["Fields"], {H, Hbar, Psi, Psibar}];
  test["two scalars: response fields", th["ResponseFields"], {Hbar, Psibar}];
  test["response field of Psi", ResponseField[th, Psi], Psibar];
  test["stage key", th["Stage"], "Theory"]];

test["one scalar (KPZ)", quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]]["Fields"], {h, hbar}];

expectFail["N = 3 aborts (deg Delta <= 2)", DefineTheory["Scalars" -> {f1, f2, f3}, "Dimension" -> Dim], "Theory.N"];
expectFail["reserved name k", DefineTheory["Scalars" -> {k}, "Dimension" -> Dim], "Theory.scalars.reserved"];
expectFail["library name collision", DefineTheory["Scalars" -> {scal}, "Dimension" -> Dim], "Theory.scalars.library"];
expectFail["System symbol", DefineTheory["Scalars" -> {Pi}, "Dimension" -> Dim], "Theory.scalars.system"];
hh = 3;
expectFail["scalar with a value", DefineTheory["Scalars" -> {hh}, "Dimension" -> Dim], "Theory.scalars.symbols"];
ClearAll[hh];
hf[x_] := x;
expectFail["scalar with a definition", DefineTheory["Scalars" -> {hf}, "Dimension" -> Dim], "Theory.scalars.values"];
ClearAll[hf];
expectFail["missing dimension", DefineTheory["Scalars" -> {H}], "Theory.dimension"];
expectFail["dimension equal to a field", DefineTheory["Scalars" -> {H}, "Dimension" -> Hbar], "Theory.dimensionDistinct"];
expectFail["duplicate scalars", DefineTheory["Scalars" -> {H, H}, "Dimension" -> Dim], "Theory.scalars.duplicates"];
expectFail["unknown option", DefineTheory["Scalars" -> {H}, "Dimension" -> Dim, "Foo" -> 1], "Theory.options"];

(* declared dimension (jac, design/LEVICIVITA.md) *)
test["FixedDimension[2] accepted: the stage stores 2", quiet[DefineTheory["Scalars" -> {H}, "Dimension" -> FixedDimension[2]]]["Dimension"], 2];
expectFail["FixedDimension[3] aborts", DefineTheory["Scalars" -> {H}, "Dimension" -> FixedDimension[3]], "Theory.dimension"];
expectFail["a literal 2 aborts", DefineTheory["Scalars" -> {H}, "Dimension" -> 2], "Theory.dimension"];
Dim = 2;
expectFail["a leftover Dim = 2 still aborts", DefineTheory["Scalars" -> {H}, "Dimension" -> Dim], "Theory.dimension"];
ClearAll[Dim];

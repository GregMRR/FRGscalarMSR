(* ::Package:: *)

(* Dimensions.wl -- DeriveDimensions: the scaling dimension of every coupling, field and of time, derived from the ansatz under
   the assumption that Gamma_k = Integral dt d^dx density has no dimension, and written in a basis chosen by the user:
   the RG scale k and N+1 reference couplings.
   Objective, assumptions and checks: docs/MODULES.md section "Dimensions". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

DeriveDimensions::usage = "DeriveDimensions[ansatz, \"References\" -> {nk -> \"nu\", Dhk -> \"h\", Dphik -> \"phi\"}, \
\"Dimensionless\" -> {Ak, ...}] derives the dimension of every coupling, field and of time as k^x ref1^y1 ... ref_{N+1}^y_{N+1} \
(N+1 references, N = number of scalars) and returns the \"Dimensions\" stage. A reference may be written nk -> {\"nu\", \"Mode\" -> \"Running\"} \
or nk -> {\"nu\", \"Mode\" -> \"Running\", \"Equation\" -> nk ank} (see AdimensionalEquations). Options: \"Names\" -> <|gk -> ghat|> (dimensionless names), \"Scale\" -> k.";
ShowDimensions::usage = "ShowDimensions[dims] prints the dimension table: coupling | monomials it appears in | dimension | dimensionless \
name | linear part of its flow; then the fields, time and the frequency exponent z.";
CheckDimensions::usage = "CheckDimensions[dims, <|gk -> k Dphik^(1/2) Dhk^(-1/2) nk, H -> ..., \"Frequency\" -> k^2 nk|>] compares \
hand-derived dimensions with the derived ones exactly; any mismatch prints a FAIL block and aborts. Returns the list of compared keys.";
DimensionOf::usage = "DimensionOf[dims, x] is the dimension of the coupling or field x (or \"Frequency\", \"Time\") as k^x ref^y.";
LinearPart::usage = "LinearPart[dims, c] is the canonical part of the flow of the dimensionless coupling: -x + Sum_j y_j eta_j \
for c of dimension k^x Prod_j ref_j^y_j, so that d_s chat = chat LinearPart + (d_s c)/[c] (d_s ref_j = -eta_j ref_j).";
AnsatzOf::usage = "AnsatzOf[stage] returns the Ansatz stage embedded in a later stage (e.g. an Equations stage read by LoadStage).";
EmbeddedStage::usage = "EmbeddedStage[stage, name] returns the stage called name embedded (at any depth) in stage, or Missing.";

Begin["`Private`"];

(* ---------- embedded stages ---------- *)

EmbeddedStage[s_Association, name_String] := If[Lookup[s, "Stage", None] === name, s,
   FirstCase[EmbeddedStage[#, name] & /@ Select[Values[s], AssociationQ[#] && KeyExistsQ[#, "Stage"] &], _Association, Missing["NotEmbedded", name]]];
AnsatzOf[s_Association] := Block[{$FRGModule = "Dimensions"}, Module[{a = EmbeddedStage[s, "Ansatz"]},
    FRGAssert["Dimensions.ansatzOf", ! MissingQ[a], Lookup[s, "Stage", s], "this stage does not embed an Ansatz stage"]; a]];

(* ---------- dimension vectors ----------
   A dimension is a vector v over the basis (k, ref_1, ..., ref_{N+1}): the quantity scales as k^v1 ref_1^v2 ...
   [x] = 1/k, so a spatial derivative (one index: SpaceIndex or an EpsIndex slot) adds +e_k; a time derivative adds -[t]. *)

basisForm[v_List, basis_List] := Times @@ MapThread[Power, {basis, v}];
linearPart[v_List, etas_List] := -First[v] + Rest[v].etas;

(* the power of each coupling in one term of a coefficient: n_c = c dterm/dc / term, a number for a monomial *)
couplingPowers[term_, cps_List] := Module[{n = Together[# D[term, #]/term] & /@ cps},
   FRGAssert["Dimensions.monomialCoefficient", AllTrue[n, NumberQ], term,
    "each term of a coefficient must be a number times a product of powers of couplings"];
   n];
termsOf[c_] := With[{e = Expand[c]}, If[Head[e] === Plus, List @@ e, {e}]];

(* one equation per monomial term: [term] + Sum[fields] - nt [t] + n_x e_k + [t] - Dim e_k = 0 (Gamma has no dimension) *)
monomialEquations[mons_, cps_, fields_] := Flatten[Function[m, Function[t,
      <|"Monomial" -> MonomialString[m], "Term" -> t, "CouplingPowers" -> couplingPowers[t, cps],
       "FieldCounts" -> (Count[m["Fields"], #] & /@ fields), "TimeOrder" -> m["TimeOrder"], "MomentumOrder" -> m["MomentumOrder"]|>] /@
     termsOf[m["Coefficient"]]] /@ mons, 1];

(* unknown columns: time, fields, couplings; each equation is  row . u + const = 0, one copy per basis component *)
equationRow[eq_, cols_, cps_, fields_] := Join[{1 - eq["TimeOrder"]}, eq["FieldCounts"], eq["CouplingPowers"]][[cols]];

(* ---------- independent check: rescale the raw density ----------
   New units: every quantity of dimension v is multiplied by s(v) = Prod_j lambda_j^v_j. With t_new = s([t]) t, x_new = x/lambda_k,
   phi_new(t, x) = s([phi]) phi(t/s([t]), lambda_k x), c_new = s([c]) c, the density must obey
   L_new(t, x) = s(Dim e_k - [t]) L(t/s([t]), lambda_k x), i.e. the action Integral dt d^dx L is unchanged. Evaluated on explicit
   polynomial test fields with exact rationals, from the density as written (not from the parsed monomial table). *)
evalDensity[e_, fieldVals_, tt_, xs_] := Module[{ev},
   ev[s_Symbol] := Lookup[fieldVals, s, s];
   ev[dt[a_]] := D[ev[a], tt];
   ev[fn[f_, a_]] := ev[a];   (* LO: f is dimensionless with f(0) = 1, so fn[f, X] scales as X *)
   ev[grad[a_]] := D[ev[a], #] & /@ xs;
   ev[div[a_]] := Total[MapThread[D, {ev[a], xs}]];
   ev[lap[a_]] := Total[D[ev[a], {#, 2}] & /@ xs];
   ev[jac[a_, b_]] := With[{A = ev[a], B = ev[b]}, D[A, xs[[1]]] D[B, xs[[2]]] - D[A, xs[[2]]] D[B, xs[[1]]]];
   ev[Dot[a_, b_]] := ev[a] . ev[b];
   ev[x_Plus] := Plus @@ (ev /@ List @@ x);
   ev[x_Times] := Times @@ (ev /@ List @@ x);
   ev[Power[a_, n_]] := ev[a]^n;
   ev[x_] := x;
   ev[e]];

rescalingCheckQ[dens_, fields_, cps_, dimOf_, timeDim_, dimSym_, d_, seed_] := Module[{tt, xs, vec, den, lam, s, cv, fv, fvNew, cvNew, pt, lhs, rhs, deg = 3},
   xs = Array[Unique["x"] &, d]; tt = Unique["t"];
   vec[v_] := v /. dimSym -> d;
   den = LCM @@ Denominator[Flatten[vec /@ Join[Values[dimOf], {timeDim}]]];
   SeedRandom[seed];
   lam = (RandomInteger[{2, 9}, Length[timeDim]]/RandomInteger[{2, 9}, Length[timeDim]])^den;
   s[v_] := Times @@ (lam^vec[v]);
   cv = AssociationThread[cps -> RandomInteger[{1, 30}, Length[cps]]/RandomInteger[{1, 30}, Length[cps]]];
   fv = AssociationThread[fields -> Table[Sum[RandomInteger[{-5, 5}]/RandomInteger[{1, 7}] tt^a Times @@ (xs^RandomInteger[{0, deg}, d]), {a, 0, deg}, {3}],
       Length[fields]]];
   fvNew = AssociationThread[fields -> (s[dimOf[#]] (fv[#] /. Prepend[Thread[xs -> First[lam] xs], tt -> tt/s[timeDim]]) & /@ fields)];
   cvNew = AssociationThread[cps -> (s[dimOf[#]] cv[#] & /@ cps)];
   pt = Thread[Prepend[xs, tt] -> RandomInteger[{1, 9}, d + 1]/RandomInteger[{1, 9}, d + 1]];
   lhs = evalDensity[dens /. Normal[cvNew], fvNew, tt, xs] /. pt;
   rhs = s[UnitVector[Length[timeDim], 1] dimSym - timeDim] (evalDensity[dens /. Normal[cv], fv, tt, xs] /. Prepend[Thread[xs -> First[lam] xs], tt -> tt/s[timeDim]] /. pt);
   Together[lhs - rhs] === 0];

(* ---------- DeriveDimensions ---------- *)

Options[DeriveDimensions] = {"References" -> {}, "Dimensionless" -> {}, "Names" -> <||>, "Scale" :> Symbol["Global`k"]};

DeriveDimensions[ansatz_, opts : OptionsPattern[]] := Block[{$FRGModule = "Dimensions"},
  Module[{t0 = AbsoluteTime[], theory, nS, cps, fields, mons, dimSym, refSpec, refs, refNames, refModes, refEqs, dl, kk, basis, nb,
    eqs, unknownCps, cols, colNames, mFull, nsFull, mU, consts, sol, vec, dimOf, timeDim, checks = {}, warns = {}, bad, names,
    etaSyms, stage, ns, refRows, undet, conflicts, residual, nUnk, generated, aq, nameOverride},
   validateStage[ansatz, "Ansatz"];
   FRGAssert["Dimensions.options", MatchQ[{opts}, {(Rule | RuleDelayed)["References" | "Dimensionless" | "Names" | "Scale", _] ...}], {opts},
    "only the options \"References\", \"Dimensionless\", \"Names\", \"Scale\" are accepted"];
   theory = ansatz["Theory"]; nS = Length[theory["Scalars"]]; cps = ansatz["Couplings"]; fields = theory["Fields"];
   mons = ansatz["Monomials"]; dimSym = theory["Dimension"]; kk = OptionValue["Scale"];
   dl = OptionValue["Dimensionless"]; nameOverride = OptionValue["Names"];

   (* references: nk -> "nu" or nk -> {"nu", "Mode" -> "Running"} *)
   refSpec = OptionValue["References"];
   FRGAssert["Dimensions.referenceSyntax", ListQ[refSpec] &&
     AllTrue[refSpec, MatchQ[#, _Symbol -> (_String | {_String, ((("Mode" -> ("Running" | "Fixed")) | ("Equation" -> _)) ...)})] &], refSpec,
    "\"References\" -> {nk -> \"nu\", Dhk -> \"h\", ...}; a reference times an amplitude may be nk -> {\"nu\", \"Mode\" -> \"Running\", \"Equation\" -> nk ank}"];
   refs = refSpec[[All, 1]];
   refNames = Replace[refSpec[[All, 2]], {s_String :> s, {s_String, ___} :> s}, {1}];
   refModes = Replace[refSpec[[All, 2]], {_String -> Automatic, {_String, o___} :> Lookup[{o}, "Mode", Automatic]}, {1}];
   refEqs = Replace[refSpec[[All, 2]], {_String -> Automatic, {_String, o___} :> Lookup[{o}, "Equation", Automatic]}, {1}];
   AppendTo[checks, assertCheck["Dimensions.referenceCount", Length[refs] == nS + 1, refs,
     "choose exactly N+1 = " <> ToString[nS + 1] <> " reference couplings: N for the field dimensions, 1 for the time dimension",
     ToString[nS + 1] <> " references for N = " <> ToString[nS] <> ": " <> StringRiffle[ToString /@ refs, ", "]]];
   FRGAssert["Dimensions.referencesAreCouplings", SubsetQ[cps, refs] && DuplicateFreeQ[refs], refs, "references must be distinct couplings of the ansatz"];
   FRGAssert["Dimensions.dimensionlessAreCouplings", ListQ[dl] && SubsetQ[cps, dl], dl, "\"Dimensionless\" -> {couplings of the ansatz}"];
   FRGAssert["Dimensions.referenceDimensionless", Intersection[refs, dl] === {}, Intersection[refs, dl], "a reference cannot be dimensionless"];
   FRGAssert["Dimensions.etaNames", DuplicateFreeQ[refNames] && AllTrue[refNames, StringMatchQ[#, RegularExpression["[A-Za-z][A-Za-z0-9]*"]] &], refNames,
    "eta names must be distinct alphanumeric strings (\"nu\" gives the symbol etanu and the FlowRG name :nu)"];
   FRGAssert["Dimensions.scale", Head[kk] === Symbol && SymbolName[kk] === "k", kk, "\"Scale\" -> k (the RG scale symbol used in the regulator)"];

   basis = Prepend[refs, kk]; nb = Length[basis];
   eqs = monomialEquations[mons, cps, fields];
   unknownCps = Complement[cps, refs, dl];
   cols = Join[{1}, 1 + Range[Length[fields]], 1 + Length[fields] + Flatten[Position[cps, #] & /@ unknownCps]];
   colNames = Join[{"time"}, fields, unknownCps];
   nUnk = Length[cols];

   (* freedoms: the homogeneous system with the references as unknowns too; its solutions are the changes of units that leave
      Gamma invariant besides k; there must be exactly N+1 of them, fixed by the references *)
   mFull = equationRow[#, Join[cols, 1 + Length[fields] + Flatten[Position[cps, #] & /@ refs]], cps, fields] & /@ eqs;
   nsFull = NullSpace[mFull];
   mU = equationRow[#, cols, cps, fields] & /@ eqs;
   If[Length[nsFull] =!= nS + 1,
    (* with the references fixed, the unknowns that still move in the null space are undetermined *)
    undet = If[Length[nsFull] > nS + 1, Pick[colNames, Map[AnyTrue[#, # =!= 0 &] &, Transpose[Join[NullSpace[mU], {ConstantArray[0, nUnk]}]]], True], {}];
    FRGAssert["Dimensions.freedoms", False, <|"freedoms" -> Length[nsFull], "expected" -> nS + 1, "undetermined" -> undet|>,
     If[Length[nsFull] > nS + 1,
      "the ansatz leaves more than N+1 units free: declare the amplitudes in \"Dimensionless\" (e.g. Ak in nk Ak), or check the undetermined ones",
      "the ansatz fixes more units than expected: is each phibar dt[phi] written with coefficient 1? Is a \"Dimensionless\" coupling really dimensionless?"]]];
   AppendTo[checks, assertCheck["Dimensions.freedoms", True, Null, "", ToString[nS + 1] <> " free units besides k (Gamma dimensionless fixes the rest)"]];
   refRows = nsFull[[All, nUnk + Range[nS + 1]]];
   AppendTo[checks, assertCheck["Dimensions.referencesIndependent", MatrixRank[refRows] == nS + 1, refs,
     "the references are not independent: their dimensions are related by the ansatz, choose another set (e.g. one per field pair and one for time)",
     "references independent"]];

   (* solve each basis component: M u = -(const), const = (n_x - Dim, powers of the references); M has full column rank here
      (N+1 freedoms, all fixed by the references), so the normal equations give the unique solution if one exists *)
   consts = Join[{#MomentumOrder - dimSym}, #CouplingPowers[[Flatten[Position[cps, #] & /@ refs]]]] & /@ eqs;
   sol = LinearSolve[Transpose[mU] . mU, -Transpose[mU] . consts];
   residual = Expand[mU . sol + consts];
   If[! AllTrue[Flatten[residual], # === 0 &],
    conflicts = MapThread[If[AllTrue[#2, # === 0 &], Nothing, #1["Monomial"] <> " (" <> inputForm[#1["Term"]] <> ")"] &, {eqs, residual}];
    FRGAssert["Dimensions.consistent", False, conflicts,
     "these monomials give incompatible dimensions (a coupling appears in several monomials with different field or derivative content); check the ansatz"]];
   AppendTo[checks, assertCheck["Dimensions.consistent", True, Null, "", ToString[Length[eqs]] <> " monomial equations solved exactly (overdetermined system consistent)"]];

   vec = AssociationThread[colNames -> (Expand /@ sol)];
   timeDim = vec["time"];
   dimOf = Join[AssociationThread[fields -> (vec /@ fields)],
     AssociationThread[cps -> (Which[MemberQ[refs, #], UnitVector[nb, 1 + First[FirstPosition[refs, #]]], MemberQ[dl, #], ConstantArray[0, nb], True, vec[#]] & /@ cps)]];

   (* independent route: the action is invariant under the derived rescaling (test fields, raw density) *)
   AppendTo[checks, assertCheck["Dimensions.rescaling",
     AllTrue[Flatten[Table[rescalingCheckQ[ansatz["Density"], fields, cps, dimOf, timeDim, dimSym, d, sd],
        {d, If[dimSym === 2, {2}, {2, 3}]}, {sd, {1, 2}}]], TrueQ], Null,
     "the derived dimensions do not leave the action invariant (internal error in the dimension equations)",
     "action invariant under the derived rescaling (raw density, test fields, d = " <> If[dimSym === 2, "2", "2, 3"] <> ")"]];

   (* names of the dimensionless couplings and of the etas *)
   names = Association[# -> Lookup[nameOverride, #, If[AllTrue[dimOf[#], # === 0 &], #, Symbol[Context[#] <> "a" <> SymbolName[#]]]] & /@
      Select[cps, ! MemberQ[refs, #] &]];
   etaSyms = With[{ctx = Context @@ {First[refs]}}, Symbol[ctx <> "eta" <> #] & /@ refNames];
   aq = Symbol["Global`aq"];
   generated = Join[Cases[KeyValueMap[List, names], {c_, n_} /; c =!= n :> n], etaSyms];
   bad = Union[Select[generated, MemberQ[Join[cps, fields, {kk, aq, dimSym}], #] || OwnValues[#] =!= {} || DownValues[#] =!= {} &],
     Select[generated, Count[Join[Values[names], etaSyms], #] > 1 &]];
   FRGAssert["Dimensions.names", bad === {}, bad,
    "a generated name collides with a coupling, a field, k, aq or another generated name, or has a value: rename with \"Names\" -> <|c -> newName|>"];
   AppendTo[checks, assertCheck["Dimensions.names", True, Null, "", "dimensionless names and eta symbols collision-free"]];
   With[{undeclared = Select[Complement[unknownCps, Keys[nameOverride]], AllTrue[dimOf[#], # === 0 &] &]},
    If[undeclared =!= {}, AppendTo[warns, infoCheck["INFO", "Dimensions.derivedDimensionless",
       "dimensionless without being declared (name kept): " <> StringRiffle[ToString /@ undeclared, ", "]]]]];

   (* LO: each function's reference sets its eta through f(0) = 1, so it must be a running reference *)
   If[ansatz["Pipeline"] === "LO",
    bad = Select[ansatz["Functions"], ! MemberQ[refs, #[[2]]] || refModes[[First[FirstPosition[refs, #[[2]]]]]] === "Fixed" &];
    AppendTo[checks, assertCheck["Dimensions.functionReferences", bad === {}, bad,
      "the reference of each function (DefineAnsatz \"Functions\") must be a running reference of DeriveDimensions: f(0) = 1 fixes its eta",
      "function references are running references: " <> StringRiffle[(ToString[#[[1]]] <> " -> " <> ToString[#[[2]]]) & /@ ansatz["Functions"], ", "]]]];
   stage = stampStage[<|"Stage" -> "Dimensions", "Pipeline" -> ansatz["Pipeline"], "Upstream" -> {ansatz["Fingerprint"]},
      Sequence @@ If[ansatz["Pipeline"] === "LO", {"Functions" -> ansatz["Functions"]}, {}],
      "Ansatz" -> ansatz, "Basis" -> basis, "Scale" -> kk, "Dimension" -> dimSym,
      "References" -> refs, "EtaNames" -> AssociationThread[refs -> refNames], "Etas" -> AssociationThread[refs -> etaSyms],
      "ReferenceModes" -> AssociationThread[refs -> refModes], "ReferenceEquations" -> AssociationThread[refs -> refEqs], "Dimensionless" -> dl,
      "Vectors" -> dimOf, "Time" -> timeDim, "Frequency" -> -timeDim, "Names" -> names, "Equations" -> eqs,
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowDimensions[stage];
   PrintCard[stage];
   stage]];

(* ---------- inspection ---------- *)

dimForm[dims_, v_List] := basisForm[v, dims["Basis"]];
(* plain-text table, columns as wide as their widest cell *)
textTable[rows_] := With[{w = Max /@ Transpose[Map[StringLength, rows, {2}]]},
   StringRiffle[StringRiffle[MapThread[StringPadRight, {#, w}], " | "] & /@ rows, "\n"]];
DimensionOf[dims_Association, x_] := Block[{$FRGModule = "Dimensions"},
   validateStage[dims, "Dimensions"];
   Which[x === "Frequency" || x === "Time", dimForm[dims, dims[x]],
    KeyExistsQ[dims["Vectors"], x], dimForm[dims, dims["Vectors"][x]],
    True, FRGAssert["Dimensions.unknownSymbol", False, x, "not a coupling or field of this ansatz (or \"Frequency\", \"Time\")"]]];
LinearPart[dims_Association, c_] := Block[{$FRGModule = "Dimensions"},
   validateStage[dims, "Dimensions"];
   FRGAssert["Dimensions.unknownCoupling", MemberQ[dims["Ansatz"]["Couplings"], c], c, "not a coupling of this ansatz"];
   Expand[linearPart[dims["Vectors"][c], Values[dims["Etas"]]]]];

(* z-like exponent of a quantity of dimension v: d_s ln(value in units k^x ref^y) = x - Sum y_j eta_j *)
scalingExponent[dims_, v_] := Expand[First[v] - Rest[v] . Values[dims["Etas"]]];

ShowDimensions[dims_Association] := Module[{cps, rows, refs, header, fieldRows, mons},
   validateStage[dims, "Dimensions"];
   cps = dims["Ansatz"]["Couplings"]; refs = dims["References"];
   mons[c_] := StringRiffle[Union[Select[dims["Equations"], ! FreeQ[#Term, c] &][[All, "Monomial"]]], "; "];
   rows = {ToString[#], mons[#], displayString[dimForm[dims, dims["Vectors"][#]]],
       If[MemberQ[refs, #], "reference (eta" <> dims["EtaNames"][#] <> ")", ToString[dims["Names"][#]]],
       If[MemberQ[refs, #], "d_s " <> ToString[#] <> " = -eta" <> dims["EtaNames"][#] <> " " <> ToString[#],
        displayString[linearPart[dims["Vectors"][#], Values[dims["Etas"]]]]]} & /@ cps;
   header = {"coupling", "monomials", "dimension", "dimensionless", "linear part"};
   Print["Dimensions in the basis ", displayString[dims["Basis"]], "   (Gamma_k = Integral dt d^dx density is dimensionless, [x] = 1/k, d_s ref = -eta ref)"];
   Print[If[$FrontEnd === Null,
     textTable[Prepend[rows, header]],
     Grid[Prepend[rows, header], Frame -> All, Alignment -> Left]]];
   fieldRows = (ToString[#] <> ": " <> displayString[dimForm[dims, dims["Vectors"][#]]]) & /@ dims["Ansatz"]["Theory"]["Fields"];
   Print["fields: ", StringRiffle[fieldRows, ",  "]];
   Print["frequency [w] = ", displayString[dimForm[dims, dims["Frequency"]]], "   so  z = ", displayString[scalingExponent[dims, dims["Frequency"]]]];];

CheckDimensions[dims_Association, hand_Association] := Block[{$FRGModule = "Dimensions"},
  Module[{basis, toVec, derived, keys, mism},
   validateStage[dims, "Dimensions"]; basis = dims["Basis"];
   (* exponent vector of a product of powers of the basis symbols *)
   toVec[e_] := Module[{lg = Expand[PowerExpand[Log[e]]]},
     FRGAssert["Dimensions.handForm", FreeQ[lg /. Thread[Log /@ basis -> 0], Alternatives @@ basis] && FreeQ[lg /. Thread[Log /@ basis -> 0], Log], e,
      "a hand value must be a product of powers of " <> displayString[basis]];
     Coefficient[lg, Log[#]] & /@ basis];
   keys = Keys[hand];
   derived[x_] := Which[x === "Frequency" || x === "Time", dims[x], KeyExistsQ[dims["Vectors"], x], dims["Vectors"][x],
     True, FRGAssert["Dimensions.unknownSymbol", False, x, "not a coupling or field of this ansatz (or \"Frequency\", \"Time\")"]];
   mism = Select[keys, ! AllTrue[Expand[toVec[hand[#]] - derived[#]], # === 0 &] &];
   FRGAssert["Dimensions.handCheck", mism === {},
    <|# -> <|"hand" -> hand[#], "derived" -> dimForm[dims, derived[#]]|> & /@ mism|>,
    "hand-derived and derived dimensions differ (shown above)"];
   Print["CheckDimensions: ", Length[keys], " hand values agree exactly: ", StringRiffle[ToString /@ keys, ", "]];
   keys]];

End[];
EndPackage[];

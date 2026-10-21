(* ::Package:: *)

(* Levi-Civita operator jac in a declared d = 2 (design/LEVICIVITA.md): two scalars coupled by the cubic terms
   lk Hbar eps_ij d_i Psi d_j H and kap Hbar grad[Psi].grad[H]. The kap term creates single-wedge cross terms (wedge x scal) in the
   loop, so the F4 drop is not trivial; the jac^2 terms test the F2 reduction. Checked against the numerical oracle in d = 2. *)

jacTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> FixedDimension[2]]];
jacAnsatz = quiet[DefineAnsatz[jacTheory, "Density" ->
      Hbar (dt[H] + lk jac[Psi, H] + kap grad[Psi].grad[H] - nk lap[H]) - dk Hbar^2 - dh2 grad[Hbar].grad[Hbar] + Psibar (dt[Psi] - mk lap[Psi]) - ek Psibar^2,
    "Couplings" -> {lk, kap, nk, dk, dh2, mk, ek}]];
jacSymmetries = quiet[DeclareSymmetries[jacAnsatz,
    Neglected[_, "Reason" -> "beyond the ansatz"],
    NonRenormalised[lk, "Reason" -> "cubic coupling, not projected (test model)"],
    NonRenormalised[kap, "Reason" -> "cubic coupling, not projected (test model)"],
    "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>]];
jacVertices = quiet[BuildVertices[jacAnsatz]];
jacRegulator = quiet[DefineRegulator[jacAnsatz,
    {Hbar, H} -> nk q^2 jrh[q^2/k^2], {Hbar, Hbar} -> -2 dk jdh[q^2/k^2], {Psibar, Psi} -> mk q^2 jrp[q^2/k^2], {Psibar, Psibar} -> -2 ek jdp[q^2/k^2],
    "Shapes" -> {jrh, jdh, jrp, jdp}]];
jacPropagator = quiet[BuildPropagator[jacVertices, jacRegulator]];
jacLoops = quiet[BuildLoopIntegrands[jacVertices, jacPropagator, jacRegulator]];
jacProjected = quiet[ProjectOnMonomials[jacLoops, jacSymmetries]];
jacClosure = quiet[CheckClosure[IntegrateFrequencies[jacProjected], jacSymmetries]];
jacEquations = quiet[AssembleEquations[jacClosure]];
jacOracle = quiet[RunOracle[jacEquations]];

test["jac 2D: Epsilon key on the ansatz", jacAnsatz["Epsilon"], True];
test["jac 2D: Gamma3 Hbar H Psi has the wedge term", ! FreeQ[VertexValue[jacVertices, {Hbar, H, Psi}, {{w1, q1}, {w2, q2}}], wedge], True];
test["jac 2D: F4 drop exercised (nonzero single-wedge part in some integrand)",
  StringMatchQ[SelectFirst[jacProjected["Checks"], #Check === "Projection.wedgeOdd" &]["Detail"], ___ ~~ "dropped (F4) in " ~~ Except["0"] ~~ ___], True];
test["jac 2D: no wedge or vector left in the equations", FreeQ[jacEquations["Equations"][[All, "Integrand"]], wedge | scal], True];
(* at p = 0 the jac vertex vanishes (wedge[q, -q] = 0): lk enters through the p^2 projections only *)
test["jac 2D: lk^2 enters the Hbar Hbar p^2 equation (wedge^2 reduced by F2)", ! FreeQ[EquationOf[jacEquations, dh2]["Integrand"], lk], True];
test["jac 2D: lk enters the Hbar H p^2 equation", ! FreeQ[EquationOf[jacEquations, nk]["Integrand"], lk], True];
test["jac 2D: Oracle ran in d = 2 only", Union[jacOracle["Results"][[All, "Dimension"]]], {2}];
test["jac 2D: Oracle agrees within 10^-8 (" <> ToString[Length[jacOracle["Results"]]] <> " comparisons)",
  Length[jacOracle["Results"]] > 0 && AllTrue[jacOracle["Results"], #Error < 10^-8 &], True];

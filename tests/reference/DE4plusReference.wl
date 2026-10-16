(* ::Package:: *)

(* DE4plusReference.wl -- the shallow-water DE4plus Gamma2, Gamma3 and ansatz of
   ~/sync_cnrs/Projects/shallow_water/DE4plus_fast.wl (lines 44-46, 65-86, 98-113), copied mechanically.
   Only change: the heads are renamed so that they cannot collide with test symbols
   (Γ2 -> refGamma2, Γ3 -> refGamma3, Γ -> refV, η -> refEta, γ -> refGam, k[..] -> refK, l -> refL, h -> refH).
   The scale k inside the form* functions is unchanged. scal is FRGscalarMSR`scal (same rules as the script). *)

refGamma2[ω_, q_] = {{0, refEta[1, 1][ω, scal[q, q]], 0, refL[ω, scal[q, q]]}, {refEta[1, 1][-ω, scal[q, q]], refEta[0, 2][ω, scal[q, q]], refH[ω, scal[q, q]],
    refK[ω, scal[q, q]]}, {0, refH[-ω, scal[q, q]], 0, refGam[1, 1][ω, scal[q, q]]}, {refL[-ω, scal[q, q]], refK[-ω, scal[q, q]], refGam[1, 1][-ω, scal[q, q]],
    refGam[0, 2][ω, scal[q, q]]}};

(* Gamma3[i]: external leg i at (nu, p), row leg at (w, q), column leg at -(nu + w, p + q);
   refV[nH, nHbar, nPsi, nPsibar][mom1, mom2]: mom1, mom2 = momenta of the first two legs in field order *)
refGamma3[1][ν_, p_, ω_, q_] := {{refV[3, 0, 0, 0][ν, p, ω, q], refV[2, 1, 0, 0][ν, p, ω, q], refV[2, 0, 1, 0][ν, p, ω, q], refV[2, 0, 0, 1][ν, p, ω, q]},
  {refV[2, 1, 0, 0][ν, p, -ν - ω, -p - q], refV[1, 2, 0, 0][ν, p, ω, q], refV[1, 1, 1, 0][ν, p, ω, q], refV[1, 1, 0, 1][ν, p, ω, q]},
  {refV[2, 0, 1, 0][ν, p, -ν - ω, -p - q], refV[1, 1, 1, 0][ν, p, -ν - ω, -p - q], refV[1, 0, 2, 0][ν, p, ω, q], refV[1, 0, 1, 1][ν, p, ω, q]},
  {refV[2, 0, 0, 1][ν, p, -ν - ω, -p - q], refV[1, 1, 0, 1][ν, p, -ν - ω, -p - q], refV[1, 0, 1, 1][ν, p, -ν - ω, -p - q], refV[1, 0, 0, 2][ν, p, ω, q]}};
refGamma3[2][ν_, p_, ω_, q_] := {{refV[2, 1, 0, 0][ω, q, -ν - ω, -p - q], refV[1, 2, 0, 0][ω, q, ν, p], refV[1, 1, 1, 0][ω, q, ν, p],
   refV[1, 1, 0, 1][ω, q, ν, p]}, {refV[1, 2, 0, 0][-ν - ω, -p - q, ν, p], refV[0, 3, 0, 0][ν, p, ω, q], refV[0, 2, 1, 0][ν, p, ω, q],
   refV[0, 2, 0, 1][ν, p, ω, q]}, {refV[1, 1, 1, 0][-ν - ω, -p - q, ν, p], refV[0, 2, 1, 0][-ν - ω, -p - q, ν, p], refV[0, 1, 2, 0][ν, p, ω, q],
   refV[0, 1, 1, 1][ν, p, ω, q]}, {refV[1, 1, 0, 1][-ν - ω, -p - q, ν, p], refV[0, 2, 0, 1][-ν - ω, -p - q, ν, p], refV[0, 1, 1, 1][ν, p, -ν - ω, -p - q],
   refV[0, 1, 0, 2][ν, p, ω, q]}};
refGamma3[3][ν_, p_, ω_, q_] := {{refV[2, 0, 1, 0][ω, q, -ν - ω, -p - q], refV[1, 1, 1, 0][ω, q, -ν - ω, -p - q], refV[1, 0, 2, 0][ω, q, ν, p],
   refV[1, 0, 1, 1][ω, q, ν, p]}, {refV[1, 1, 1, 0][-ν - ω, -p - q, ω, q], refV[0, 2, 1, 0][ω, q, -ν - ω, -p - q], refV[0, 1, 2, 0][ω, q, ν, p],
   refV[0, 1, 1, 1][ω, q, ν, p]}, {refV[1, 0, 2, 0][-ν - ω, -p - q, ν, p], refV[0, 1, 2, 0][-ν - ω, -p - q, ν, p], refV[0, 0, 3, 0][ν, p, ω, q],
   refV[0, 0, 2, 1][ν, p, ω, q]}, {refV[1, 0, 1, 1][-ν - ω, -p - q, ν, p], refV[0, 1, 1, 1][-ν - ω, -p - q, ν, p], refV[0, 0, 2, 1][ν, p, -ν - ω, -p - q],
   refV[0, 0, 1, 2][ν, p, ω, q]}};
refGamma3[4][ν_, p_, ω_, q_] := {{refV[2, 0, 0, 1][ω, q, -ν - ω, -p - q], refV[1, 1, 0, 1][ω, q, -ν - ω, -p - q], refV[1, 0, 1, 1][ω, q, -ν - ω, -p - q],
   refV[1, 0, 0, 2][ω, q, -ν - ω, -p - q]}, {refV[1, 1, 0, 1][-ν - ω, -p - q, ω, q], refV[0, 2, 0, 1][ω, q, -ν - ω, -p - q],
   refV[0, 1, 1, 1][ω, q, -ν - ω, -p - q], refV[0, 1, 0, 2][ω, q, -ν - ω, -p - q]}, {refV[1, 0, 1, 1][-ν - ω, -p - q, ω, q],
   refV[0, 1, 1, 1][-ν - ω, -p - q, ω, q], refV[0, 0, 2, 1][ω, q, -ν - ω, -p - q], refV[0, 0, 1, 2][ω, q, -ν - ω, -p - q]},
  {refV[1, 0, 0, 2][-ν - ω, -p - q, ω, q], refV[0, 1, 0, 2][-ν - ω, -p - q, ω, q], refV[0, 0, 1, 2][-ν - ω, -p - q, ω, q],
   refV[0, 0, 0, 3][ω, q, -ν - ω, -p - q]}};

AnsatzBMW = {refV[1, 1, 1, 0][ν_, p_, ω_, q_] -> lk*scal[p, p + q] - lk*scal[p + q, p + q], refV[1, 1, 0, 1][ν_, p_, ω_, q_] -> 0,
   refV[1, 0, 1, 1][ν_, p_, ω_, q_] -> 0, refV[0, 1, 1, 1][ν_, p_, ω_, q_] -> 0, refV[3, 0, 0, 0][ν_, p_, ω_, q_] -> 0, refV[2, 1, 0, 0][ν_, p_, ω_, q_] -> 0,
   refV[2, 0, 1, 0][ν_, p_, ω_, q_] -> 0, refV[2, 0, 0, 1][ν_, p_, ω_, q_] -> 0, refV[1, 2, 0, 0][ν_, p_, ω_, q_] -> 0, refV[0, 3, 0, 0][ν_, p_, ω_, q_] -> 0,
   refV[0, 2, 1, 0][ν_, p_, ω_, q_] -> 0, refV[0, 2, 0, 1][ν_, p_, ω_, q_] -> 0, refV[1, 0, 2, 0][ν_, p_, ω_, q_] -> 0, refV[0, 1, 2, 0][ν_, p_, ω_, q_] -> 0,
   refV[0, 0, 3, 0][ν_, p_, ω_, q_] -> 0, refV[0, 0, 2, 1][ν_, p_, ω_, q_] -> (-lk)*scal[p, q], refV[1, 0, 0, 2][ν_, p_, ω_, q_] -> 0,
   refV[0, 1, 0, 2][ν_, p_, ω_, q_] -> 0, refV[0, 0, 1, 2][ν_, p_, ω_, q_] -> 0, refV[0, 0, 0, 3][ν_, p_, ω_, q_] -> 0, refK[ω_, s_] -> fk[s],
   refEta[1, 1][ω_, s_] -> I*ω + feta[s], refEta[0, 2][ω_, s_] -> -fDh[s], refGam[1, 1][ω_, s_] -> I*ω + fgamma[s], refGam[0, 2][ω_, s_] -> -fDphi[s],
   refL[ω_, s_] -> fl[s], refH[ω_, s_] -> -fh[s]};

formeta[s_] = nk*s*(Ak + arnh[s/k^2]);
formDh[s_] = 2*Dhk*s*(aDhk + armh[s/k^2]);
formgamma[s_] = nk*s*(ank + arnphi[s/k^2]);
formDphi[s_] = 2*Dphik*(aDphik + armphi[s/k^2] + s*D2phik);
formk[s_] = Kk*s;
forml[s_] = gk + g2k*s;
formh[s_] = Hk*lk*s;

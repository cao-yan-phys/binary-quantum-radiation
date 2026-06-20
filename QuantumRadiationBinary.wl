BeginPackage["QuantumRadiationBinary`"];

SingleGravitonCoefficient::usage = "SingleGravitonCoefficient[] gives 32/5.";
PairPowerCoefficients::usage = "PairPowerCoefficients[] gives scalar, Dirac fermion and photon coefficients.";
PairPower::usage = "PairPower[species, mu, a, Omega] returns pair-radiation power.";
AllPowerResults::usage = "AllPowerResults[mu, a, Omega, omegaMin] returns the five selected powers.";
TwoBodyProjectionNorms::usage = "TwoBodyProjectionNorms[] gives spin-2/spin-0 norms of the two-body final state.";
SpectralCoefficients::usage = "SpectralCoefficients[] gives {spin2, spin0} spectral coefficients.";
SourceSpinIntegrals::usage = "SourceSpinIntegrals[] gives the circular-source spin integrals.";
DoubleGravitonIRRawPolynomial::usage = "DoubleGravitonIRRawPolynomial[x] gives the angular polynomial for the IR pole graph.";
DoubleGravitonIRSpectrumCoefficient::usage = "DoubleGravitonIRSpectrumCoefficient[x] gives dP/dx coefficient for the IR pole graph.";
DoubleGravitonIRFiniteCoefficient::usage = "DoubleGravitonIRFiniteCoefficient[] gives the finite term.";
DoubleGravitonIRLogCoefficient::usage = "DoubleGravitonIRLogCoefficient[] gives the soft-log coefficient.";
DoubleGravitonIRPower::usage = "DoubleGravitonIRPower[mu, a, Omega, omegaMin] returns the cutoff-regulated pole-graph power.";

Begin["`Private`"];

eta = DiagonalMatrix[{1, -1, -1, -1}];
dim = 4;

MDot[u_, v_] := Expand[u . eta . v];
Lower[v_] := eta . v;
Lower2[t_] := eta . t . eta;
MTraceUU[t_] := Expand[Total[Flatten[eta t]]];
Rank4Contract[a_, b_] := Sum[
  a[[mu, nu, rho, sig]] b[[mu, nu, rho, sig]],
  {mu, dim}, {nu, dim}, {rho, dim}, {sig, dim}
];

TraceReverseLower[TUU_] := Lower2[TUU] - eta MTraceUU[TUU]/2;

eps3 = {1, I, 0}/Sqrt[2];
SourceUU[q_, amp_: 1] := Module[{e4},
  e4 = Join[{eps3 . q[[2 ;; 4]]/q[[1]]}, eps3];
  amp Outer[Times, e4, e4]
];

SpinProjectors[q_] := Module[{s = MDot[q, q], pi, p2, p0},
  pi = eta - Outer[Times, q, q]/s;
  p2 = Table[
    (pi[[mu, rho]] pi[[nu, sig]] + pi[[mu, sig]] pi[[nu, rho]])/2 -
      pi[[mu, nu]] pi[[rho, sig]]/3,
    {mu, 4}, {nu, 4}, {rho, 4}, {sig, 4}
  ];
  p0 = Table[pi[[mu, nu]] pi[[rho, sig]]/3,
    {mu, 4}, {nu, 4}, {rho, 4}, {sig, 4}];
  {p2, p0}
];

TauScalarUU[k1_, k2_] :=
  Outer[Times, k1, k2] + Outer[Times, k2, k1] - eta MDot[k1, k2];

s1 = {{0, 1}, {1, 0}};
s2 = {{0, -I}, {I, 0}};
s3 = {{1, 0}, {0, -1}};
id2 = IdentityMatrix[2];
z2 = ConstantArray[0, {2, 2}];
gamma0 = ArrayFlatten[{{id2, z2}, {z2, -id2}}];
gammaI[s_] := ArrayFlatten[{{z2, s}, {-s, z2}}];
gammaU = {gamma0, gammaI[s1], gammaI[s2], gammaI[s3]};

Slash[k_] := Sum[Lower[k][[mu]] gammaU[[mu]], {mu, 4}];
DiracAdjoint[m_] := gamma0 . ConjugateTranspose[m] . gamma0;
GammaFermionUU[k1_, k2_] := Module[{d = k1 - k2},
  Table[(gammaU[[mu]] d[[nu]] + gammaU[[nu]] d[[mu]])/4, {mu, 4}, {nu, 4}]
];
FermionSpinTensor[k1_, k2_] := Module[{g = GammaFermionUU[k1, k2]},
  FullSimplify@Table[
    Tr[Slash[k1] . g[[mu, nu]] . Slash[k2] . DiracAdjoint[g[[rho, sig]]]],
    {mu, 4}, {nu, 4}, {rho, 4}, {sig, 4}
  ]
];

Fuu[k_, eps_] := Outer[Times, k, eps] - Outer[Times, eps, k];
LorentzContract2[A_, B_] := Total[Flatten[A (eta . B . eta)]];
TauPhotonUU[k1_, eps1_, k2_, eps2_] := Module[{f1, f2, ff},
  f1 = Fuu[k1, eps1];
  f2 = Fuu[k2, eps2];
  ff = LorentzContract2[f1, f2];
  -f1 . eta . Transpose[f2] - f2 . eta . Transpose[f1] + eta ff/2
];
PhotonPhysicalK[k1_, k2_] := Module[{eps, terms},
  eps = {{0, 1, 0, 0}, {0, 0, 1, 0}};
  terms = Flatten[Table[
    With[{tau = TauPhotonUU[k1, e1, k2, e2]},
      Table[tau[[mu, nu]] tau[[rho, sig]], {mu, 4}, {nu, 4}, {rho, 4}, {sig, 4}]
    ],
    {e1, eps}, {e2, eps}], 1];
  Total[terms]
];

TwoBodyProjectionNorms[] := TwoBodyProjectionNorms[] = Module[
  {q, k1, k2, p2, p0, tauS, kS, kF, kP},
  q = {1, 0, 0, 0};
  k1 = {1/2, 0, 0, 1/2};
  k2 = {1/2, 0, 0, -1/2};
  {p2, p0} = SpinProjectors[q];
  tauS = TauScalarUU[k1, k2];
  kS = Table[tauS[[mu, nu]] tauS[[rho, sig]], {mu, 4}, {nu, 4}, {rho, 4}, {sig, 4}];
  kF = FermionSpinTensor[k1, k2];
  kP = PhotonPhysicalK[k1, k2];
  <|
    "Scalar" -> {FullSimplify[Rank4Contract[kS, p2]], FullSimplify[Rank4Contract[kS, p0]]},
    "Fermion" -> {FullSimplify[Rank4Contract[kF, p2]], FullSimplify[Rank4Contract[kF, p0]]},
    "Photon" -> {FullSimplify[Rank4Contract[kP, p2]], FullSimplify[Rank4Contract[kP, p0]]}
  |>
];

SpectralCoefficients[] := SpectralCoefficients[] = Module[{n = TwoBodyProjectionNorms[]},
  <|
    "Scalar" -> {FullSimplify[n["Scalar"][[1]]/(80 Pi)], FullSimplify[n["Scalar"][[2]]/(16 Pi)]},
    "Fermion" -> {FullSimplify[n["Fermion"][[1]]/(40 Pi)], FullSimplify[n["Fermion"][[2]]/(8 Pi)]},
    "Photon" -> {FullSimplify[n["Photon"][[1]]/(80 Pi)], FullSimplify[n["Photon"][[2]]/(16 Pi)]}
  |>
];

SourceProjectionExpressions[] := SourceProjectionExpressions[] = Module[{q, h, p2, p0, ass, q2, q0},
  q = {2, r Sqrt[1 - u^2] Cos[ph], r Sqrt[1 - u^2] Sin[ph], r u};
  h = TraceReverseLower[SourceUU[q, 1]];
  {p2, p0} = SpinProjectors[q];
  ass = 0 <= r < 2 && -1 <= u <= 1 && Element[{r, u, ph}, Reals];
  q2 = FullSimplify[
    ComplexExpand[Sum[h[[mu, nu]] Conjugate[h[[rho, sig]]] p2[[mu, nu, rho, sig]],
      {mu, 4}, {nu, 4}, {rho, 4}, {sig, 4}]],
    Assumptions -> ass
  ];
  q0 = FullSimplify[
    ComplexExpand[Sum[h[[mu, nu]] Conjugate[h[[rho, sig]]] p0[[mu, nu, rho, sig]],
      {mu, 4}, {nu, 4}, {rho, 4}, {sig, 4}]],
    Assumptions -> ass
  ];
  <|"Q2" -> q2, "Q0" -> q0|>
];

SourceSpinIntegrals[] := SourceSpinIntegrals[] = Module[{p = SourceProjectionExpressions[]},
  <|
    "I2" -> FullSimplify[2 Pi Integrate[r^2 p["Q2"], {r, 0, 2}, {u, -1, 1}]],
    "I0" -> FullSimplify[2 Pi Integrate[r^2 p["Q0"], {r, 0, 2}, {u, -1, 1}]]
  |>
];

PairPowerCoefficients[] := PairPowerCoefficients[] = Module[{c = SpectralCoefficients[], ints = SourceSpinIntegrals[], coeff},
  coeff[{c2_, c0_}] := FullSimplify[16/Pi (c2 ints["I2"] + c0 ints["I0"])];
  <|"Scalar" -> coeff[c["Scalar"]], "Fermion" -> coeff[c["Fermion"]], "Photon" -> coeff[c["Photon"]]|>
];

PairPower[species_String, mu_, a_, Omega_] := PairPowerCoefficients[][species] mu^2 a^4 Omega^8;

SingleGravitonAngularPolynomial[u_] := 1 + 6 u^2 + u^4;
SingleGravitonCoefficient[] := FullSimplify[Integrate[SingleGravitonAngularPolynomial[u], {u, -1, 1}]];
SingleGravitonPower[mu_, a_, Omega_] := SingleGravitonCoefficient[] mu^2 a^4 Omega^6;

DoubleGravitonIRRawPolynomial[x_] := Pi/75 (
  2108 - 836 (x - 1)^2 + 1124 (x - 1)^4 + 148 (x - 1)^6 + 16 (x - 1)^8
);
DoubleGravitonIRRegularPolynomial[x_] := Expand[DoubleGravitonIRRawPolynomial[x] - DoubleGravitonIRRawPolynomial[0]];
DoubleGravitonIRSpectrumCoefficient[x_] := FullSimplify[
  DoubleGravitonIRRawPolynomial[x]/(64 Pi^4 x (2 - x))
];
DoubleGravitonIREndpointCoefficient[] := FullSimplify[
  Limit[x DoubleGravitonIRSpectrumCoefficient[x], x -> 0]
];
DoubleGravitonIRLogCoefficient[] := FullSimplify[2 DoubleGravitonIREndpointCoefficient[]];
DoubleGravitonIRFiniteCoefficient[] := FullSimplify[
  Integrate[DoubleGravitonIRRegularPolynomial[x]/(64 Pi^4 x (2 - x)), {x, 0, 2}]
];
DoubleGravitonIRPower[mu_, a_, Omega_, omegaMin_] := FullSimplify[
  (DoubleGravitonIRFiniteCoefficient[] +
    DoubleGravitonIRLogCoefficient[] Log[(2 Omega - omegaMin)/omegaMin]) mu^2 a^4 Omega^8,
  Assumptions -> 0 < omegaMin < Omega && Omega > 0
];

AllPowerResults[mu_, a_, Omega_, omegaMin_] := <|
  "ScalarPair" -> PairPower["Scalar", mu, a, Omega],
  "DiracFermionPair" -> PairPower["Fermion", mu, a, Omega],
  "PhotonPair" -> PairPower["Photon", mu, a, Omega],
  "SingleGraviton" -> SingleGravitonPower[mu, a, Omega],
  "DoubleGravitonIRPoleGraph" -> DoubleGravitonIRPower[mu, a, Omega, omegaMin]
|>;

End[];
EndPackage[];

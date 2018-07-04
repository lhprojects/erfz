(* ::Package:: *)

N[Table[Erf[4 Exp[2 Pi n/ 10 I]],{n, 0, 9}],20]
N[Table[Erf[ n(1+I)],{n, 0, 9}],20]
N[Table[Erf[10 n(1+I)],{n, 0, 9}],20]

Plot[Re[Erf[x (1+I)]],{x, 0, 5}, PlotRange->Full]
Plot[
{Re[Erf[4 Exp[theta I]]],Im[Erf[4 Exp[theta I]]],Abs[Erf[4 Exp[theta I]]]},
{theta, 0, 2 Pi}
, PlotRange->Full]





#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Modeling expected change in fraction of infected hosts (modified from spatial model)*)


fbar::usage = "fbar[iP, iQ, fPU, fPI, fQU, fQI] finds the average host fecundity given the fraction of infected hosts in each patch (iP and iQ) and the fecundities of all possible hosts types (fPU, fPI, fQU, fQI, where P and Q indicate patch and U and I indicate uninfected or infected hosts).";
mbar::usage = "mbar[iP, mPU, mPI] finds the average mortality in a patch given the fraction of infected hosts (iP) and the mortality of uninfected (mPU) and infected (mPI) hosts.";
didt::usage = "didt[iP, iQ, h, v, d, c, fPU, fPI, fQI, sPU, sPI, sQU, sQI, mPU, mPI, mQU, mQI, n, spont] gives the change in the expected fraction of infected hosts in a patch. It takes as input the fraction of infected hosts in each patch (iP for the focal patch and iQ for the other patch), the horizontal and vertical transmission probabilities (h and v), the dispersal rate (d) and number of potentially infectious contacts (c), the fecundities (fPU, fPI, fQU, fQI), establishment probabilities (sPU, sPI, sQU, sQI), and mortalities (mPU, mPI, mQU, mQI) of uninfected and infected (indicated by U and I) in patches P and Q, the population size (n), and the probability of spontaneous infection (spont).";
di1dt::usage = "di1dt[i1, i2, t, environ, h, v, d, f, s, m, n] finds the expected change in the fraction of infected hosts in patch 1 over a single time step. It takes as input the fraction of infected hosts in patches 1 (i1) and 2 (i2), the environmental states of patches 1 (env1) and 2 (env2), the horizontal and vertical transmission probabilities (h and v), the dispersal rate (d), the fecundity (f) and establishment probability (s) of hosts with the harmful infection status, the mortality of hosts with the beneficial infection status (m), the population size (n), and the probability of spontaneous infection (spont).";
di2dt::usage = "di2dt[i1, i2, t, environ, h, v, d, f, s, m, n] finds the expected change in the fraction of infected hosts in patch 2 over a single time step. It takes as input the fraction of infected hosts in patches 1 (i1) and 2 (i2), the environmental states of patches 1 (env1) and 2 (env2), the horizontal and vertical transmission probabilities (h and v), the dispersal rate (d), the fecundity (f) and establishment probability (s) of hosts with the harmful infection status, the mortality of hosts with the beneficial infection status (m), the population size (n), and the probability of spontaneous infection (spont).";
guaranteedNoInfectionQ::usage = "guaranteedNoInfectionQ[h, v, d, f, s, m, n, spont] tests whether infection can be maintained when both patches are in State M. If the only possible equilibrium has all individuals UNinfected, it returns True. Otherwise, it False. It takes as input the horizontal and vertical probabilities (h and v), the dispersal rate (d), the fecundity (f) and establishment probability (s) of hosts with the harmful infection status, the mortality (m) of hists with the beneficial infection status, the population size (n), and the probability of spontaneous infection (spont).";
getStateList::usage = "getStateList[TS, offset, cycles] gives a list of the environmental states in each patch over time.\nInput: TS = time scale of environmental change; offset = fraction of TS by which the state in Patch 2 precedes Patch 1; cycles = # cycles of environmental change over which to make the list\nOutput: A list with entries of the form {tStart, tEnd, state in Patch 1, state in Patch 2}, where tStart and tEnd are the start and end points of the environmental state";
forecastInfection::usage = "forecastInfection[timescale, offset, h, v, d, f, s, m, n, spont, initialInfected, cycles, checkOutOfBounds] forecasts the fraction of infected hosts over a given number of cycles of environmental change (cycles), started from a given fraction of infected hosts (initialInfected).\n Parameters: the timescale of environmental change (timescale) , the fraction of the timescale of environmental change by which the state in Patch 2 precedes Patch 1, the horizontal and vertical transmission probabilities (h and v), the dispersal rate (d), the fecundity (f) and establishment probability (s) of hosts with the harmful infection status, the mortality of hosts with the beneficial infection status (m), the population size (n), and the probability of spontaneous infection (spont). Optional parameters: checkOutOfBounds, if True (default) returns {} if the fraction of infected hosts in either patch becomes <0 or > 1; also takes as input options to pass to NDSolve.";
findAvgInfection::usage = "findAvgInfection[timescale, offset, h, v, d, f, s, m, n, spont, initialInfected, cycles, burninCycles, checkOutOfBounds] finds the average fraction of infected hosts in each environmental state in the population as a whole as well as in each patch.\nForecasts the fraction of infected hosts over a number of cycles of environmental change (given by cycles).\nAverages the fraction of infected hosts over time starting from the end of the (burninCycles)th cycle of environmental change. To average over all cycles forecast, set burninCycles to 0.\nInput: the timescale of environmental change (timescale) , the fraction of the timescale of environmental change by which the state in Patch 2 precedes Patch 1, the horizontal and vertical transmission probabilities (h and v), the dispersal rate (d), the fecundity (f) and establishment probability (s) of hosts with the harmful infection status, the mortality of hosts with the beneficial infection status (m), the population size (n), the probability of spontaneous infection (spont), the initial fraction of infected hosts in each patch (initialInfected), the cycles of environmental change over which to forecast infection dynamics (cycles) and the cycle after which to begin calculating the average infection (burninCycles).  Optional parameters: checkOutOfBounds, if True (default) returns {} if the fraction of infected hosts in either patch becomes <0 or > 1; also takes as input options to pass to NDSolve. \nOutput: Average fraction of infected hosts in... {State M (whole population), State P (whole pop.), State M (Patch 1), State P (Patch 1), State M (Patch 2), State P (Patch 2)}";


(* ::Subsubsection::Closed:: *)
(*fbar and mbar*)


(* ::Text:: *)
(*The functions fbar and mbar finds the average host fecundity and mortality, respectively. *)
(*fbar takes as input the fraction of infected hosts in each patch (iP and iQ) as well as the fecundity of uninfected and infected hosts in each patch (fPU, fPI, fQU, fQI), where P or Q indicate the patch, and U and I indicate uninfected or infected hosts, respectively. *)
(*mbar takes as input the fraction of infected hosts in one patch (iP for a generic patch) and the mortality of uninfected and infected hosts in that patch (mPU and mPI).*)


fbar[iP_, iQ_, fPU_, fPI_, fQU_, fQI_] := (fPU * (1 - iP) + fPI * iP + fQU * (1 - iQ) + fQI * iQ)/2


mbar[iP_, mPU_, mPI_] := mPU * (1 - iP) + mPI * iP


(* ::Subsubsection::Closed:: *)
(*didt*)


(* ::Text:: *)
(*didt gives the change in infection in a generic patch in a finite population (size n). Assumes 2 patches each of size n/2. Allows spontaneous infection.*)
(*Input:*)
(*     (P indicates focal patch; Q the other patch. U and I refer to uninfected and infected hosts, respectively.)*)
(*     iP, iQ: fraction of infected hosts in each patch*)
(*     h, v: horizontal and vertical transmission probabilities, respectively (assumes entire population has same transmission probabilities)*)
(*     d: dispersal rate of newborn hosts*)
(*     c: number of potentially infectious contacts a newborn host has (assumed to be 1 for entire project)*)
(*     fPU, fPI, fQU, fQI: fecundity of uninfected and infected hosts in each patch*)
(*     sPU, sPI, sQU, sQI: establishment probabilities of newborn hosts*)
(*     mPU, mPI, mQU, mQI: adult mortality*)
(*     n: population size (used for scaling answers by number of hosts; calculations allow for fractional hosts, so setting to 1 will give answers in terms of proportions of an infinite population)*)
(*     spont: spontaneous infection probability*)
(**)
(*Output:*)
(*     Expected change in fraction of infected hosts in a patch ("P") over a single time step*)


didt[iP_, iQ_, h_, v_,d_, c_, fPU_, fPI_, fQU_, fQI_, sPU_, sPI_, sQU_, sQI_, mPU_, mPI_, mQU_, mQI_, n_, spont_]:=
(1/(n * fbar[iP, iQ, fPU, fPI, fQU, fQI] * mbar[iP, mPU, mPI])) * ( 
(* fraction infected goes up *)
(v * (1 - d) * fPI * iP     +    v * d * fQI * iQ    + (* vertical transmission *)
((1 - v) * (1 - d) * fPI * iP    +    (1 - v) * d * fQI * iQ    +    (1 - d) * fPU * (1 - iP)    +   d * fQU * (1 - iQ)) *
( 1 - (1 - spont) * (1 - h * iP)^c) (* horizontal transmission or spontaneous infection *) ) * 
sPI * (1 - iP) * mPU - 
(* fraction infected goes down *)
((1 - v) * (1 - d) * fPI * iP + (1 - v) * d * fQI * iQ + (1 - d) * fPU * (1 - iP) + d * fQU * (1 - iQ)) * 
(1 - spont) * (1 - h * iP)^c * sPU * iP * mPI )


(* ::Subsubsection::Closed:: *)
(*di1dt and di2dt*)


(* ::Text:: *)
(*di1dt and di2dt give the change in the fraction of infected hosts in patches 1 and 2, allowing spontaneous infection and assuming a fixed environmental state. These are simpler versions of didt that reflect some of the additional assumptions made in the simulations.*)
(**)
(*Input:*)
(*     i1, i2: fraction of infected hosts in patches 1 and 2*)
(*     env1, env2: environmental states in Patches 1 and 2 (0 for State M, 1 for State P)*)
(*     h, v: horizontal and vertical transmission probabilities*)
(*     d: dispersal rate of newborn hosts*)
(*     f: fecundity of hosts with harmful infection status*)
(*     s: establishment probability of hosts with harmful infection status*)
(*     m: adult mortality of hosts with beneficial infection status*)
(*     n: population size (set to 1 for infinite population)*)
(*     spont: spontaneous infection probability*)
(**)
(*Output:*)
(*     Expected change in fraction of infected hosts in patch 1 (di1dt) or 2 (di2dt) over a single time step*)
(*     *)
(*Assumptions:*)
(*     c = 1: Each newborn only has one potentially infectious contact*)
(*     Infection status is either beneficial or harmful. Whether it is beneficial or harmful is determined by the environmental state (environ).*)
(*     When the environmental state is 1, being uninfected is beneficial. When it is 0, being infected is beneficial.*)
(*     Fecundity and establishment probability are 1 if a host has a beneficial infection status, f and s otherwise*)
(*     Mortality is m if a host has a beneficial infection status, 1 otherwise*)


di1dt[i1_, i2_, env1_, env2_, h_, v_, d_, f_, s_, m_, n_, spont_] := didt[i1, i2, h, v, d, 1, 
env1 + (1 - env1) * f, env1 * f + (1 - env1), env2 + (1 - env2) * f, env2 * f + (1 - env2),
env1 + (1 - env1) * s, env1 * s + (1 - env1), env2 + (1 - env2) * s, env2 * s + (1 - env2),
env1 * m + (1 - env1), env1 + (1 - env1) * m, env2 * m + (1 - env2), env2 + (1 - env2) * m, n, spont]

di2dt[i1_, i2_, env1_, env2_, h_, v_, d_, f_, s_, m_, n_, spont_] := didt[i2, i1, h, v, d, 1, 
env2 + (1 - env2) * f, env2 * f + (1 - env2), env1 + (1 - env1) * f, env1 * f + (1 - env1),
env2 + (1 - env2) * s, env2 * s + (1 - env2), env1 + (1 - env1) * s, env1 * s + (1 - env1),
env2 * m + (1 - env2), env2 + (1 - env2) * m, env1 * m + (1 - env1), env1 + (1 - env1) * m, n, spont]


(* ::Subsubsection::Closed:: *)
(*Identifying regions where no infection is expected (guaranteedNoInfectionQ)*)


(* ::Text:: *)
(*guaranteedNoInfectionQ[h, v, d, f, s, m, n, spont] tests whether infection can be maintained when both patches are in State M.*)
(*If infection cannot be maintained in this case, infection should not be maintained in any environment, so one can assume the containment is 0.*)
(*Input:*)
(*	h, v: horizontal and vertical transmission probabilities*)
(*	d: dispersal rate of newborn hosts*)
(*	f: fecundity of hosts with harmful infection status*)
(*	s: establishment probability of hosts with harmful infection status*)
(*	m: adult mortality of hosts with beneficial infection status*)
(*	n: population size*)
(*	spont: spontaneous infection probability*)
(*Output:*)
(*	Returns TRUE if the only possible equilibrium has all individuals UNinfected*)
(*	Returns FALSE otherwise*)


guaranteedNoInfectionQ[h_, v_, d_, f_, s_, m_, n_, spont_] :=
Module[{equilList},
equilList =
	Solve[{di1dt[i1, i2, 0, 0, h, v, d, f, s, m, n, spont] == 0,
			di2dt[i1, i2, 0, 0, h, v, d, f, s, m, n, spont] == 0, 
			i1 <= 1, i1 >= 0, i2 <= 1, i2 >= 0},
		{i1, i2}];

If[Length[equilList] == 1 && equilList[[1]] == {i1 -> 0, i2 -> 0}, True, False]]


(* ::Subsubsection::Closed:: *)
(*getStateList*)


(* ::Text:: *)
(*getStateList[TS, offset, cycles] gives a list of the environmental states in each patch over time.*)
(*Input: TS = time scale of environmental change; offset = fraction of TS by which the state in Patch 2 precedes Patch 1; cycles = # cycles of environmental change over which to make the list*)
(*Output: A list with entries of the form {tStart, tEnd, state in Patch 1, state in Patch 2}, where tStart and tEnd are the start and end points of the environmental state*)
(**)
(*Possible combinations of environmental states in each patch (in order):*)
(*offset := Mod(inputOffset, 1)*)
(*For offset = 0, states are (0, 0) and (1, 1)*)
(*For offset = 1/2, states are (0, 1) and (1, 0)*)
(*For offset in (0, 1/2), states are (0, 0), (0, 1), (1, 1), (1, 0)*)
(*For offset in (1/2, 1), states are (0, 1), (0, 0), (1, 0), (1, 1)*)
(**)
(*Intervals on which the environment is in each state:*)
(*For offset = 0, state (0, 0) on [0, TS/2), state (1, 1) on [TS/2, TS)*)
(*For offset = 1/2, state (0, 1) on [0, TS/2), state (1, 0) on [TS/2, TS)*)
(*For offset in (0, 1/2), state (0, 0) on [0, TS * (0.5 - offset) ), state (0, 1) on  [TS * (0.5 - offset) , TS/2), state (1, 1) on [TS/2, TS * (1 - offset) ), state (1, 0) on [TS * (1 - offset), TS)*)
(*For offset in (1/2, 1), state (0, 1) on [0, TS * (1 - offset) ), state (0, 0) on [TS * (1 - offset), TS/2], state (1, 0) on [TS/2, TS * (1.5 - offset) ), state (1, 1) on [TS * (1.5 - offset), TS)*)


getStateList[TS_, offset_, cycles_] := (
(* Check input *)
If[Mod[cycles, 1] != 0 || cycles < 1, Print["cycles must be a whole, positive number of cycles of environmental change"]; Abort[]];
If[TS <= 0, Print["Time scale must be a positive number"]; Abort[]];
If[Abs[offset] > 1, Print["Warning: offset should be a fraction of time scale."]];

Module[{off = Mod[offset, 1], cycleCounter = 0, t0 = 0, stateList = {}},
(* Keep adding environmental states and times to the list until have reached "cycle" environmental cycles *)
While[cycleCounter < cycles,
AppendTo[stateList, #]& /@
Switch[off,
       (* Patches are completely synchronized *)
       x_/;x == 0 || x == 1,
       {{t0, t0 + TS/2, 0, 0}, {t0 + TS/2, t0 + TS, 1, 1}},
       
       (* Patches are completely anti-synchronized *)
       1/2,
       {{t0, t0 + TS/2, 0, 1}, {t0 + TS/2, t0 + TS, 1, 0}},
       
       (* Patches are intermediately synchronized; in the same state at time 0 *)
       x_/;0 < x < 1/2,
       {{t0, t0 + TS * (1/2 - off), 0, 0}, {t0 + TS * (1/2 - off), t0 + TS/2, 0, 1}, {t0 + TS/2, t0 + TS * (1 - off), 1, 1}, {t0 + TS * (1 - off), t0 + TS, 1, 0}},
       
       (* Patches are intermediately synchronized; in opposite states at time 0 *)
       x_/;1/2 < x < 1,
       {{t0, t0 + TS * (1 - off), 0, 1}, {t0 + TS * (1 - off), t0 + TS/2, 0, 0}, {t0 + TS/2, t0 + TS * (3/2 - off), 1, 0}, {t0 + TS * (3/2 - off), t0 + TS, 1, 1}}];
       
t0 += TS; (* Update the time at the start of the cycle *)
cycleCounter ++;];
Return[stateList]])


(* ::Subsubsection::Closed:: *)
(*forecastInfection*)


(* ::Text:: *)
(*Forecasts the fraction of infected hosts over a given number of cycles of environmental change (cycles), started from a given fraction of infected hosts in each patch (initialInfected). To speed things up, checks to see if infection can be maintained in the best-case scenario (both patches in State M). If not, it returns 0 for the fraction of infected hosts in both patches.*)
(*Input:*)
(*     timescale: timescale of environmental change*)
(*     offset:  fraction of the timescale of environmental change by which the state in Patch 2 precedes Patch 1*)
(*     h, v: horizontal and vertical transmission probabilities*)
(*     d: dispersal rate (d)*)
(*     f: fecundity of hosts with the harmful infection status*)
(*     s: establishment probability of hosts with the harmful infection status*)
(*     m: mortality of hosts with the beneficial infection status *)
(*     n: population size*)
(*     spont: probability of spontaneous infection*)
(*     initialInfected: initial fraction of infected hosts (same in each patch)*)
(*     cycles: cycles of environmental change over which to forecast infection*)
(*     checkOutOfBounds: if true (default), checks if the fractions of infected hosts are > 1 or < 0; returns {} if so*)
(*     Also takes options for NDSolve*)
(*Output:*)
(*     {i1[t], i2[t]}, where i1[t] and i2[t] are piecewise functions giving the fraction of infected hosts in patches 1 (i1[t]) and 2 (i2[t]) as a function of t*)


forecastInfection[timescale_, offset_, h_, v_, d_, f_, s_, m_, n_, spont_, initialInfected_, cycles_, cutoff_, checkOutOfBounds:_?BooleanQ:True, opts:OptionsPattern[]] :=
Module[{forecastList, i1piecewise, i2piecewise, stateList = getStateList[timescale, offset, cycles], i1start = initialInfected, i2start = initialInfected},
(* t must be symbolic for the part that writes the piecewise functions to work. Check that it is *)
If[ValueQ[t], Print["Please run ClearAll[t]. t must be symbolic."]; Abort[]];

(* Check if infection can be maintained in the best-case scenario (both patches in State M). 
If not, skip further calculations and output 0 for the fraction of infected hosts in both patches. *)
If[guaranteedNoInfectionQ[h, v, d, f, s, m, n, spont], Return[{0, 0}]];

forecastList =
Table[Module[{soln},
soln = NDSolve[{i1'[t] == di1dt[i1[t], i2[t], state[[3]], state[[4]], h, v, d, f, s, m, n, spont],
                i2'[t] == di2dt[i1[t], i2[t], state[[3]], state[[4]], h, v, d, f, s, m, n, spont],
                WhenEvent[i1[t] < cutoff, i1[t] -> 0, "DetectionMethod"->"Interpolation"], 
                WhenEvent[i2[t] < cutoff, i2[t] -> 0, "DetectionMethod"->"Interpolation"],
                WhenEvent[i1[t] > 1 - cutoff, i1[t] -> 1, "DetectionMethod"->"Interpolation"], 
                WhenEvent[i2[t] > 1 - cutoff, i2[t] -> 1, "DetectionMethod"->"Interpolation"],
                i1[state[[1]]] == i1start,
                i2[state[[1]]] == i2start},
               {i1, i2},
               {t, state[[1]], state[[2]]}, FilterRules[{opts}, Options[NDSolve]]];
i1start = i1[state[[2]]] /. First[soln];
i2start = i2[state[[2]]] /. First[soln];

If[checkOutOfBounds,
If[Max[Table[i1[t]/. First[soln], {t, state[[1]], state[[2]]}]] > 1 ||
Max[Table[i2[t]/. First[soln], {t, state[[1]], state[[2]]}]] > 1 ||
Min[Table[i1[t]/. First[soln], {t, state[[1]], state[[2]]}]] < 0 ||
Min[Table[i2[t]/. First[soln], {t, state[[1]], state[[2]]}]] < 0,
Print["Something's gone wrong with integration. Impossible fraction of infected hosts. Try re-running NDSolve."]; Return[{}, Table]]];
soln],
{state, stateList}];

If[forecastList == {}, Return[{}]];

i1piecewise = Piecewise[MapThread[{(i1[t]/.#1[[1]]), #2[[1]] <= t < #2[[2]]}&, {forecastList, stateList}]];
i2piecewise = Piecewise[MapThread[{(i2[t]/.#1[[1]]), #2[[1]] <= t < #2[[2]]}&, {forecastList, stateList}]];

{i1piecewise, i2piecewise}];


(* ::Subsubsection::Closed:: *)
(*findAvgInfection*)


(* ::Text:: *)
(*Finds the average fraction of infected hosts in each environmental state in the population as a whole as well as in each patch. Forecasts the fraction of infected hosts over a number of cycles of environmental change (given by cycles). Averages the fraction of infected hosts over time starting from the end of the (burninCycles)th cycle of environmental change. To average over all cycles forecast, set burninCycles to 0.*)
(*Input:*)
(*     timescale: timescale of environmental change*)
(*     offset:  fraction of the timescale of environmental change by which the state in Patch 2 precedes Patch 1*)
(*     h, v: horizontal and vertical transmission probabilities*)
(*     d: dispersal rate (d)*)
(*     f: fecundity of hosts with the harmful infection status*)
(*     s: establishment probability of hosts with the harmful infection status*)
(*     m: mortality of hosts with the beneficial infection status *)
(*     n: population size*)
(*     spont: probability of spontaneous infection*)
(*     initialInfected: initial fraction of infected hosts (same in each patch)*)
(*     cycles: cycles of environmental change over which to forecast infection*)
(*     burninCycles: the cycle after which to begin calculating the average infection (burninCycles)*)
(*     checkOutOfBounds: if true (default), checks if the fractions of infected hosts are > 1 or < 0; returns {} if so*)
(*     Also takes options for NDSolve*)
(*Output: *)
(*     Average fraction of infected hosts in... {State M (whole population), State P (whole pop.), State M (Patch 1), State P (Patch 1), State M (Patch 2), State P (Patch 2)}*)


findAvgInfection[timescale_, offset_, h_, v_, d_, f_, s_, m_, n_, spont_, initialInfected_, cycles_, burninCycles_, cutoff_, checkOutOfBounds:_?BooleanQ:True, opts:OptionsPattern[]] := 
Module[{environ, forecast, timePointsPatch1stateM, timePointsPatch1stateP, timePointsPatch2stateM, timePointsPatch2stateP,
rawInfectionPatch1stateM, rawInfectionPatch1stateP, rawInfectionPatch2stateM, rawInfectionPatch2stateP,
avgInfectionStateM, avgInfectionStateP, avgInfectionPatch1stateM, avgInfectionPatch1stateP, avgInfectionPatch2stateM, avgInfectionPatch2stateP},
environ = getStateList[timescale, offset, cycles];
forecast = forecastInfection[timescale, offset, h, v, d, f, s, m, n, spont, initialInfected, cycles, cutoff, checkOutOfBounds, opts];
If[forecast == {}, Return[{}]];

(* time points in each state *)
timePointsPatch1stateM = Select[environ, #[[3]] == 0 && #[[1]] >= burninCycles * timescale&];
timePointsPatch1stateP = Select[environ, #[[3]] == 1 && #[[1]] >= burninCycles * timescale&];
timePointsPatch2stateM = Select[environ, #[[4]] == 0 && #[[1]] >= burninCycles * timescale&];
timePointsPatch2stateP = Select[environ, #[[4]] == 1 && #[[1]] >= burninCycles * timescale&];

(* {infection, time spent in state} *)
rawInfectionPatch1stateM = Map[{NIntegrate[forecast[[1]], {t, #[[1]], #[[2]]}], #[[2]] - #[[1]]}&, timePointsPatch1stateM];
rawInfectionPatch1stateP = Map[{NIntegrate[forecast[[1]], {t, #[[1]], #[[2]]}], #[[2]] - #[[1]]}&, timePointsPatch1stateP];
rawInfectionPatch2stateM = Map[{NIntegrate[forecast[[2]], {t, #[[1]], #[[2]]}], #[[2]] - #[[1]]}&, timePointsPatch2stateM];
rawInfectionPatch2stateP = Map[{NIntegrate[forecast[[2]], {t, #[[1]], #[[2]]}], #[[2]] - #[[1]]}&, timePointsPatch2stateP];

(* Find average infection in each state and in each patch x state *)
avgInfectionStateM = Join[rawInfectionPatch1stateM, rawInfectionPatch2stateM] // Total[#[[All, 1]]]/Total[#[[All, 2]]]&;
avgInfectionStateP = Join[rawInfectionPatch1stateP, rawInfectionPatch2stateP] // Total[#[[All, 1]]]/Total[#[[All, 2]]]&;
avgInfectionPatch1stateM = rawInfectionPatch1stateM // Total[#[[All, 1]]]/Total[#[[All, 2]]]&;
avgInfectionPatch1stateP = rawInfectionPatch1stateP // Total[#[[All, 1]]]/Total[#[[All, 2]]]&;
avgInfectionPatch2stateM = rawInfectionPatch2stateM // Total[#[[All, 1]]]/Total[#[[All, 2]]]&;
avgInfectionPatch2stateP = rawInfectionPatch2stateP // Total[#[[All, 1]]]/Total[#[[All, 2]]]&;

{avgInfectionStateM, avgInfectionStateP, avgInfectionPatch1stateM, avgInfectionPatch1stateP, avgInfectionPatch2stateM, avgInfectionPatch2stateP}]


(* ::Section:: *)
(*Manually adjusting settings to find containment*)


(* ::Text:: *)
(*Find the average infection using some default settings. This will fail for some parameter values, and then we will manually adjust settings in NDSolve and forecastInfection to get those solutions.*)


justAvgs = Table[Module[{myavg = findAvgInfection[160, off, h, v, d, fsm[[1]], fsm[[2]], fsm[[3]], 1, 0, 1/2, 5, 4, 0, True, Method -> "StiffnessSwitching"]},
{off, d, fsm, h, v, myavg}], 
{off, {0, 1/8, 1/4, 3/8, 1/2}}, {d, {1/200, 1/20, 1/2}}, {fsm, {{1/2, 1, 1}, {1, 1/2, 1}, {1, 1, 1/2}}}, {h, 0, 1, 1/10}, {v, 0, 1, 1/10}];


Export["justAvgs.wdx", justAvgs]


(* ::Text:: *)
(*justAvgsPlusManual will hold the averages found using the default settings, plus those found by manually adjusting the solver/forecast settings*)


justAvgsPlusManual = justAvgs;


(* ::Subsection::Closed:: *)
(*Offset = 0, dispersal = 1/200, fecundity effects*)


(* ::Subsubsection::Closed:: *)
(*Problem: offset = 0, disp = 1/200, fecundity effects, h = 1/10, v = 6/10*)


justAvgs[[1, 1, 1, 2, 7]]


(* ::Text:: *)
(*Switching to the default method for NDSolve fixes the problem*)


justAvgsPlusManual[[1, 1, 1, 2, 7, 6]] = findAvgInfection[160, 0, 1/10, 6/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 0, True];


(* ::Subsubsection::Closed:: *)
(*Problem: offset = 0, disp = 1/200, fecundity effects, h = 1/10, v = 7/10*)


justAvgs[[1, 1, 1, 2, 8]]


(* ::Text:: *)
(*Default method for NDSolve + a cutoff of 1/10000 does the trick*)


justAvgsPlusManual[[1, 1, 1, 2, 8, 6]] = findAvgInfection[160, 0, 1/10, 7/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 1/10000, True];


(* ::Subsubsection::Closed:: *)
(*Problem: offset = 0, disp = 1/200, fecundity effects, h = 1/10, v = 8/10*)


justAvgs[[1, 1, 1, 2, 9]]


(* ::Text:: *)
(*Use default method for NDSolve + a cutoff of 1/10000*)


justAvgsPlusManual[[1, 1, 1, 2, 9, 6]] = findAvgInfection[160, 0, 1/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 1/10000, True];


(* ::Subsubsection::Closed:: *)
(*Problem: offset = 0, disp = 1/200, fecundity effects, h = 1/10, v = 9/10*)


justAvgs[[1, 1, 1, 2, 10]]


(* ::Text:: *)
(*Use default method for NDSolve*)


justAvgsPlusManual[[1, 1, 1, 2, 10, 6]] = findAvgInfection[160, 0, 1/10, 9/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 0, True];


(* ::Subsubsection::Closed:: *)
(*Problem: offset = 0, disp = 1/200, fecundity effects, h = 2/10, v = 7/10*)


justAvgs[[1, 1, 1, 3, 8]]


(* ::Text:: *)
(*Use NDSolve defaults*)


justAvgsPlusManual[[1, 1, 1, 3, 8, 6]] = findAvgInfection[160, 0, 2/10, 7/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 0, True];


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 2/10, v = 8/10*)


justAvgs[[1, 1, 1, 3, 9]]


(* ::Text:: *)
(*Use NDSolve defaults*)


justAvgsPlusManual[[1, 1, 1, 3, 9, 6]] = findAvgInfection[160, 0, 2/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 0, True];


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 2/10, v = 9/10*)


justAvgs[[1, 1, 1, 3, 10]]


(* ::Text:: *)
(*Use NDSolve defaults + cutoff = 1/10000*)


justAvgsPlusManual[[1, 1, 1, 3, 10, 6]] = findAvgInfection[160, 0, 2/10, 9/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 1/10000, True];


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 3/10, v = 6/10*)


justAvgs[[1, 1, 1, 4, 7]]


(* ::Text:: *)
(*Use NDSolve defaults*)


justAvgsPlusManual[[1, 1, 1, 4, 7, 6]] = findAvgInfection[160, 0, 3/10, 6/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 0, True];


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 3/10, v = 7/10*)


justAvgs[[1, 1, 1, 4, 8]]


(* ::Text:: *)
(*Use NDSolve defaults (had to re-run once)*)


justAvgsPlusManual[[1, 1, 1, 4, 8, 6]] = findAvgInfection[160, 0, 3/10, 7/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 0, True];


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 3/10, v = 8/10*)


justAvgs[[1, 1, 1, 4, 9]]


(* ::Text:: *)
(*NDSolve defaults produces a solution, but the fraction infected is still changing*)


findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 5, 4, 0, True]


findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 1/2, 6, 5, 0, True]


(* ::Text:: *)
(*Starting from a lower initial fraction of infected hosts produces exactly the same pattern(i.e. things aren't necessarily "trying" to decline to 0).*)


findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 1/10, 5, 4, 0, True]
findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 1/10, 6, 5, 0, True]


(* ::Text:: *)
(* Averaging over a longer time may be helpful.*)


findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 6, 4, 0, True]
findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 8, 6, 0, True]
findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 10, 4, 0, True]
findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 16, 10, 0, True]


(* ::Text:: *)
(*Because I'm a bit concerned about introducing numerical errors after running these forecasts for a while, I'll average for longer (6 cycles) starting from the same point, as well as using the default NDSolve method:*)


justAvgsPlusManual[[1, 1, 1, 4, 9, 6]] = findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 10, 4, 0, True];


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 4/10, v = 6/10*)


justAvgs[[1, 1, 1, 5, 7]]


(* ::Text:: *)
(*Changing the cutoff to 1/100000 worked.*)


justAvgsPlusManual[[1, 1, 1, 5, 7, 6]] = findAvgInfection[160, 0, 4/10, 6/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/100000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 4/10, v = 1*)


justAvgs[[1, 1, 1, 5, 11]]


(* ::Text:: *)
(*Changing the cutoff to 1/100000 worked.*)


justAvgsPlusManual[[1, 1, 1, 5, 11, 6]] = findAvgInfection[160, 0, 4/10, 10/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/100000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection:: *)
(*Problem : offset = 0, disp = 1/200, fecundity effects, h = 4/10, v = 6/10*)


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[1, 1, 1]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 0, dispersal = 1/200, mortality effects*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 0, v = 7/10*)


justAvgs[[1, 1, 3, 1, 8]]


(* ::Text:: *)
(*Changing the cutoff to 1/100000 worked*)


justAvgsPlusManual[[1, 1, 3, 1, 8, 6]] = findAvgInfection[160, 0, 0/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/100000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 0, v = 8/10*)


justAvgs[[1, 1, 3, 1, 9]]


(* ::Text:: *)
(*Changing the cutoff to 1/100000 worked*)


justAvgsPlusManual[[1, 1, 3, 1, 9, 6]] = findAvgInfection[160, 0, 0/10, 8/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/100000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 0, v = 9/10*)


justAvgs[[1, 1, 3, 1, 10]]


(* ::Text:: *)
(*Changing to default method worked*)


justAvgsPlusManual[[1, 1, 3, 1, 10, 6]] = findAvgInfection[160, 0, 0/10, 9/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 1/10, v = 6/10*)


justAvgs[[1, 1, 3, 2, 7]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 2, 7, 6]] = findAvgInfection[160, 0, 1/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 1/10, v = 7/10*)


justAvgs[[1, 1, 3, 2, 8]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 2, 8, 6]] = findAvgInfection[160, 0, 1/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 1/10, v = 8/10*)


justAvgs[[1, 1, 3, 2, 9]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 2, 9, 6]] = findAvgInfection[160, 0, 1/10, 8/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 1/10, v = 9/10*)


justAvgs[[1, 1, 3, 2, 10]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 2, 10, 6]] = findAvgInfection[160, 0, 1/10, 9/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 2/10, v = 4/10*)


justAvgs[[1, 1, 3, 3, 5]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 3, 5, 6]] = findAvgInfection[160, 0, 2/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 2/10, v = 7/10*)


justAvgs[[1, 1, 3, 3, 8]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 3, 8, 6]] = findAvgInfection[160, 0, 2/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 2/10, v = 9/10*)


justAvgs[[1, 1, 3, 3, 10]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 3, 10, 6]] = findAvgInfection[160, 0, 2/10, 9/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 2/10, v = 1*)


justAvgs[[1, 1, 3, 3, 11]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 3, 11, 6]] = findAvgInfection[160, 0, 2/10, 10/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 3/10, v = 3/10*)


justAvgs[[1, 1, 3, 4, 4]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 4, 4,6]] = findAvgInfection[160, 0, 3/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 3/10, v = 5/10*)


justAvgs[[1, 1, 3, 4, 6]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 4, 6, 6]] = findAvgInfection[160, 0, 3/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 3/10, v = 7/10*)


justAvgs[[1, 1, 3, 4, 8]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 4, 8, 6]] = 
findAvgInfection[160, 0, 3/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 3/10, v = 8/10*)


justAvgs[[1, 1, 3, 4, 9]]


(* ::Text:: *)
(*Changing to the default method worked (had to re-run once)*)


justAvgsPlusManual[[1, 1, 3, 4, 9,6]] =
findAvgInfection[160, 0, 3/10, 8/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 4/10, v = 3/10*)


justAvgs[[1, 1, 3, 5, 4]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 5, 4, 6]] = 
findAvgInfection[160, 0, 4/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 4/10, v = 4/10*)


justAvgs[[1, 1, 3, 5, 5]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 5, 5, 6]] = 
findAvgInfection[160, 0, 4/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 4/10, v = 8/10*)


justAvgs[[1, 1, 3, 5, 9]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 5, 9, 6]] =
findAvgInfection[160, 0, 4/10, 8/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 4/10, v = 1*)


(* ::Text:: *)
(*(This one looks like something weird happened with evaluating the function, not necessarily related to the parameters.)*)


justAvgs[[1, 1, 3, 5, 11]]


(* ::Text:: *)
(*Re-evaluating worked*)


justAvgsPlusManual[[1, 1, 3, 5, 11, 6]] =
findAvgInfection[160, 0, 4/10, 10/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 5/10, v = 1/10*)


justAvgs[[1, 1, 3, 6, 2]]


(* ::Text:: *)
(*Changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 1, 3, 6, 2, 6]] =
findAvgInfection[160, 0, 5/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 5/10, v = 2/10*)


justAvgs[[1, 1, 3, 6, 3]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 6, 3, 6]] =
findAvgInfection[160, 0, 5/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 5/10, v = 3/10*)


justAvgs[[1, 1, 3, 6, 4]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 6, 4, 6]] = 
findAvgInfection[160, 0, 5/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 5/10, v = 4/10*)


justAvgs[[1, 1, 3, 6, 5]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 6, 5, 6]] =
findAvgInfection[160, 0, 5/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 5/10, v = 5/10*)


justAvgs[[1, 1, 3, 6, 6]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 6, 6, 6]] =
findAvgInfection[160, 0, 5/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 5/10, v = 6/10*)


justAvgs[[1, 1, 3, 6, 7]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 6, 7, 6]] =
findAvgInfection[160, 0, 5/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 5/10, v = 9/10*)


justAvgs[[1, 1, 3, 6, 10]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 6, 10, 6]] =
findAvgInfection[160, 0, 5/10, 9/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 6/10, v = 1/10*)


justAvgs[[1, 1, 3, 7, 2]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 7, 2, 6]] =
findAvgInfection[160, 0, 6/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 6/10, v = 2/10*)


justAvgs[[1, 1, 3, 7, 3]]


(* ::Text:: *)
(*Changing to the default method worked (had to re-run once)*)


justAvgsPlusManual[[1, 1, 3, 7, 3, 6]] =
findAvgInfection[160, 0, 6/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 7/10, v = 0*)


justAvgs[[1, 1, 3, 8, 1]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 8, 1, 6]] = 
findAvgInfection[160, 0, 7/10, 0/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 7/10, v = 1/10*)


justAvgs[[1, 1, 3, 8, 2]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 8, 2, 6]] =
findAvgInfection[160, 0, 7/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 7/10, v = 2/10*)


justAvgs[[1, 1, 3, 8, 3]]


(* ::Text:: *)
(*Changing to the default method worked (had to re-run once)*)


justAvgsPlusManual[[1, 1, 3, 8, 3, 6]] =
findAvgInfection[160, 0, 7/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 7/10, v = 3/10*)


justAvgs[[1, 1, 3, 8, 4]]


(* ::Text:: *)
(*Changing to the default method worked*)


justAvgsPlusManual[[1, 1, 3, 8, 4, 6]] =
findAvgInfection[160, 0, 7/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 8/10, v = 3/10*)


justAvgs[[1, 1, 3, 9, 4]]


(* ::Text:: *)
(*Changing cutoff to 1/10000 worked (had to re-run once)*)


justAvgsPlusManual[[1, 1, 3, 9, 4, 6]] =
findAvgInfection[160, 0, 8/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 8/10, v = 4/10*)


justAvgs[[1, 1, 3, 9, 5]]


(* ::Text:: *)
(*Changing to default method worked*)


justAvgsPlusManual[[1, 1, 3, 9, 5, 6]] =
findAvgInfection[160, 0, 8/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 8/10, v = 5/10*)


justAvgs[[1, 1, 3, 9, 6]]


(* ::Text:: *)
(*Changing to default method worked*)


justAvgsPlusManual[[1, 1, 3, 9, 6, 6]] =
findAvgInfection[160, 0, 8/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 8/10, v = 6/10*)


justAvgs[[1, 1, 3, 9, 7]]


(* ::Text:: *)
(*Changing to default method worked*)


justAvgsPlusManual[[1, 1, 3, 9, 7, 6]] =
findAvgInfection[160, 0, 8/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 8/10, v = 7/10*)


justAvgs[[1, 1, 3, 9, 8]]


(* ::Text:: *)
(*Changing to default method worked*)


justAvgsPlusManual[[1, 1, 3, 9, 8, 6]] =
findAvgInfection[160, 0, 8/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 9/10, v = 3/10*)


justAvgs[[1, 1, 3, 10, 4]]


(* ::Text:: *)
(*Changing cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 1, 3, 10, 4, 6]] =
findAvgInfection[160, 0, 9/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 9/10, v = 4/10*)


justAvgs[[1, 1, 3, 10, 5]]


(* ::Text:: *)
(*Changing  to default method worked*)


justAvgsPlusManual[[1, 1, 3, 10, 5, 6]] =
findAvgInfection[160, 0, 9/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/200, mortality effects, h = 9/10, v = 5/10*)


justAvgs[[1, 1, 3, 10, 6]]


justAvgsPlusManual[[1, 1, 3, 10, 6, 6]] =
findAvgInfection[160, 0, 9/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[1, 1, 3]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 0, dispersal = 1/2, fecundity effects *** (contains important problem)*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 0, v = 6/10*)


justAvgs[[1, 3, 1, 1, 7]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 1, 1, 7, 6]] =
findAvgInfection[160, 0, 0/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 1/10, v = 5/10*)


justAvgs[[1, 3, 1, 2, 6]]


(* ::Text:: *)
(*Changing cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 2, 6, 6]] =
findAvgInfection[160, 0, 1/10, 5/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 1/10, v = 6/10*)


justAvgs[[1, 3, 1, 2, 7]]


(* ::Text:: *)
(*Changing method to BDF and changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 2, 7, 6]] = 
findAvgInfection[160, 0, 1/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, 
Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 1/10, v = 8/10*)


justAvgs[[1, 3, 1, 2, 9]]


(* ::Text:: *)
(*Changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 2, 9, 6]] = 
findAvgInfection[160, 0, 1/10, 8/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 1/10, v = 9/10*)


justAvgs[[1, 3, 1, 2, 10]]


(* ::Text:: *)
(*Using default method and changing cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 2, 10, 6]] =
findAvgInfection[160, 0, 1/10, 9/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 1/10, v = 1*)


justAvgs[[1, 3, 1, 2, 11]]


(* ::Text:: *)
(*Using default method worked*)


justAvgsPlusManual[[1, 3, 1, 2, 11, 6]] =
findAvgInfection[160, 0, 1/10, 10/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 2/10, v = 5/10*)


justAvgs[[1, 3, 1, 3, 6]]


(* ::Text:: *)
(*Changing method to BDF worked*)


justAvgsPlusManual[[1, 3, 1, 3, 6, 6]] =
findAvgInfection[160, 0, 2/10, 5/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 2/10, v = 6/10*)


justAvgs[[1, 3, 1, 3, 7]]


(* ::Text:: *)
(*Changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 3, 7, 6]] =
findAvgInfection[160, 0, 2/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 2/10, v = 7/10*)


justAvgs[[1, 3, 1, 3, 8]]


(* ::Text:: *)
(*Changing the method to BDF and the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 3, 8, 6]] =
findAvgInfection[160, 0, 2/10, 7/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 2/10, v = 8/10*)


justAvgs[[1, 3, 1, 3, 9]]


(* ::Text:: *)
(*Changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 3, 9, 6]] =
findAvgInfection[160, 0, 2/10, 8/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 2/10, v = 9/10*)


justAvgs[[1, 3, 1, 3, 10]]


(* ::Text:: *)
(*Changing method to BDF worked (had to re-run once)*)


justAvgsPlusManual[[1, 3, 1, 3, 10, 6]] =
findAvgInfection[160, 0, 2/10, 9/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 3/10, v = 4/10*)


justAvgs[[1, 3, 1, 4, 5]]


(* ::Text:: *)
(*Changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 4, 5, 6]] =
findAvgInfection[160, 0, 3/10, 4/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 3/10, v = 6/10*)


justAvgs[[1, 3, 1, 4, 7]]


(* ::Text:: *)
(*Changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 4, 7, 6]] =
findAvgInfection[160, 0, 3/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 3/10, v = 7/10*)


justAvgs[[1, 3, 1, 4, 8]]


(* ::Text:: *)
(*Changing the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 1, 4, 8, 6]] =
findAvgInfection[160, 0, 3/10, 7/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, fecundity effects, h = 3/10, v = 8/10*)


justAvgs[[1, 3, 1, 4, 9]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 1, 4, 9, 6]] =
findAvgInfection[160, 0, 3/10, 8/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Important problem : offset = 0, disp = 1/2, fecundity effects, h = 4/10, v = 6/10 *** READ THIS*)


justAvgs[[1, 3, 1, 5, 7]]


(* ::Text:: *)
(*Using the BDF method predicts a loss of infection*)


justAvgsPlusManual[[1, 3, 1, 5, 7, 6]] =
findAvgInfection[160, 0, 4/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 9/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Text:: *)
(*Points with both higher and lower vertical transmission probabilities show maintenance of infection.*)


justAvgs[[1, 3, 1, 5, 6]] (* h = 4/10, v = 5/10 *)
justAvgs[[1, 3, 1, 5, 8]] (* h = 4/10, v = 6/10 *)


(* ::Text:: *)
(*Turning off the check for out-of-bounds (<1 or <0) infection values produces an answer where infection is contained. You can see that the line just *slightly* crosses the x-axis at the end of each cycle.*)


findAvgInfection[160, 0, 4/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, False, Method -> "StiffnessSwitching"]


forecastInfection[160, 0, 4/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 1/2, 5, 0, False, Method -> "StiffnessSwitching"];

Row[{Plot[%, {t, 0, 5*160}, PlotLabel -> "Offset = 0, dispersal rate = 1/2, fecundity effects,\n h = 4/10, v = 6/10",
AxesLabel -> {"Time", "Fraction infected"}, PlotLegends -> {"Patch 1", "Patch 2"}, ImageSize -> Medium], 
Plot[%, {t, 0, 5*160}, PlotRange -> {{0, 5*160}, {-0.00001, 0.00001}}, PlotLabel -> "Same plot, zoomed in on x-axis",
AxesLabel -> {"Time", "Fraction infected"}, PlotLegends -> {"Patch 1", "Patch 2"}, ImageSize -> Medium]}]


(* ::Text:: *)
(*What happens if we try a very small amount of spontaneous infection? Infection is maintained, and the fraction of infected is different from the case where we allow the fraction infected to go negative. It's hard to tell whether this is an important difference, though, since letting the fraction of infected hosts become negative doesn't really make any sense.*)


findAvgInfection[160, 0, 4/10, 6/10, 1/2, 1/2, 1, 1, 1, 1/1000000000, 5/10, 5, 4, 0, True, Method -> "StiffnessSwitching"]
findAvgInfection[160, 0, 4/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, False, Method -> "StiffnessSwitching"]


(* ::Text:: *)
(*Try allowing spontaneous infection in some cases where we can get answers (in the no spontaneous infection case) without allowing the fraction of infected hosts to become negative. This has a pretty big effect on the fraction of infected hosts in M-states, and a small effect on the fraction of infected hosts in P-states.*)


findAvgInfection[160, 0, 4/10, 5/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "StiffnessSwitching"]
findAvgInfection[160, 0, 4/10, 5/10, 1/2, 1/2, 1, 1, 1, 1/1000000000, 5/10, 5, 4, 0, True, Method -> "StiffnessSwitching"]
findAvgInfection[160, 0, 4/10, 5/10, 1/2, 1/2, 1, 1, 1, 1/10000, 5/10, 5, 4, 0, True, Method -> "StiffnessSwitching"]
findAvgInfection[160, 0, 4/10, 5/10, 1/2, 1/2, 1, 1, 1, 1/10000, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Text:: *)
(*Questions:*)
(*1. Would allowing spontaneous infection loss make the fraction of infected hosts closer to what it is in the no spontaneous infection case? Or would it make containment larger by pushing the fraction of infected hosts lower in P-states while having little effect on M-states (the opposite of what spontaneous infection does)?*)
(*2. Does this whole thing matter when things are not completely synchronized, and the other patch would be a potential reserve of infection?*)
(*3. What looks most like the simulation results? Does that matter?*)
(**)
(*I think for right now, I'll just try to finish the plots not allowing for spontaneous infection. A priority should be looking at intermediate offsets, though.*)


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[1, 3, 1]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 0, dispersal = 1/2, mortality effects*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 0, v = 7/10*)


justAvgs[[1, 3, 3, 1, 8]]


(* ::Text:: *)
(*Using default method worked*)


justAvgsPlusManual[[1, 3, 3, 1, 8, 6]] =
findAvgInfection[160, 0, 0/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 0, v = 8/10*)


justAvgs[[1, 3, 3, 1, 9]]


(* ::Text:: *)
(*Setting cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 1, 9, 6]] = 
findAvgInfection[160, 0, 0/10, 8/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 1/10, v = 6/10*)


justAvgs[[1, 3, 3, 2, 7]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 2, 7, 6]] =
findAvgInfection[160, 0, 1/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 1/10, v = 7/10*)


justAvgs[[1, 3, 3, 2, 8]]


(* ::Text:: *)
(*Setting the cutoff to 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 2, 8, 6]] =
findAvgInfection[160, 0, 1/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 1/10, v = 8/10*)


justAvgs[[1, 3, 3, 2, 9]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 2, 9, 6]] =
findAvgInfection[160, 0, 1/10, 8/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 1/10, v = 9/10*)


justAvgs[[1, 3, 3, 2, 10]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 2, 10, 6]] = 
findAvgInfection[160, 0, 1/10, 9/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 1/10, v = 1*)


justAvgs[[1, 3, 3, 2, 11]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 2, 11, 6]] = 
findAvgInfection[160, 0, 1/10, 10/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 2/10, v = 4/10*)


justAvgs[[1, 3, 3, 3, 5]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 3, 5, 6]] = 
findAvgInfection[160, 0, 2/10, 4/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 2/10, v = 7/10*)


justAvgs[[1, 3, 3, 3, 8]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 3, 8, 6]] =
findAvgInfection[160, 0, 2/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 2/10, v = 9/10*)


justAvgs[[1, 3, 3, 3, 10]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 3, 10, 6]] =
findAvgInfection[160, 0, 2/10, 9/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 3/10, v = 7/10*)


justAvgs[[1, 3, 3, 4, 8]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 4, 8, 6]] =
findAvgInfection[160, 0, 3/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 3/10, v = 8/10*)


justAvgs[[1, 3, 3, 4, 9]]


(* ::Text:: *)
(*Using BDF for the method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 4, 9, 6]] =
findAvgInfection[160, 0, 3/10, 8/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 3/10, v = 9/10*)


justAvgs[[1, 3, 3, 4, 10]]


(* ::Text:: *)
(*Using BDF for the method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 4, 10, 6]] =
findAvgInfection[160, 0, 3/10, 9/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 4/10, v = 2/10*)


justAvgs[[1, 3, 3, 5, 3]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 5, 3, 6]] =
findAvgInfection[160, 0, 4/10, 2/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 4/10, v = 5/10*)


justAvgs[[1, 3, 3, 5, 6]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 5, 6, 6]] =
findAvgInfection[160, 0, 4/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 4/10, v = 6/10*)


justAvgs[[1, 3, 3, 5, 7]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 5, 7, 6]] =
findAvgInfection[160, 0, 4/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 4/10, v = 7/10*)


justAvgs[[1, 3, 3, 5, 8]]


justAvgsPlusManual[[1, 3, 3, 5, 8, 6]] =
findAvgInfection[160, 0, 4/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 4/10, v = 8/10*)


justAvgs[[1, 3, 3, 5, 9]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/1000 worked (had to re-run 3 times)*)


justAvgsPlusManual[[1, 3, 3, 5, 9, 6]] =
findAvgInfection[160, 0, 4/10, 8/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 4/10, v = 9/10*)


justAvgs[[1, 3, 3, 5, 10]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/1000 worked*)


justAvgsPlusManual[[1, 3, 3, 5, 10, 6]] =
findAvgInfection[160, 0, 4/10, 9/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 5/10, v = 2/10*)


justAvgs[[1, 3, 3, 6, 3]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 6, 3, 6]] =
findAvgInfection[160, 0, 5/10, 2/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 5/10, v = 3/10*)


justAvgs[[1, 3, 3, 6, 4]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 6, 4, 6]] =
findAvgInfection[160, 0, 5/10, 3/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 5/10, v = 5/10*)


justAvgs[[1, 3, 3, 6, 6]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 6, 6, 6]] =
findAvgInfection[160, 0, 5/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 5/10, v = 6/10*)


justAvgs[[1, 3, 3, 6, 7]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 6, 7, 6]] =
findAvgInfection[160, 0, 5/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 5/10, v = 7/10*)


justAvgs[[1, 3, 3, 6, 8]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 6, 8, 6]] =
findAvgInfection[160, 0, 5/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 5/10, v = 9/10*)


justAvgs[[1, 3, 3, 6, 10]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 6, 10, 6]] =
findAvgInfection[160, 0, 5/10, 9/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 6/10, v = 0*)


justAvgs[[1, 3, 3, 7, 1]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 7, 1, 6]] =
findAvgInfection[160, 0, 6/10, 0/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 6/10, v = 1/10*)


justAvgs[[1, 3, 3, 7, 2]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 7, 2, 6]] =
findAvgInfection[160, 0, 6/10, 1/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 6/10, v = 5/10*)


justAvgs[[1, 3, 3, 7, 6]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked (had to re-run once)*)


justAvgsPlusManual[[1, 3, 3, 7, 6, 6]] =
findAvgInfection[160, 0, 6/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 6/10, v = 6/10*)


justAvgs[[1, 3, 3, 7, 7]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 7, 7, 6]] =
findAvgInfection[160, 0, 6/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 6/10, v = 7/10*)


justAvgs[[1, 3, 3, 7, 8]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 7, 8, 6]] =
findAvgInfection[160, 0, 6/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 6/10, v = 8/10*)


justAvgs[[1, 3, 3, 7, 9]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 7, 9, 6]] =
findAvgInfection[160, 0, 6/10, 8/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 6/10, v = 9/10*)


justAvgs[[1, 3, 3, 7, 10]]


(* ::Text:: *)
(*Using BDF as the method worked (had to re-run once)*)


justAvgsPlusManual[[1, 3, 3, 7, 10, 6]] =
findAvgInfection[160, 0, 6/10, 9/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 7/10, v = 0*)


justAvgs[[1, 3, 3, 8, 1]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 8, 1, 6]] =
findAvgInfection[160, 0, 7/10, 0/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 7/10, v = 3/10*)


justAvgs[[1, 3, 3, 8, 4]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 8, 4, 6]] =
findAvgInfection[160, 0, 7/10, 3/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 7/10, v = 5/10*)


justAvgs[[1, 3, 3, 8, 6]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 8, 6, 6]] =
findAvgInfection[160, 0, 7/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 7/10, v = 6/10*)


justAvgs[[1, 3, 3, 8, 7]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked (had to re-run once)*)


justAvgsPlusManual[[1, 3, 3, 8, 7, 6]] =
findAvgInfection[160, 0, 7/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 7/10, v = 7/10*)


justAvgs[[1, 3, 3, 8, 8]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked (had to re-run once)*)


justAvgsPlusManual[[1, 3, 3, 8, 8, 6]] =
findAvgInfection[160, 0, 7/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 8/10, v = 2/10*)


justAvgs[[1, 3, 3, 9, 3]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 9, 3, 6]] =
findAvgInfection[160, 0, 8/10, 2/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 8/10, v = 3/10*)


justAvgs[[1, 3, 3, 9, 4]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 9, 4, 6]] =
findAvgInfection[160, 0, 8/10, 3/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 8/10, v = 5/10*)


justAvgs[[1, 3, 3, 9, 6]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 9, 6, 6]] =
findAvgInfection[160, 0, 8/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 8/10, v = 6/10*)


justAvgs[[1, 3, 3, 9, 7]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 9, 7, 6]] =
findAvgInfection[160, 0, 8/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 9/10, v = 1/10*)


justAvgs[[1, 3, 3, 10, 2]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[1, 3, 3, 10, 2, 6]] =
findAvgInfection[160, 0, 9/10, 1/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 9/10, v = 4/10*)


justAvgs[[1, 3, 3, 10, 5]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 10, 5, 6]] =
findAvgInfection[160, 0, 9/10, 4/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 9/10, v = 5/10*)


justAvgs[[1, 3, 3, 10, 6]]


(* ::Text:: *)
(*Using BDF as the method worked*)


justAvgsPlusManual[[1, 3, 3, 10, 6, 6]] =
findAvgInfection[160, 0, 9/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 0, disp = 1/2, mortality effects, h = 9/10, v = 6/10*)


justAvgs[[1, 3, 3, 10, 7]]


(* ::Text:: *)
(*Using the default method and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[1, 3, 3, 10, 7, 6]] =
findAvgInfection[160, 0, 9/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True]


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[1, 3, 3]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 1/2, dispersal = 1/200, fecundity effects*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, fecundity effects, h = 5/10, v = 1*)


justAvgs[[5, 1, 1, 6, 11]]


(* ::Text:: *)
(*Original method worked but had to re-run multiple times (somewhere between 5 and 10)*)


justAvgsPlusManual[[5, 1, 1, 6, 11, 6]] =
findAvgInfection[160, 1/2, 5/10, 10/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, fecundity effects, h = 6/10, v = 1*)


justAvgs[[5, 1, 1, 7, 11]]


(* ::Text:: *)
(*Changing the method to ImplicitRungeKutta worked*)


justAvgsPlusManual[[5, 1, 1, 7, 11, 6]] =
findAvgInfection[160, 1/2, 6/10, 10/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, fecundity effects, h = 7/10, v = 1*)


justAvgs[[5, 1, 1, 8, 11]]


(* ::Text:: *)
(*Changing the method to ImplicitRungeKutta worked*)


justAvgsPlusManual[[5, 1, 1, 8, 11, 6]] =
findAvgInfection[160, 1/2, 7/10, 10/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, fecundity effects, h = 9/10, v = 1*)


justAvgs[[5, 1, 1, 10, 11]]


(* ::Text:: *)
(*Re-running with the original parameters worked*)


justAvgsPlusManual[[5, 1, 1, 10, 11, 6]] =
findAvgInfection[160, 1/2, 9/10, 10/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, fecundity effects, h = 1, v = 5/10*)


justAvgs[[5, 1, 1, 11, 6]]


(* ::Text:: *)
(*Using Implicit Runge-Kutta as the method worked*)


justAvgsPlusManual[[5, 1, 1, 11, 6, 6]] =
findAvgInfection[160, 1/2, 10/10, 5/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, fecundity effects, h = 1, v = 1*)


justAvgs[[5, 1, 1, 11, 11]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 1, 11, 11, 6]] =
findAvgInfection[160, 1/2, 10/10, 10/10, 1/200, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[5, 1, 1]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 1/2, dispersal = 1/200, mortality effects*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 0, v = 6/10*)


justAvgs[[5, 1, 3, 1, 7]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 1, 7, 6]] =
findAvgInfection[160, 1/2, 0/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 1/10, v = 5/10*)


justAvgs[[5, 1, 3, 2, 6]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 2, 6, 6]] =
findAvgInfection[160, 1/2, 1/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 2/10, v = 4/10*)


justAvgs[[5, 1, 3, 3, 5]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 3, 5, 6]] =
findAvgInfection[160, 1/2, 2/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 3/10, v = 3/10*)


justAvgs[[5, 1, 3, 4, 4]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 4, 4, 6]] = 
findAvgInfection[160, 1/2, 3/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 4/10, v = 2/10*)


justAvgs[[5, 1, 3, 5, 3]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 5, 3, 6]] = 
findAvgInfection[160, 1/2, 4/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 5/10, v = 1/10*)


justAvgs[[5, 1, 3, 6, 2]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 6, 2, 6]] =
findAvgInfection[160, 1/2, 5/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 5/10, v = 2/10*)


justAvgs[[5, 1, 3, 6, 3]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 6, 3, 6]] =
findAvgInfection[160, 1/2, 5/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 5/10, v = 1*)


justAvgs[[5, 1, 3, 6, 11]]


(* ::Text:: *)
(*This one was really hard to find -- I couldn't get anything to work consistently. However, using the default method occasionally worked, so I just tried it over and over until it produced a result.*)


Module[{count = 0},
While[justAvgsPlusManual[[5, 1, 3, 6, 11, 6]] == {} && count < 100,
justAvgsPlusManual[[5, 1, 3, 6, 11, 6]] =
findAvgInfection[160, 1/2, 5/10, 10/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True];
count++;]]


justAvgsPlusManual[[5, 1, 3, 6, 11, 6]]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 6/10, v = 1/10*)


justAvgs[[5, 1, 3, 7, 2]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 7, 2, 6]] =
findAvgInfection[160, 1/2, 6/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 6/10, v = 1*)


justAvgs[[5, 1, 3, 7, 11]]


(* ::Text:: *)
(*Using Implicit Runge-Kutta for the method worked*)


justAvgsPlusManual[[5, 1, 3, 7, 11, 6]] =
findAvgInfection[160, 1/2, 6/10, 10/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 7/10, v = 1*)


justAvgs[[5, 1, 3, 8, 11]]


(* ::Text:: *)
(*Using Implicit Runge-Kutta for the method worked (had to re-run once)*)


justAvgsPlusManual[[5, 1, 3, 8, 11, 6]] =
findAvgInfection[160, 1/2, 7/10, 10/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 8/10, v = 1*)


justAvgs[[5, 1, 3, 9, 11]]


(* ::Text:: *)
(*Using implicit Runge-Kutta and a cutoff of 1/10000 worked*)


justAvgsPlusManual[[5, 1, 3, 9, 11, 6]] =
findAvgInfection[160, 1/2, 8/10, 10/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 9/10, v = 1*)


justAvgs[[5, 1, 3, 10, 11]]


(* ::Text:: *)
(*Using implicit Runge-Kutta for the method worked*)


justAvgsPlusManual[[5, 1, 3, 10, 11, 6]] =
findAvgInfection[160, 1/2, 9/10, 10/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/200, mortality effects, h = 1, v = 4/10*)


justAvgs[[5, 1, 3, 11, 5]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 1, 3, 11, 5, 6]] =
findAvgInfection[160, 1/2, 10/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[5, 1, 3]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 1/2, dispersal = 1/2, fecundity effects*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h =1/10, v = 6/10*)


justAvgs[[5, 3, 1, 2, 7]]


(* ::Text:: *)
(*Using BDF for the method worked*)


justAvgsPlusManual[[5, 3, 1, 2, 7, 6]] =
findAvgInfection[160, 1/2, 1/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h =1/10, v = 7/10*)


justAvgs[[5, 3, 1, 2, 8]]


(* ::Text:: *)
(*Using BDF for the method worked*)


justAvgsPlusManual[[5, 3, 1, 2, 8, 6]] =
findAvgInfection[160, 1/2, 1/10, 7/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h =1/10, v = 8/10*)


justAvgs[[5, 3, 1, 2, 9]]


(* ::Text:: *)
(*Using a cutoff of 1/10000 worked*)


justAvgsPlusManual[[5, 3, 1, 2, 9, 6]] =
findAvgInfection[160, 1/2, 1/10, 8/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 1/10000, True, Method -> "StiffnessSwitching"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h =1/10, v = 1*)


justAvgs[[5, 3, 1, 2, 11]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 3, 1, 2, 11, 6]] =
findAvgInfection[160, 1/2, 1/10, 10/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h =2/10, v = 5/10*)


justAvgs[[5, 3, 1, 3, 6]]


(* ::Text:: *)
(*Using BDF as the method worked (had to re-run once)*)


justAvgsPlusManual[[5, 3, 1, 3, 6, 6]] =
findAvgInfection[160, 1/2, 2/10, 5/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h =2/10, v = 6/10*)


justAvgs[[5, 3, 1, 3, 7]]


(* ::Text:: *)
(*Using BDF as the method worked*)


justAvgsPlusManual[[5, 3, 1, 3, 7, 6]] =
findAvgInfection[160, 1/2, 2/10, 6/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h =2/10, v = 1*)


justAvgs[[5, 3, 1, 3, 11]]


(* ::Text:: *)
(*Using BDF as the method worked*)


justAvgsPlusManual[[5, 3, 1, 3, 11, 6]] =
findAvgInfection[160, 1/2, 2/10, 10/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, fecundity effects, h = 8/10, v = 1*)


justAvgs[[5, 3, 1, 9, 11]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 3, 1, 9, 11, 6]] =
findAvgInfection[160, 1/2, 8/10, 10/10, 1/2, 1/2, 1, 1, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[5, 3, 1]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 1/2, dispersal = 1/2, mortality effects*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 0, v = 6/10*)


justAvgs[[5, 3, 3, 1, 7]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity (to allow the solution to be accurately followed near 0) worked.*)


justAvgsPlusManual[[5, 3, 3, 1, 7, 6]] =
findAvgInfection[160, 1/2, 0/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 0, v = 7/10*)


justAvgs[[5, 3, 3, 1, 8]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[5, 3, 3, 1, 8, 6]] =
findAvgInfection[160, 1/2, 0/10, 7/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 0, v = 8/10*)


justAvgs[[5, 3, 3, 1, 9]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 3, 3, 1, 9, 6]] =
findAvgInfection[160, 1/2, 0/10, 8/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 1/10, v = 5/10*)


justAvgs[[5, 3, 3, 2, 6]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[5, 3, 3, 2, 6, 6]] =
findAvgInfection[160, 1/2, 1/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 1/10, v = 6/10*)


justAvgs[[5, 3, 3, 2, 7]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[5, 3, 3, 2, 7, 6]] =
findAvgInfection[160, 1/2, 1/10, 6/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 2/10, v = 4/10*)


justAvgs[[5, 3, 3, 3, 5]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[5, 3, 3, 3, 5, 6]] =
findAvgInfection[160, 1/2, 2/10, 4/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 2/10, v = 5/10*)


justAvgs[[5, 3, 3, 3, 6]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 3, 3, 3, 6, 6]] =
findAvgInfection[160, 1/2, 2/10, 5/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 2/10, v = 1*)


justAvgs[[5, 3, 3, 3, 11]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 3, 3, 3, 11, 6]] =
findAvgInfection[160, 1/2, 2/10, 10/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 3/10, v = 3/10*)


justAvgs[[5, 3, 3, 4, 4]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[5, 3, 3, 4, 4, 6]] =
findAvgInfection[160, 1/2, 3/10, 3/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 3/10, v = 1*)


justAvgs[[5, 3, 3, 4, 11]]


(* ::Text:: *)
(*Using the default method many times (44) eventually worked*)


Module[{count = 0},
While[justAvgsPlusManual[[5, 3, 3, 4, 11, 6]] == {} && count < 100,
justAvgsPlusManual[[5, 3, 3, 4, 11, 6]] = 
findAvgInfection[160, 1/2, 3/10, 10/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True];
count++];
Print[count]]
justAvgsPlusManual[[5, 3, 3, 4, 11, 6]]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 4/10, v = 2/10*)


justAvgs[[5, 3, 3, 5, 3]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 3, 3, 5, 3, 6]] =
findAvgInfection[160, 1/2, 4/10, 2/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 4/10, v = 1*)


justAvgs[[5, 3, 3, 5, 11]]


(* ::Text:: *)
(*Using implicit Runge-Kutta for the method worked (had to re-run once)*)


justAvgsPlusManual[[5, 3, 3, 5, 11, 6]] =
findAvgInfection[160, 1/2, 4/10, 10/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 6/10, v = 1*)


justAvgs[[5, 3, 3, 7, 11]]


(* ::Text:: *)
(*Using implicit Runge-Kutta for the method worked (had to re-run once)*)


justAvgsPlusManual[[5, 3, 3, 7, 11, 6]] =
findAvgInfection[160, 1/2, 6/10, 10/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "ImplicitRungeKutta"]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/2, disp = 1/2, mortality effects, h = 1, v = 1*)


justAvgs[[5, 3, 3, 11, 11]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[5, 3, 3, 11, 11, 6]] =
findAvgInfection[160, 1/2, 10/10, 10/10, 1/2, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[5, 3, 3]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 1/4, dispersal = 1/200, fecundity effects*)


(* ::Text:: *)
(*No missing values in justAvgs for these parameters*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[3, 1, 1]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


(* ::Subsection::Closed:: *)
(*Offset = 1/4, dispersal = 1/200, mortality effects*)


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 0, v = 6/10*)


justAvgs[[3, 1, 3, 1, 7]]


(* ::Text:: *)
(*Setting the accuracy goal to something large but not infinity (10000) worked. (I got an error message when setting it to infinity sometimes, though that worked once.)*)


justAvgsPlusManual[[3, 1, 3, 1, 7, 6]] =
findAvgInfection[160, 1/4, 0/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> 10000]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 0, v = 7/10*)


justAvgs[[3, 1, 3, 1, 8]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 1, 8, 6]] =
findAvgInfection[160, 1/4, 0/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 0, v = 8/10*)


justAvgs[[3, 1, 3, 1, 9]]


(* ::Text:: *)
(*Using BDF for the method and setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 1, 9, 6]] =
findAvgInfection[160, 1/4, 0/10, 8/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, Method -> "BDF",
AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 0, v = 9/10*)


justAvgs[[3, 1, 3, 1, 10]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked (had to re-run once)*)


justAvgsPlusManual[[3, 1, 3, 1, 10, 6]] =
findAvgInfection[160, 1/4, 0/10, 9/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 1/10, v = 5/10*)


justAvgs[[3, 1, 3, 2, 6]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 2, 6, 6]] =
findAvgInfection[160, 1/4, 1/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 1/10, v = 6/10*)


justAvgs[[3, 1, 3, 2, 7]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 2, 7, 6]] =
findAvgInfection[160, 1/4, 1/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 1/10, v = 7/10*)


justAvgs[[3, 1, 3, 2, 8]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity and re-running several times worked. (Setting the precision goal to 10 also worked, but I don't totally understand what that means, so I wanted to try to get an answer with a method I understand better.)*)


justAvgsPlusManual[[3, 1, 3, 2, 8, 6]] =
findAvgInfection[160, 1/4, 1/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 1/10, v = 8/10*)


justAvgs[[3, 1, 3, 2, 9]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 2, 9, 6]] =
findAvgInfection[160, 1/4, 1/10, 8/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 1/10, v = 9/10*)


justAvgs[[3, 1, 3, 2, 10]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 2, 10, 6]] =
findAvgInfection[160, 1/4, 1/10, 9/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 2/10, v = 4/10*)


justAvgs[[3, 1, 3, 3, 5]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 3, 5, 6]] =
findAvgInfection[160, 1/4, 2/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 2/10, v = 5/10*)


justAvgs[[3, 1, 3, 3, 6]]


(* ::Text:: *)
(*Setting the accuracy goal to something large but not infinity (10000) worked*)


justAvgsPlusManual[[3, 1, 3, 3, 6, 6]] =
findAvgInfection[160, 1/4, 2/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> 10000]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 2/10, v = 6/10*)


justAvgs[[3, 1, 3, 3, 7]]


(* ::Text:: *)
(*Using the default method worked (had to re-run once)*)


justAvgsPlusManual[[3, 1, 3, 3, 7, 6]] = 
findAvgInfection[160, 1/4, 2/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 2/10, v = 7/10*)


justAvgs[[3, 1, 3, 3, 8]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity and re-running once worked*)


justAvgsPlusManual[[3, 1, 3, 3, 8, 6]] =
findAvgInfection[160, 1/4, 2/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 2/10, v = 8/10*)


justAvgs[[3, 1, 3, 3, 9]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 3, 9, 6]] =
findAvgInfection[160, 1/4, 2/10, 8/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 3/10, v = 3/10*)


justAvgs[[3, 1, 3, 4, 4]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 4, 4, 6]] =
findAvgInfection[160, 1/4, 3/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 3/10, v = 4/10*)


justAvgs[[3, 1, 3, 4, 5]]


(* ::Text:: *)
(*Using BDF as the method and setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 4, 5, 6]] =
findAvgInfection[160, 1/4, 3/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "BDF", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 3/10, v = 5/10*)


justAvgs[[3, 1, 3, 4, 6]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 4, 6, 6]] =
findAvgInfection[160, 1/4, 3/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 3/10, v = 6/10*)


justAvgs[[3, 1, 3, 4, 7]]


(* ::Text:: *)
(*Using BDF as the method and setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 4, 7, 6]] =
findAvgInfection[160, 1/4, 3/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "BDF", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 3/10, v = 7/10*)


justAvgs[[3, 1, 3, 4, 8]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 4, 8, 6]] =
findAvgInfection[160, 1/4, 3/10, 7/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True, 
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 4/10, v = 2/10*)


justAvgs[[3, 1, 3, 5, 3]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 5, 3, 6]] =
findAvgInfection[160, 1/4, 4/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 4/10, v = 3/10*)


justAvgs[[3, 1, 3, 5, 4]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 5, 4, 6]] =
findAvgInfection[160, 1/4, 4/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 4/10, v = 4/10*)


justAvgs[[3, 1, 3, 5, 5]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 5, 5, 6]] =
findAvgInfection[160, 1/4, 4/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 4/10, v = 5/10*)


justAvgs[[3, 1, 3, 5, 6]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 5, 6, 6]] =
findAvgInfection[160, 1/4, 4/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method ->"StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 4/10, v = 6/10*)


justAvgs[[3, 1, 3, 5, 7]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 5, 7, 6]] =
findAvgInfection[160, 1/4, 4/10, 6/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method ->"StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 5/10, v = 1/10*)


justAvgs[[3, 1, 3, 6, 2]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 6, 2, 6]] =
findAvgInfection[160, 1/4, 5/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 5/10, v = 2/10*)


justAvgs[[3, 1, 3, 6, 3]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 6, 3, 6]] =
findAvgInfection[160, 1/4, 5/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 5/10, v = 3/10*)


justAvgs[[3, 1, 3, 6, 4]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 6, 4, 6]] =
findAvgInfection[160, 1/4, 5/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 5/10, v = 4/10*)


justAvgs[[3, 1, 3, 6, 5]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 6, 5, 6]] =
findAvgInfection[160, 1/4, 5/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 5/10, v = 5/10*)


justAvgs[[3, 1, 3, 6, 6]]


(* ::Text:: *)
(*Using BDF for the method and setting the accuracy to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 6, 6, 6]] =
findAvgInfection[160, 1/4, 5/10, 5/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "BDF", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 6/10, v = 1/10*)


justAvgs[[3, 1, 3, 7, 2]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 7, 2, 6]] =
findAvgInfection[160, 1/4, 6/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 6/10, v = 2/10*)


justAvgs[[3, 1, 3, 7, 3]]


(* ::Text:: *)
(*Using BDF and setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 7, 3, 6]] =
findAvgInfection[160, 1/4, 6/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "BDF", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 6/10, v = 3/10*)


justAvgs[[3, 1, 3, 7, 4]]


(* ::Text:: *)
(*Setting the accuracy goal to something large but not infinity (10000) worked*)


justAvgsPlusManual[[3, 1, 3, 7, 4, 6]] =
findAvgInfection[160, 1/4, 6/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> 10000]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 6/10, v = 4/10*)


justAvgs[[3, 1, 3, 7, 5]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 7, 5, 6]] =
findAvgInfection[160, 1/4, 6/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 7/10, v = 0*)


justAvgs[[3, 1, 3, 8, 1]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 8, 1, 6]] =
findAvgInfection[160, 1/4, 7/10, 0/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 7/10, v = 1/10*)


justAvgs[[3, 1, 3, 8, 2]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked (had to re-run once)*)


justAvgsPlusManual[[3, 1, 3, 8, 2, 6]] =
findAvgInfection[160, 1/4, 7/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 7/10, v = 2/10*)


justAvgs[[3, 1, 3, 8, 3]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 8, 3, 6]] =
findAvgInfection[160, 1/4, 7/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 7/10, v = 3/10*)


justAvgs[[3, 1, 3, 8, 4]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 8, 4, 6]] =
findAvgInfection[160, 1/4, 7/10, 3/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 8/10, v = 0*)


justAvgs[[3, 1, 3, 9, 1]]


(* ::Text:: *)
(*Using the default method worked*)


justAvgsPlusManual[[3, 1, 3, 9, 1, 6]] =
findAvgInfection[160, 1/4, 8/10, 0/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 8/10, v = 1/10*)


justAvgs[[3, 1, 3, 9, 2]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked (had to re-run once)*)


justAvgsPlusManual[[3, 1, 3, 9, 2, 6]] =
findAvgInfection[160, 1/4, 8/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 8/10, v = 2/10*)


justAvgs[[3, 1, 3, 9, 3]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked (had to re-run once)*)


justAvgsPlusManual[[3, 1, 3, 9, 3, 6]] =
findAvgInfection[160, 1/4, 8/10, 2/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 8/10, v = 4/10*)


justAvgs[[3, 1, 3, 9, 5]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 9, 5, 6]] =
findAvgInfection[160, 1/4, 8/10, 4/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection::Closed:: *)
(*Problem : offset = 1/4, disp = 1/200, mortality effects, h = 9/10, v = 1/10*)


justAvgs[[3, 1, 3, 10, 2]]


(* ::Text:: *)
(*Setting the accuracy goal to infinity worked*)


justAvgsPlusManual[[3, 1, 3, 10, 2, 6]] =
findAvgInfection[160, 1/4, 9/10, 1/10, 1/200, 1, 1, 1/2, 1, 0, 5/10, 5, 4, 0, True,
Method -> "StiffnessSwitching", AccuracyGoal -> Infinity]


(* ::Subsubsection:: *)
(*Plot with problems resolved*)


ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[3, 1, 3]], {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False, PlotLegends -> Automatic,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
FrameLabel -> {{"Vertical transmission probability", None}, {"Horizontal transmission probability", None}},
PlotLabel -> "Containment", LabelStyle -> {FontSize -> 14}]


justAvgsPlusManual[[3, 1, 3, 11]] //MatrixForm


(* ::Subsection::Closed:: *)
(*All plots*)


numContPlot[avgs_] :=
ListDensityPlot[Flatten[Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, avgs, {2}], 1],
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, ColorFunctionScaling -> False,
ColorFunction -> (If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&),
Frame -> True, FrameTicks -> {{{0, 0.25, 0.5, 0.75, 1.0}, None}, {{0, 0.25, 0.5, 0.75, 1.0}, None}}, AspectRatio -> 1,
FrameLabel -> {{"Pr(VT)", None}, {"Pr(HT)", None}},
LabelStyle -> {Black, FontSize -> 12, FontFamily -> "Times"},
FrameTicksStyle -> {Black, FontSize -> 12, FontFamily -> "Times"}]


?Grid


Grid[{
{""                                      , ""         , "Dispersal rate ="                        , SpanFromLeft                               , ""                       },
{""                                      , ""         , "0.005"                                   , "0.5"                                      , ""                       },
{Rotate["Symbiont affects...", 90 Degree], "lifespan" , numContPlot[justAvgsPlusManual[[1, 1, 3]]], numContPlot[justAvgsPlusManual[[1, 3, 3]]], "Complete synchronicity" },
{SpanFromAbove                           , "fecundity", numContPlot[justAvgsPlusManual[[1, 1, 1]]], numContPlot[justAvgsPlusManual[[1, 3, 1]]], SpanFromAbove            },
{SpanFromAbove                           , ""         , ""                                        , ""                                         , ""                       },
{SpanFromAbove                           , "lifespan" , numContPlot[justAvgsPlusManual[[5, 1, 3]]], numContPlot[justAvgsPlusManual[[5, 3, 3]]], "Complete asynchronicity"},
{SpanFromAbove                           , "fecundity", numContPlot[justAvgsPlusManual[[5, 1, 1]]], numContPlot[justAvgsPlusManual[[5, 3, 1]]], SpanFromAbove            },
{""                                      ,
            BarLegend[{(If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&), {0, 1}}, LegendLayout -> "Row", LegendLabel -> Placed["Containment", Bottom], LabelStyle -> {Black, FontSize -> 12, FontFamily -> "Times"}], 
					                                                                                SpanFromLeft                              , SpanFromLeft                               , ""                       }},
Frame -> {None, None, {{{3, 4}, {3, 5}} -> True, {{6, 7}, {3, 5}} -> True}},
ItemStyle -> {Black, FontSize -> 12, FontFamily -> "Times", TextAlignment -> Center},
Alignment -> {Center, Center}]



Export["numerical_containment_2020_01_06.pdf", %87]


(* ::Text:: *)
(*Compare with simulation containment values*)


intermedTSnoEvoContainmentData = Import["C:/Users/Alexandra/Desktop/host_control_2019_02_06_summary_statistics.csv"];


gatheredNoEvoContainmentByTM = GatherBy[intermedTSnoEvoContainmentData[[2;;All, {2, 3, 4, 5, 6, 30, 35}]], {#[[1]]&, #[[2]]&, #[[3]]&, #[[4]]&, #[[5]]&}];


(* Color scheme for symbiont containment when it is averaged across replicate simulations. Negative containment (where symbiont is more prevalent where it is parasitic) is shown in black. *)
averageContainmentColors[c_] := (*If[Round[c, 0.06] < 0, ColorData["GrayTones"][0.5+0.5*c], 
If[c <0, ColorData["Candy"][11/20], ColorData["Candy"][(11/20) + (9c/20)]]];*)
(*If[c < 0, If[c >= -0.02, ColorData["Candy"][3/20], ColorData["Candy"][4/20]], ColorData["Candy"][(4/20) + (16c/20)]];*)
If[c < 0, ColorData["Candy"][4/20], ColorData["Candy"][(4/20) + (16 c/20)]];


(* Loop through symbiont effects; will generate a grid of plots for each symbiont effect *)
containmentPlots = Table[
Module[{containment, averageContainment}, 
containment = Map[{#[[4]], #[[5]], #[[7]] - #[[6]]}&, gatheredNoEvoContainmentByTM[[offset, fx, disp]], {3}];
averageContainment = Map[Mean, containment, {2}];
ListDensityPlot[Flatten[averageContainment, 1], ColorFunction -> averageContainmentColors, ColorFunctionScaling -> False, 
PlotRange -> {{0, 1}, {0, 1}, {-1, 1}}, AspectRatio -> 1, FrameTicks -> {{{0, 0.25, 0.5, 0.75, 1}, None}, {{0, 0.25, 0.5, 0.75, 1}, None}}, LabelStyle -> {FontSize -> 12, FontFamily -> "Times", FontColor -> Black}]], 
{offset, {5, 2, 4}}, {disp, {1, 3, 5}}, {fx, {1, 3}}];


gatheredNoEvoContainmentByTM[[2, 3, 1, 1, 1]]


Grid[{
{""                                      , ""         , "Dispersal rate ="          , SpanFromLeft               , ""                       },
{""                                      , ""         , "0.005"                     , "0.5"                      , ""                       },
{Rotate["Symbiont affects...", 90 Degree], "lifespan" , containmentPlots[[1, 1, 2]], containmentPlots[[1, 3, 2]], "Complete synchronicity" },
{SpanFromAbove                           , "fecundity", containmentPlots[[1, 1, 1]], containmentPlots[[1, 3, 1]], SpanFromAbove            },
{SpanFromAbove                           , ""         , ""                         , ""                          , ""                       },
{SpanFromAbove                           , "lifespan" , containmentPlots[[3, 1, 2]], containmentPlots[[3, 3, 2]], "Complete asynchronicity"},
{SpanFromAbove                           , "fecundity", containmentPlots[[3, 1, 1]], containmentPlots[[3, 3, 1]], SpanFromAbove            },
{""                                      ,
            BarLegend[{(averageContainmentColors[#]&), {0, 1}}, LegendLayout -> "Row", LegendLabel -> Placed["Containment", Bottom], LabelStyle -> {Black, FontSize -> 12, FontFamily -> "Times"}], 
					                                                                 SpanFromLeft                , SpanFromLeft                               , ""                       }},
Frame -> {None, None, {{{3, 4}, {3, 5}} -> True, {{6, 7}, {3, 5}} -> True}},
ItemStyle -> {Black, FontSize -> 12, FontFamily -> "Times", TextAlignment -> Center},
Alignment -> {Center, Center}]



Export["simulation_containment_2019_02_06.pdf",%]


Map[{#[[4]], #[[5]], #[[6, 1]] - #[[6, 2]]}&, justAvgsPlusManual[[5, 1, 1]], {2}] //MatrixForm


Module[{containment, averageContainment}, 
containment = Map[{#[[4]], #[[5]], #[[7]] - #[[6]]}&, gatheredNoEvoContainmentByTM[[4, 1, 1]], {3}];
averageContainment = Map[Mean, containment, {2}]] //MatrixForm


(* ::Text:: *)
(*Intermediate containment, low dispersal, both analytical and simulation results*)


Grid[{
{""                                      , ""         , "Analytical model"                                                  , "Simulations"               , ""},
{Rotate["Symbiont affects...", 90 Degree], "lifespan" , Show[numContPlot[justAvgsPlusManual[[3, 1, 3]]], FrameLabel -> None], containmentPlots[[2, 1, 2]], Rotate["Pr(VerticalTransmission)", 90 Degree]},
{SpanFromAbove                           , "fecundity", Show[numContPlot[justAvgsPlusManual[[3, 1, 1]]], FrameLabel -> None], containmentPlots[[2, 1, 1]], SpanFromAbove},
{"", "Pr(Horizontal transmission)", SpanFromLeft},
{"\nAnalytical vs. simulation model. Low dispersal (d = 0.005).\nIntermediate synchronicity (patches in same state 50% of the time).", SpanFromLeft, SpanFromLeft, SpanFromLeft, SpanFromLeft},
{"", BarLegend[{(If[# < -1/1000000, Green, ColorData["Candy"][(4/20) + (16 #/20)]]&), {0, 1}}, LegendLayout -> "Row", LegendLabel -> Placed["Containment", Bottom], LabelStyle -> {Black, FontSize -> 12, FontFamily -> "Times"}], 
					                      SpanFromLeft                                             , SpanFromLeft}},
ItemStyle -> {Black, FontSize -> 12, FontFamily -> "Times", TextAlignment -> Center},
Alignment -> {Center, Center}]


Export["2020_01_22_TH_numerical_vs_sim_intermediate_synchronicity.pdf", %]


(* ::Subsection::Closed:: *)
(*Export results*)


Export["justAvgsPlusManual.wdx", justAvgsPlusManual];

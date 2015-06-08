(* Mathematica Package         *)
(* Created by IntelliJ IDEA    *)

(* :Title: sumHillsFofT     *)
(* :Context: sumHillsFofT`  *)
(* :Author: Tommy            *)
(* :Date: 6/8/15              *)

(* :Package Version: 1.0       *)
(* :Mathematica Version:       *)
(* :Copyright: (c) 2015 Tommy *)
(* :Keywords:                  *)
(* :Discussion:                *)

BeginPackage["sumHillsFofT`"]
(* Exported symbols added here with SymbolName::usage *)

sumHills::usage = "sumHills[HILLS_file, options] returns a list of 2D arrays that
  are the time steps of the growth of the height of the negative of the free energy
  surface from a PLUMED metadynamics calculation.
  Only made for 2 collective variables currently, but that can't be changed by
  making it check for the length of a row in the input data."
(* todo if we ever do well tempered metadynamics, will need to change this *)

plotHills::usage = "plotHills[list of matrices, options] Takes output of sumHills
  and plots time steps."

plotHillsPoint::usage = "plotHllsPoint[list of matrices, {x, y}, options] takes output of
  sumHills and plots the selected point as a function of time."

Begin["`Private`"]  Begin Private Context

Options[sumHills] =
    {
      GridSize -> 0.1
    }

sumHills[hillsFileName_, OptionsPattern[]]:=
  Module[
    {
      gridsize,
      rawdata,
      sigmaCV1, sigmaCV2,
      minMaxCV1, minMaxCV2,
      gridLengthCV1, gridLengthCV2,
      gaussianMatrix,
      scaledRotatedGaussMat,
    },
    (* Import data, first two lines are comments *)
    rawdata = Import[hillsFileName, "Table"][[3;;]];
    (* Check for empty element at end from \n presumably*)
    rawdata = If[rawdata[[-1]] == {},
      rawdata[[1;;-2]],
      rawdata];
    sigmaCV1 = rawdata[[1,4]];
    sigmaCV2 = rawdata[[1,5]];
    minMaxCV1 = MinMax[rawdata[[All, 2]]];
    minMaxCV2 = MinMax[rawdata[[All, 3]]];
    gridSize = OptionValue[GridSize];
    (* Find size (dimensions) of grid needed. *)
    gridLengthCV1 = Ceiling[(minMaxCV1[[2]] - minMaxCV1[[1]]) / gridSize];
    gridLengthCV2 = Ceiling[(minMaxCV2[[2]] - minMaxCV2[[1]]) / gridSize];
    gaussianMatrix = GaussianMatrix[
      {{gridLengthCV1, gridLengthCV2},
        {sigmaCV1 / gridSize, sigmaCV2 / gridSize}},
      Standardized -> False,
      Method -> "Gaussian"]
        * 2 Pi (sigmaCV1 * sigmaCV2) / gridSize^2;
    scaledRotatedGaussMat[row_] := Return[
      RotateLeft[
        gaussianMatrix * row[[6]],
        Round[
          {gridLengthCV1 - (row[[2]] - minMaxCV1[[1]])/gridSize,
          gridLengthCV2 - (row[[3]] - minMaxCV2[[1]])/gridSize}
        ]][[1 ;; gridLengthCV1, 1 ;; gridLengthCV2]]
      ];
    timedGaussians = ParallelMap[scaledRotatedGaussMat,rawdata];
    accumulatedGaussians = Accumulate[timedGaussians]
  ]



End[]  End Private Context

EndPackage[]
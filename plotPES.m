(* Mathematica Package *)

BeginPackage["plotPES`"]
(* Exported symbols added here with SymbolName::usage *)  

plotPES::usage = 
	"plotPES[file name], Options: export->False, convToeV->True"

plotArrayPES::usage = 
	"plotPES[file name], Options: export->False, gausTDDFT->False, 
	convToeV->True"

Begin["`Private`"] (* Begin Private Context *) 

debugPrint[args___]:=$debugPrintFunction[args]

$debugPrintFunction=Null&;
(*$debugPrintFunction=Print[##]&*)
(*Uncomment the line above to get more (maybe) useful output during evaluation*)

importPES[fileName_] := 
 Flatten[Import[ToString[fileName], "Data"]]
 
importArrayPES[fileName_] := 
 Module[{imported}, 
 	debugPrint["importing array"];
  imported = 
   ImportString[
    StringReplace[
     Import[ToString[fileName], "String"], {"[" -> "", "]" -> "", 
      "'" -> ""}], "Data"];
  (*Check for empty element at end and remove if necessary*)
  imported = 
   If[imported[[-1]] == {""}, imported[[1 ;; -2]], imported]]
   
makePlotPES[data_, opts:OptionsPattern[{makePlotPES,ListLinePlot}]]:=
	Module[{numPoints, plotdata, range, rangeSize, graph, midpoint},
		debugPrint["Attempting to plot with opts:"];
		debugPrint[opts];
		If[OptionValue[array],
			numPoints = Length[data[[All, 1]]];
			plotdata = Transpose[data],
			numPoints = Length[data];
			plotdata = data
		];
		midpoint = Ceiling[numPoints / 2];
		rangeSize = (Min[data] - Max[data])/10;
		range = {Min[data] + rangeSize, Max[data] - rangeSize};
		debugPrint["making plot with options"];
		debugPrint[FilterRules[{opts}, Options[ListLinePlot]]];
		graph = Show[
			ListLinePlot[plotdata, 
				FilterRules[{opts}, Options[ListLinePlot]], 
				PlotRange -> {{1, numPoints}, range}
			], 
			Graphics[Line[Transpose[{{midpoint, midpoint}, range}]]]];
		debugPrint["made Plot"];
		If[StringQ[OptionValue[export]],
			Export[OptionValue[export] <> ".eps", graph],
			""
		];
		graph
	]

convTDDFT[data_, opts:OptionsPattern[]]:=
	Module[{eVGSEnergies, zeroedeVGSEnergies,tempdata},
		debugPrint["Converting Gaussian TDDFT output"];
		tempdata = data;
		eVGSEnergies = tempdata[[All, 1]] * 27.21139;
		zeroedeVGSEnergies = eVGSEnergies - Min[eVGSEnergies];
		tempdata[[All, 1]] = 0;
		tempdata = tempdata + zeroedeVGSEnergies
	]

convertToeV[data_, opts:OptionsPattern[]]:=
	Module[{tempdata},
		(*Need to use some temp name here because "data" is just the array given,
		not the name of the array so you can't assign to it.*)
		debugPrint["Converting data to eV"];
		debugPrint[data];
		tempdata = data;
		tempdata = tempdata * 27.21139;
		debugPrint["Zeroing data to lowest point"];
		tempdata = tempdata - Min[tempdata];
		debugPrint["Data converted to eV"];
		tempdata
	]

plotPES[fileName_, opts:OptionsPattern[
			{plotPES, ListLinePlot, convertToeV, convTDDFT, makePlotPES}
		]]:=
	Module[{data},
		data = importPES[fileName];
		If[OptionValue[convToeV],
			data = convertToeV[data, FilterRules[{opts}, Options[converToeV]]],
			""
		];
		debugPrint["Did data conversions, sending to makePlot"];
		makePlotPES[data, 
			Sequence @@ FilterRules[{opts} ~Join~ Options[plotPES],
				{Options[makePlotPES], Options[ListLinePlot]}
			]
		]
	]

plotArrayPES[fileName_String, opts:OptionsPattern[
			{plotArrayPES, ListLinePlot, convertToeV, convTDDFT, makePlotPES}
		]]:=
	Module[{data},
		data = importArrayPES[fileName];
		If[OptionValue[gausTDDFT],
			data = convTDDFT[data, FilterRules[{opts}, Options[convTDDFT]]],
			If[OptionValue[convToeV],
				data = convertToeV[data, FilterRules[{opts}, Options[convertToeV]]]
			]
		];
		debugPrint["Did data conversions, sending to makePlot with opts"];
		debugPrint[FilterRules[{opts}, 
				{Options[makePlotPES], Options[ListLinePlot]}
			]];
		debugPrint["All opts:", opts];
		makePlotPES[data, array->True,
			Sequence @@ FilterRules[{opts} ~Join~ Options[plotArrayPES],
				{Options[makePlotPES], Options[ListLinePlot]}
			]
		]
	]

Options[makePlotPES] = {export->False, array->False}
Options[plotArrayPES] = 
	{
		gausTDDFT->False, 
		convToeV->True,
		PlotStyle -> Black, 
		LabelStyle -> Black, 
		InterpolationOrder -> 2, 
		Frame -> True, 
		FrameTicks -> {{True, False}, {Array[{#,""}&,100], False}}
	}
Options[plotPES] = 
	{
		convToeV->True,
		PlotStyle -> Black, 
		LabelStyle -> Black,
		InterpolationOrder -> 2, 
		Frame -> True, 
		FrameTicks -> {{True, False}, {Array[{#,""}&,100], False}}
	}



End[] (* End Private Context *)

EndPackage[]
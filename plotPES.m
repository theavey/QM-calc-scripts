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

Options[plotArrayPES] =
    {
      gausTDDFT -> False,
      convToeV -> True,
      PlotStyle -> Black,
      LabelStyle -> Black,
      InterpolationOrder -> 2,
      Frame -> True,
      FrameTicks -> {{True, False}, {Array[{#,""}&,100], False}}
    }

Options[makePlotPES] =
(*I think options may need to be declared before the function is? For OptionsPattern to work?*)
    {
      export -> False,
      dimensions -> 1,
      array -> False,
      ignoreFirst -> 0,
      ignoreLast -> 0
    }

Options[plotPES] =
    {
      convToeV -> True,
      PlotStyle -> Black,
      LabelStyle -> Black,
      InterpolationOrder -> 2,
      Frame -> True,
      FrameTicks -> {{True, False}, {Array[{#,""}&,100], False}}
    }

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

   
makePlotPES[data_, opts:OptionsPattern[{makePlotPES, ListLinePlot, ListPlot3D}]]:=
	Module[{numPoints, plotdata, range, rangeSize, graph, midpoint, dims},
		debugPrint["Attempting to plot with opts:"];
		debugPrint[opts];
		If[OptionValue[array],
			numPoints = Length[data[[All, 1]]];
			plotdata = Transpose[data],
			numPoints = Length[data];
			plotdata = data
		];
    dims = OptionValue[dimensions];
    If[dims > 1,
      numPoints = Length[data]^(1 / dims);
      plotdata = Partition[#, numPoints] & /@ plotdata;
      If[dims > 2, Print["I don't know how to handle dimensions > 2"]]
    ];
    plotdata = plotdata[[
        1 + OptionValue[ignoreFirst];;
        -(1 + OptionValue[ignoreLast])
        ]];
		midpoint = Ceiling[numPoints / 2];
		rangeSize = (Min[plotdata] - Max[plotdata])/10;
		range = {Min[plotdata] + rangeSize, Max[plotdata] - rangeSize};
		debugPrint["making plot with options"];
		debugPrint[FilterRules[{opts}, Options[ListLinePlot]]];
    If[dims == 1,
      graph = Show[
        ListLinePlot[plotdata,
          FilterRules[{opts}, Options[ListLinePlot]],
          PlotRange -> {{1, numPoints}, range}
        ],
        Graphics[Line[Transpose[{{midpoint, midpoint}, range}]]]],
      If[dims == 2,
        graph =
          ListPlot3D[plotdata,
            PlotStyle->Automatic,
            FilterRules[{opts}, Options[ListPlot3D]],
            PlotRange -> {{1, numPoints}, {1, numPoints}, range}
          ]
        ]
      ];
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
	    {plotArrayPES, ListLinePlot, convertToeV, convTDDFT, makePlotPES, ListContourPlot3D}
		]]:=
	Module[{data,tempopts},
		data = importArrayPES[fileName];
		tempopts = {opts} ~Join~ Options[plotArrayPES];
		(*tempopts is the set of options to be passed to all subfunctions
		after appropriate filtering.*)
		If[OptionValue[gausTDDFT],
			data = convTDDFT[data, FilterRules[{tempopts}, Options[convTDDFT]]],
			If[data[[1,2]] > 0,
				Print["Seems like file may be Gaussian TDDFT output."];
				If[Input["Seems like file may be Gaussian TDDFT output."<>
					"\nProcess accordingly?",True],
					data = convTDDFT[data, FilterRules[{tempopts}, Options[convTDDFT]]];
					data = data / 27.21139 (*because it still goes to convToeV after*),
					"",
					Print["input not recognized, assuming 'False'"]
				],
				""
			];
			If[OptionValue[convToeV],
				data = convertToeV[data, FilterRules[{tempopts}, Options[convertToeV]]];
				AppendTo[tempopts, {FrameLabel -> {"", "Rel. Energy / eV"}}]
			]
		];
		debugPrint["Did data conversions, sending to makePlot with opts"];
		debugPrint[FilterRules[{tempopts},
				{Options[makePlotPES], Options[ListLinePlot]}
			]];
		debugPrint["All opts:", tempopts];
		makePlotPES[data, array->True,
			FilterRules[{tempopts},
				{Options[makePlotPES], Options[ListLinePlot], Options[ListPlot3D]}
			]
		]
	]




End[] (* End Private Context *)

EndPackage[]
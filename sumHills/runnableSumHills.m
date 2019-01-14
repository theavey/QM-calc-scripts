#!/usr/local/bin/MathematicaScript -script

(* This function needs to be executable. It can be called with 0, 1, or 2 arguments.
If provided, the first argument is the name of the HILLS file which is HILLS by default.
If provided, the second argument is the name of the variable of the output data.
Again, HILLS by default.
It will save everything as mathematicaHILLS.m in the current directory. *)

If[Length[$ScriptCommandLine] > 1,
  If[Length[$ScriptCommandLine] > 2,
    If[Length[$ScriptCommandLine] > 3,
      Print["Ignoring unknown positional arguments ",
        $ScriptCommandLine[[4;;]]],
      ""
    ];
    varName = ToExpression[$ScriptCommandLine[[3]]],
    varName = Automatic
  ];
  inFileName = $ScriptCommandLine[[2]],
  inFileName = "HILLS";
  varName = Automatic
];

(* ::Package:: *)

SetOptions[$Output, FormatType -> OutputForm];

LaunchKernels[1];

(* Mathematica Package         *)
(* Created by IntelliJ IDEA    *)

(* :Title: sumHillsFofT     *)
(* :Context:                *)
(* :Author: Thomas Heavey   *)
(* :Date: 8/24/15           *)

(* :Package Version: 0.3.1.0   *)
(* :Mathematica Version: 9     *)
(* :Copyright: (c) 2015 Thomas Heavey *)
(* :Keywords:                  *)
(* :Discussion:                *)


Options[sumHills] =
    {
      GridSize -> 0.1,
    (* 1000 is chunking to tenths of nanoseconds for our original default parameters *)
      TimeChunkSize -> 1000,
      name -> Automatic
    };

sumHills[hillsFileName_, OptionsPattern[]]:=
    Module[
      {
        numofCVs,
        rawData,
        sigmaCVs,
        minMaxCVs,
        gridLengthCVs,
        gridCVs, gridAllD, flatGridAllD,
        gridSize,
        timeChunk, timeChunkwUnits,
        gaussianMatrix,
        scaledRotatedGaussMat,
        filler,
        processData, processedData,
        variableName,
        partSpec,
        compAccum,
        compAddGrid,
        compTimeChunkFunc,
        periodic
      },
    (* Assign name for output of data *)
      variableName = If[
        OptionValue[name] === Automatic,
      (* Take data file name and keep only alphanumerics. *)
        ToExpression[StringReplace[hillsFileName, Except[WordCharacter] -> ""]],
      (* Use given name. *)
        OptionValue[name]];
      Print["Data will be output as ", ToString[variableName]];
      (* Import data, checking for comments and empty elements*)
      Print["Importing data from ", hillsFileName];
      rawData = DeleteCases[#, {_String, __} | {}]& @ Import[hillsFileName, "Table"];
      If[rawData == $Failed,
        Print[
          StringForm["!! Data import failed!! Incorrect file name? Given ``", hillsFileName]];
        Return[$Failed]];
      Print["Data imported successfully"];
      numofCVs = (Length[rawData[[1]]] - 3) / 2;
      Print["Number of CVs found in HILLS file: ", numofCVs];
      sigmaCVs = Table[rawData[[1, 1 + numofCVs + i]], {i, numofCVs}];
      minMaxCVs = Table[
        {Min[rawData[[All, 1 + i]]], Max[rawData[[All, 1 + i]]]},
        {i, numofCVs}];
      gridSize = OptionValue[GridSize];
      timeChunk = OptionValue[TimeChunkSize];
      timeChunkwUnits = Round[timeChunk * rawData[[1, 1]]];
      DistributeDefinitions[timeChunkwUnits];
      (* Find size (dimensions) of grid needed.
      This should agree with what is generated for gridCVs (it would be a problem if it didn't),
      but I don't see why it wouldn't. Could add checking for this. Probs a good ideas. *)
      gridLengthCVs = Table[
        Ceiling[(minMaxCVs[[i, 2]] - minMaxCVs[[i, 1]]) / gridSize],
        {i, numofCVs}];
      (* Values along grid axes *)
      gridCVs = Table[
        Round[Table[i, Evaluate[{i, ## & @@ minMaxCVs[[j]], gridSize}]], gridSize],
        {j, numofCVs}];
      Table[If[gridLengthCVs[[i]] != Length[gridCVs[[i]]],
        Print[StringForm[sumHills::griderror, gridLengthCVs[[i]], Length[gridCVs[[i]]], i]]],
        {i, numofCVs}];
      (* Test periodicity by summing the bounds for each CV.
         If the sum of the min and max of a CV is about 0 (in the range -0.1 to 0.1),
         this will return True.
         Seem like it will only not work if tiny CVs are used, or something else that goes
         negative besides periodic CVs, but I can't think of anything like that currently.
         I'm sure it's possible though because of the crazy many ways you can define CVs
         with PLUMED. *)
      periodic = Table[
        IntervalMemberQ[Interval[{-0.1, 0.1}], Total[minMaxCVs[[i]]]],
        {i, numofCVs}];
      Print["Found grid parameters:"];
      Table[Print[StringForm["  Collective variable `` range: `` Ang or Rad", i, minMaxCVs[[i]]]], {i, numofCVs}];
      Table[Print[StringForm["  Is CV `` periodic? ``", i, periodic[[i]]]], {i, numofCVs}];
      Print[StringForm["  Grid dimensions: ``", gridLengthCVs]];
      Print[StringForm["  Size of time chunks: `` ps", timeChunkwUnits]];
      (* Create gaussian matrix that will be translated as needed later. *)
      gaussianMatrix = Chop[GaussianMatrix[
        {gridLengthCVs,
          sigmaCVs / gridSize},
        Method -> "Gaussian"]
          * ((2 Pi)^(numofCVs / 2)) Apply[Times, sigmaCVs] / gridSize^numofCVs,
        10^-100];
      (* grid in all dimensions. Expected to throw an irrelevant error,
      which is why it's wrapped in Quiet. *)
      gridAllD = Quiet[Array[
        Evaluate[Table[gridCVs[[i, Slot[i]]], {i, numofCVs}]] &,
        gridLengthCVs],
        {Part::pkspec1, Part::pspec}];
      (* Blank filler that is a list the same length as rawData elements *)
      filler = Table[0., {Evaluate[(3 + 2 * numofCVs)]}];
      (* Makes a list the length of the number of points in gaussianMatrix,
        then makes it into the same shape as gaussianMatrix, then takes the
        part that we care about, and Flattens it. This is the parts of the
        Flattened summed gaussian matrices we care about. Use as a Part
        specification. *)
      partSpec =
          Flatten[ArrayReshape[
            Range[Times @@ (2 * gridLengthCVs + 1)], (2 * gridLengthCVs +
                1)][[## & @@ Table[1 ;; gridLengthCVs[[i]], {i, numofCVs}]]]];
      (* Compiled function to act on time chunks.
      First, moves gaussianMatrix based on the CVs read from the
      rawData. Then scales the gaussian based on the height
      (second to last column of raw data). Next, totals all the
      gaussians within the time chunk, flattens the totaled gaussians,
      and takes the relevant part based on partSpec defined above.
      Being Listable, it can work on the time chunks separately
      (in parallel if that option is added). *)
      (* Because this did not work on periodic systems (because
      it didn't wrap around the edges, obviously), I know try to
      check for that. If it's not periodic, use the old definition.
      If at least one of them is periodic, do something different*)
      If[Not[Or @@ periodic],
        Print["Using function for non-periodic CVs"];
        compTimeChunkFunc =
            With[{numofCVs = numofCVs, gridLengthCVs = gridLengthCVs,
              gaussianMatrix = gaussianMatrix, minMaxCVs = minMaxCVs,
              gridSize = gridSize, partSpec = partSpec},
              Compile[{{chunkedData, _Real, 2}},
                Flatten[Total[
                (* Move the gaussian matrix (scaled appropriatly by the height
                 which is the second to last column of the rawData) so that it is
                 centered at the correct coordinates. *)
                  RotateLeft[-gaussianMatrix * #[[-2]],
                    Round[Table[
                      gridLengthCVs[[i]] - (#[[1 + i]] - minMaxCVs[[i, 1]])/
                          gridSize, {i, numofCVs}]]]
                      & /@ chunkedData]][[partSpec]],
                RuntimeAttributes -> {Listable},
                Parallelization -> True,
                CompilationTarget -> "C"]],
        If[numofCVs == 1,
          Print["Using function for 1 periodic CV"];
          compTimeChunkFunc =
              With[{numofCVs = numofCVs, gridLengthCV = gridLengthCVs[[1]],
                gaussianMatrix = gaussianMatrix, minMaxCVs = minMaxCVs,
                gridSize = gridSize, partSpec = partSpec},
                Compile[{{chunkedData, _Real, 2}},
                  Total[Partition[
                    Flatten[Total[
                    (* Move the gaussian matrix (scaled appropriatly by the height
                      which is the second to last column of the rawData) so that it is
                      centered at the correct coordinates. *)
                      RotateLeft[-gaussianMatrix * #[[-2]],
                        Round[
                          gridLengthCV - (#[[2]] - minMaxCVs[[1, 1]])/
                              gridSize]]
                          & /@ chunkedData
                    ]],
                    gridLengthCV, gridLengthCV, {1, 1}, 0
                  ]],
                  {{Partition[__], _Real, 2}},
                  RuntimeAttributes -> {Listable},
                (*Parallelization -> True,*)
                  CompilationTarget -> "C"]],
          Print[StringForm["!! I don't know how to deal with this data! \
        The number of CVs found was ``, and the periodicity found was ``",
            numofCVs, periodic]];
          Return[$Failed]
        ]
      ];
      (* Flattened full dimensional grid *)
      flatGridAllD = Flatten[gridAllD, numofCVs - 1];
      (* Compiled Accumulate function *)
      compAccum = Compile[{{totaledChunks, _Real, 2}},
        Accumulate[totaledChunks],
        CompilationTarget -> "C"];
      (* Compiled function to add grid to height data.
      Being listable, it works on time chunks separately,
      so it can be run in parallel if that option is added. *)
      compAddGrid = With[{flatGrid = flatGridAllD, numofCVs = numofCVs},
        Compile[{{accumedData, _Real, 1}},
          MapThread[Append, {flatGrid, accumedData}],
          RuntimeAttributes -> {Listable},
          Parallelization -> True,
          CompilationTarget -> "C"]];
      processData =
          Function[data,
          (* Flatten the grid to numofCVs-1 level, then append
            the Flattened list of heights to that, and do that for each time
            point. *)
            compAddGrid[
            (* Sum all previous time points for each time. *)
              compAccum[
                (compTimeChunkFunc[
                (* Partition into chunks of size timeChunk, non-overlapping,
                  no overhang, padded with a row of 0. if needed *)
                  Partition[data, timeChunk, timeChunk, {1, 1}, {filler}]
                ])
              ]
            ]
          ];
      Print["Processing data..."];
      processedData = processData[rawData];
      Print["Done processing data"];
      (* Set downvalues of output *)
      Evaluate[variableName][getData] = processedData;
      Evaluate[variableName][getMinMax] = minMaxCVs;
      (* gridSize is grid spacing *)
      Evaluate[variableName][getGridSize] = gridSize;
      Evaluate[variableName][getGrid] = gridCVs;
      Evaluate[variableName][getTimeChunk] = timeChunk;
      Evaluate[variableName][getTimeChunkwUnits] = timeChunkwUnits;
      (* Times of time chunks (only rows beginning through end by every timeChunk) *)
      Evaluate[variableName][getTimes] = rawData[[;; ;; timeChunk, 1]];
      Evaluate[variableName][getNumofCVs] = numofCVs;
      Evaluate[variableName][about] =
          {
            {"Number of CVs ", numofCVs},
            {"Number of points per time chunk ", timeChunk},
            {"Picoseconds per time chunk ", timeChunkwUnits},
            {"Number of time points (chunks) ", Length[processedData]},
            {"Grid spacing (angstroms or rad) ", gridSize},
            {"Dimensions of grid ", gridLengthCVs},
            {"Originally processed on ", Date[]},
            {"Downvalues originally assigned: ",
              {getData, getMinMax, getGridSize,
                getGrid, getTimeChunk, getTimeChunkwUnits,
                getTimes, getNumofCVs, about}}
          };
      (* Set upvalues of output *)
      (* todo change the upvalues based on numofCVs (or update functions to take different numbers of CVs) *)
      Evaluate[variableName] /: Plot[Evaluate[variableName],
        opts:OptionsPattern[plotHills]] :=
          plotHills[Evaluate[variableName], opts];
      Evaluate[variableName] /: Plot[Evaluate[variableName]] :=
          plotHills[Evaluate[variableName]];
      Evaluate[variableName] /: Plot[Evaluate[variableName],
        {a_, b_},
        opts:OptionsPattern[plotHillsPoint]] :=
          plotHillsPoint[Evaluate[variableName], {a, b}, opts];
      Evaluate[variableName] /: Plot[Evaluate[variableName], {a_, b_}] :=
          plotHillsPoint[Evaluate[variableName], {a, b}];
      Evaluate[variableName] /: Plot[Evaluate[variableName], "diff",
        opts:OptionsPattern[plotHillsDiff]] :=
          plotHillsDiff[Evaluate[variableName], opts];
      Evaluate[variableName] /: Plot[Evaluate[variableName], "diff"] :=
          plotHillsDiff[Evaluate[variableName]];
      variableName
    ] (* Don't put a semicolon here *)

sumHills::griderror = "gridLength (`1`) does not match the length of generated grid (`2`) for CV `3`";

Options[plotHills] =
    {
      manipulate -> True,
      timePoint -> Automatic,
      PlotRange -> All,
      ColorFunction -> "TemperatureMap",
      ## & @@
          Options[ListPlot3D],
      ## & @@
          Options[Manipulate]
    };

plotHills[dataName_, opts:OptionsPattern[plotHills]]:=
    Module[
      {
        data,
        tempopts,
        timepoint,
        timeLength,
        plot
      },
      tempopts = {opts} ~Join~ Options[plotHills];
      data = dataName[getData];
      timeLength = Length[data];
      timepoint = If[
        OptionValue[timePoint] === Automatic,
        -1,
        If[
          Abs[OptionValue[timePoint]] > timeLength,
          Print["Given timepoint not in data, using last point"];
          -1,
          OptionValue[timePoint]
        ]
      ];
      If[OptionValue[manipulate],
        plot = Manipulate[
          ListPlot3D[data[[i]],
            FilterRules[{tempopts}, Options[ListPlot3D]]],
          {{i, timepoint, "Time Chunk"}, 1, timeLength, 1,
            Appearance -> "Labeled"}
        (*FilterRules[{tempopts}, Options[Manipulate]]*)
        ],
        plot = ListPlot3D[data[[timepoint]],
          FilterRules[{tempopts}, Options[ListPlot3D]]]
      ];
      plot
    ]

Options[plotHillsPoint] =
    {
      dynamic -> False,
      PlotRange -> All,
      Frame -> True,
      LabelStyle -> Black,
      ImageSize -> Medium,
      ## & @@ Options[ListLinePlot],
      ## & @@ Options[ListDensityPlot]
    };

plotHillsPoint[dataName_, {x_:Null, y_:Null}, opts:OptionsPattern[]]:=
    Module[
      {
        data, lastTimePoint, times,
        tempopts,
        minMaxCV1, minMaxCV2,
        xChecked, yChecked,
        nearestFunction, nearestFunctionxy,
        location, locationxy,
        plotData,
        plot
      },
      tempopts = {opts} ~Join~ Options[plotHillsPoint];
      data = dataName[getData];
      {minMaxCV1, minMaxCV2} = dataName[getMinMax];
      times = dataName[getTimes];
      (* Check x and y values. Use Mean if invalid *)
      xChecked = If[
        x === Null || ! IntervalMemberQ[Interval[minMaxCV1], x],
        Print["Invalid x coordinate."];
        Mean[minMaxCV1],
        x];
      yChecked = If[
        y === Null || ! IntervalMemberQ[Interval[minMaxCV2], y],
        Print["Invalid y coordinate."];
        Mean[minMaxCV2],
        y];
      (* Use Nearest to find best estimate of requested location on the grid. *)
      (* Arguments are 1: the data (just taking the first time point),
    2: -> Automatic maps the data onto the integers so that it gives the
    position of the nearest, not the value, so I don't have to search again
    afterwards. 3: the requested (input) coordinates.
    If two are equally close, not quite sure what it will do. *)
      nearestFunction = Nearest[data[[1]][[All, 1;;2]] -> Automatic];
      nearestFunctionxy = Nearest[data[[1]][[All, 1;;2]]];
      location =  nearestFunction[{xChecked, yChecked}];
      locationxy = data[[1]][[location, 1;;2]][[1]];
      If[OptionValue[dynamic],
      (*Print[Dimensions[data]];*)
      (*Print[Dimensions[data[[-1]]]];*)
        lastTimePoint = data[[-1]];
        DynamicModule[
          {
            spotxy = locationxy
          },
          Column[{
            Dynamic[spotxy],
            Dynamic[ListLinePlot[
              Transpose[{times,
                Flatten[
                  data[[All,
                      nearestFunction[spotxy][[1]],
                      3
                      ]]]}],
              FrameLabel -> {"Time / ps", "Free Energy"},
              FilterRules[{tempopts}, Options[ListLinePlot]]]],
            Show[
              ListDensityPlot[lastTimePoint,
                ColorFunction -> "TemperatureMap",
                FilterRules[{tempopts}, Options[ListDensityPlot]]],
              Graphics[Locator[Dynamic[spotxy]]
              ]]}],
          UnsavedVariables -> {spotxy}
        ],
        Print["Taking data at point ", locationxy];
        (* From All times, take determined location, then just take the height there. *)
        plotData = Transpose[{times, Flatten[data[[All, location, 3]]]}];
        plot = ListLinePlot[plotData,
          FrameLabel -> {"Time", "Free Energy"},
          FilterRules[{tempopts}, Options[ListLinePlot]]];
        Return[plot]
      ]
    (* If I add anything here, need to change above because this is set up to return
    the value of the the "If" statement, so adding a semicolon will cause a dynamic
    plot not to be returned. *)
    ]

Options[plotHillsDiff] =
    {
      ColorFunction -> "TemperatureMap",
      ## & @@ Options[ListPlot3D]
    };

plotHillsDiff[dataName_, opts:OptionsPattern[]] :=
    Module[
      {tempOpts, data, numTimePoints},
      tempOpts = {opts} ~Join~ Options[plotHillsDiff];
      data = dataName[getData];
      numTimePoints = Length[data];
      Manipulate[
        ListPlot3D[
          Transpose[{
            data[[1, All, 1]],
            data[[1, All, 2]],
            data[[i, All, 3]] - data[[i + timeDiff, All, 3]]
          }],
          FilterRules[{tempOpts}, Options[ListPlot3D]]
        ],
        {{i, 1, "Ref. Time Point"},
          1, numTimePoints - timeDiff, 1, Appearance -> "Labeled"},
        {{timeDiff, 5, "Diff. in Time"},
          1, numTimePoints - i, 1, Appearance -> "Labeled"}
      ]
    ]


Print["Loaded sum hills package, applying..."];

output = sumHills[inFileName, name -> Evaluate[varName]]

Print["sumHills complete, trying to save file..."];

FullDefinition[output] >> "mathematicaHILLS.m"

Print["File saved as mathematicaHILLS.m; done; quitting..."];

Quit[]

(* Mathematica Package         *)
(* Created by IntelliJ IDEA    *)

(* :Title: sumHillsFofT     *)
(* :Context: sumHillsFofT`  *)
(* :Author: Thomas Heavey   *)
(* :Date: 09/01/15          *)

(* :Package Version: 0.3.2.0   *)
(* :Mathematica Version: 9     *)
(* :Copyright: (c) 2015 Thomas Heavey *)
(* :Keywords:                  *)
(* :Discussion:                *)

BeginPackage["sumHillsFofT`"]
(* Exported symbols added here with SymbolName::usage *)

sumHills::usage = "sumHills[HILLS_file, options] returns a list of 2D arrays that
  are the time steps of the growth of the height of the negative of the free energy
  surface from a PLUMED metadynamics calculation.
  Only made for 2 collective variables currently, but that can't be changed by
  making it check for the length of a row in the input data.";

plotHills::usage = "plotHills[name of HILLS variable, options] Takes output of sumHills
  and plots time steps.";

plotHillsPoint::usage = "plotHillsPoint[name of HILLS variable, {x, y}, options] takes output of
  sumHills and plots the selected point as a function of time.";

plotHillsDiff::usage = "plotHillsDiff[name of HILLS variable] returns a plot that can
  be manipulated of the difference between two time points along a HILLS trajectory";

importColvar::usage = "importColvar[file name] imports a COLVAR file and returns the data";

plotColvar::usage = "plotColvar[colvar data] plots a COLVAR data set with Manipulate";

plot2Colvar::usage = "plot2Colvar[colvar data 1, colvar data 2] plots 2 COLVAR data
  sets side by side with Manipulate";

plotSSR::usage = "plotSSR[{data1, data2, ...}] plots the standard deviation of the
  height of binned data. Takes imported COLVAR data and returns the plot.
  Takes options for ListLinePlot and ssr";

plotHillsSSR::usage = "plotHillsSSR[name of HILLS variable, options] plots the
  standard deviation of the difference of two time points from the summed HILLS file
  as a function of time. With dynamic -> True (default), you can Manipulate the
  amount of time between the two time points.";

plotDensityHillsSSR::usage = "plotDensityHillsSSR[name of HILLS var, options]
  plots the standard deviation of the difference of two time points as a function
  of both time and time difference as a flat colored plot.";

plot3DHillsSSR::usage = "plot3DHillsSSR[name of HILLS var, options]
  plots the standard deviation of the difference of two time points as a function
  of both time and time difference as a 3D plot.";

fixOldSummedHills::usage = "fixOldSummedHills[name of HILLS var, options]
  This will fix old Private definions of getData, etc. to make them
  not Private";

fixPySummedHills::usage = "fixPySummedHills[name of HILLS var, options]
  This will add downvalues for grid min and max, grid, and
  grid spacing (grid size)";

(* Begin Private Context *)
Begin["`Private`"];
(* todo create function to fix old processed summed Hills variables (sum...`Pr..`getData vs. just getData) use Share[]? *)
(* todo delete the import functionality in this package? only useful for testing? *)

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
      variableName[getData] = processedData;
      variableName[getMinMax] = minMaxCVs;
        (* gridSize is grid spacing *)
      variableName[getGridSize] = gridSize;
      variableName[getGrid] = gridCVs;
      variableName[getTimeChunk] = timeChunk;
      variableName[getTimeChunkwUnits] = timeChunkwUnits;
        (* Times of time chunks (only rows beginning through end by every timeChunk) *)
      variableName[getTimes] = rawData[[;; ;; timeChunk, 1]];
      variableName[getNumofCVs] = numofCVs;
      variableName[about] =
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
      variableName /: Plot[variableName,
        opts:OptionsPattern[plotHills]] :=
          plotHills[variableName, opts];
      variableName /: Plot[variableName] :=
          plotHills[variableName];
      variableName /: Plot[variableName,
        {a_, b_},
        opts:OptionsPattern[plotHillsPoint]] :=
          plotHillsPoint[variableName, {a, b}, opts];
      variableName /: Plot[variableName, {a_, b_}] :=
          plotHillsPoint[variableName, {a, b}];
      variableName /: Plot[variableName, "diff",
        opts:OptionsPattern[plotHillsDiff]] :=
          plotHillsDiff[variableName, opts];
      variableName /: Plot[variableName, "diff"] :=
          plotHillsDiff[variableName];
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
        timeLength,
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

importColvar[fileName_String] := Module[{},
  (DeleteCases[#, {_String, __} | {}] &@
      Import[fileName, "Table"])[[All, 2 ;; 3]]]

plot2Colvar[data1_, data2_] := Module[{},
  Manipulate[
    GraphicsRow[{
      ListPlot[data1[[i ;; i + number]],
        AspectRatio -> 1,
        PlotRange -> {{0, 10}, {0, 10}},
        ImageSize -> Medium,
        PlotLabel -> "first data set"],
      ListPlot[data2[[i ;; i + number]],
        AspectRatio -> 1,
        PlotRange -> {{0, 10}, {0, 10}},
        ImageSize -> Medium,
        PlotLabel -> "second data set"]
    }],
    {{i, 1, "Position of Points"},
      1, Length[data1] - number, 10000,
      Appearance -> "Labeled"},
    {{number, 50000, "Number of Points to Plot"},
      1000, 100000, 1000,
      Appearance -> "Labeled"}]]

plotColvar[data_] := Module[{},
  Manipulate[
    ListPlot[data[[i ;; i + number]],
      AspectRatio -> 1,
      PlotRange -> {{0, 10}, {0, 10}},
      ImageSize -> Medium],
    {{i, 1, "Position of Points"},
      1, Length[data] - number, 10000,
      Appearance -> "Labeled"},
    {{number, 50000, "Number of Points to Plot"},
      1000, 100000, 1000,
      Appearance -> "Labeled"}]]

Options[ssr] =
    {
      binSize -> 0.1,
      silent -> True,
      minDimension -> 100
    };

ssr[data_, opts : OptionsPattern[]] :=
    Module[{binned, mins, maxs, arraySize, flatArray, binFactor},
      binFactor = 1 / OptionValue[binSize];
      mins = Min /@ Transpose[data];
      maxs = Max /@ Transpose[data];
      (* Find array size in both dimensions *)
      arraySize =
          Max[{OptionValue[minDimension], #}] & /@
              IntegerPart[(maxs - mins) * binFactor];
      If[OptionValue[silent], "", Print["Size of array: ", arraySize]];
      flatArray =
          Flatten[Array[{##} &, arraySize, IntegerPart[mins * binFactor]], 1];
      binned =
          Tally[flatArray ~Join~ (IntegerPart /@ (binFactor * data))][[All, 2]] - 1;
      N[StandardDeviation[Flatten[binned]]]]
(* todo change minDimension to minRange or something like that that accounts for binSize/binFactor *)

Options[plotSSR] = {
  ssrSize -> 40000,
  fineness -> 1000,
  PlotRange -> All,
  ## & @@ Options[ssr],
  ## & @@ Options[ListLinePlot]};

plotSSR[datas__, opts : OptionsPattern[]] :=
    Module[
      {plotData, chunkSize, tempOpts, numberofData, dataSpacing},
      tempOpts = opts ~Join~ Options[plotSSR];
      chunkSize = OptionValue[ssrSize];
      dataSpacing = OptionValue[fineness];
      numberofData = Length[datas];
      plotData =
          Transpose[Table[
            ssr[#[[i ;; i + chunkSize]],
              FilterRules[{tempOpts}, Options[ssr]]
            ] & /@ datas,
            {i, 1, Length[datas[[1]]] - chunkSize, dataSpacing}
          ]];
      ListLinePlot[plotData,
        PlotLegends -> Range[numberofData],
        FilterRules[{tempOpts}, Options[ListLinePlot]]
      ]
    ]


hillsSSR = Compile[{{hillsTP1, _Real, 1}, {hillsTP2, _Real, 1}},
  StandardDeviation[hillsTP1 - hillsTP2]
  (*, CompilationTarget -> "C"*)
  (* At least on my computer, it's a little slower compiled to C code *)
];

Options[plotHillsSSR] = {
  timeDifference -> 10,
  fineness -> 1,
  PlotRange -> All,
  dynamic -> True,
  ## & @@ Options[hillsSSR],
  ## & @@ Options[ListLinePlot]};

plotHillsSSR[hillsVarName_, opts:OptionsPattern[]] :=
    Module[
      {data, plotData, timeDiff, tempOpts, dataSpacing},
      (* Combined input and default options *)
      tempOpts = {opts} ~Join~ Options[plotHillsSSR];
      (* Number of time chunks apart for taking the difference.
        Default is 10 time chunks, or 1 ns. *)
      timeDiff = OptionValue[timeDifference];
      (* How often to take differences. Default is every time chunk *)
      dataSpacing = OptionValue[fineness];
      (* data obtained from input variable name, but we only want the heights,
        so just take the third column from All rows of All time points,
        then Chop it (make small numbers exactly 0) so it can be passed
        so it can be passed to the compiled standard deviation function. *)
      data = Chop[
        hillsVarName[sumHillsFofT`Private`getData][[All, All, 3]]];
      If[OptionValue[dynamic],
        Manipulate[
          ListLinePlot[
            Table[
              hillsSSR[data[[i]], data[[i + timeDiffM]]],
              {i, 1, Length[data] - timeDiffM, dataSpacing}
            ],
            FilterRules[{tempOpts}, Options[ListLinePlot]]],
          {{timeDiffM, timeDiff, "Time Chunk Difference"},
            1, Length[data] - 1, 1,
            Appearance -> "Labeled"}
        ],
        (* Using the chopped data, make a table with the desired time points,
          and send them to hillsSSR to get the standard deviations *)
        plotData = Table[
          hillsSSR[data[[i]], data[[i + timeDiff]]],
          {i, 1, Length[data] - timeDiff, dataSpacing}
        ];
        ListLinePlot[plotData,
          FilterRules[{tempOpts}, Options[ListLinePlot]]
        ]]
    ]

Options[plotDensityHillsSSR] = {
  timeSpacing -> 1,
  diffSpacing -> 1,
  skipFirst -> 0,
  ColorFunction -> "TemperatureMap",
  ## & @@ Options[hillsSSR],
  ## & @@ Options[ListDensityPlot]};

plotDensityHillsSSR[hillsVarName_, opts : OptionsPattern[]] :=
    Module[{data, plotData, tempOpts, dataSpacing, diffFiness},
    (* Combined input and default options *)
      tempOpts = {opts} ~Join~ Options[plotDensityHillsSSR];
      dataSpacing = OptionValue[timeSpacing];
      diffFiness = OptionValue[diffSpacing];
      If[dataSpacing == diffFiness == 1,
        If[Input["This may take a while with timeSpacing = diffSpacing = 1\nContinue?", True],
          "",
          Return["Okay, try setting timeSpacing, diffSpacing, and skipFirst options"],
          Return["Okay, try setting timeSpacing, diffSpacing, and skipFirst options"]
        ],
        ""];
      data = Chop[
        hillsVarName[sumHillsFofT`Private`getData][[All, All, 3]]];
      plotData =
          PadRight[
            Table[
              Table[
                hillsSSR[data[[time]], data[[time + timeDiff]]],
                {time,
                  OptionValue[skipFirst] + 1,
                  Length[data] - timeDiff,
                  dataSpacing}],
              {timeDiff,
                1,
                Length[data] - OptionValue[skipFirst] - 1,
                diffFiness}
            ]];
      ListDensityPlot[plotData,
        FrameLabel -> {"Time", "\[CapitalDelta] Time"},
        FilterRules[{tempOpts}, Options[ListDensityPlot]]
      ]
    ]

Options[plot3DHillsSSR] = {
  timeSpacing -> 1,
  diffSpacing -> 1,
  skipFirst -> 0,
  ColorFunction -> "TemperatureMap",
  PlotRange -> All,
  ## & @@ Options[hillsSSR],
  ## & @@ Options[ListPlot3D]};

plot3DHillsSSR[hillsVarName_, opts : OptionsPattern[]] :=
    Module[{data, plotData, tempOpts, dataSpacing, diffFiness},
    (* Combined input and default options *)
      tempOpts = {opts} ~Join~ Options[plot3DHillsSSR];
      dataSpacing = OptionValue[timeSpacing];
      diffFiness = OptionValue[diffSpacing];
      If[dataSpacing == diffFiness == 1,
        If[Input["This may take a while with timeSpacing = diffSpacing = 1\nContinue?", True],
          "",
          Return["Okay, try setting timeSpacing, diffSpacing, and skipFirst options"],
          Return["Okay, try setting timeSpacing, diffSpacing, and skipFirst options"]
        ],
        ""];
      data = Chop[
        hillsVarName[sumHillsFofT`Private`getData][[All, All, 3]]];
      plotData =
          PadRight[
            Table[
              Table[
                hillsSSR[data[[time]], data[[time + timeDiff]]],
                {time,
                  OptionValue[skipFirst] + 1,
                  Length[data] - timeDiff,
                  dataSpacing}],
              {timeDiff,
                1,
                Length[data] - OptionValue[skipFirst] - 1,
                diffFiness}
            ]];
      ListPlot3D[plotData,
        AxesLabel -> {"Time", "\[CapitalDelta] Time", ""},
        FilterRules[{tempOpts}, Options[ListPlot3D]]
      ]
    ]

Options[fixOldSummedHills] = {doShare -> True, silent -> False};

fixOldSummedHills[hillsVarName_, opts : OptionsPattern[]] :=
    Module[{outputList},
    (* Check to see if it needs fixing *)
      If[
        OptionValue[silent],
        If[ValueQ[hillsVarName[getData]],
          Return[]],
        If[ValueQ[hillsVarName[getData]],
          Return["Seems to be fixed (based on getData)"]]];
      outputList = {};
      (* Fix the 5 known downvalues if they're assigned *)
      If[ValueQ[hillsVarName[sumHillsFofT`Private`getData]],
        hillsVarName[getData] =
            hillsVarName[sumHillsFofT`Private`getData];
        AppendTo[outputList, "getData"]];
      If[ValueQ[hillsVarName[sumHillsFofT`Private`getMinMax]],
        hillsVarName[getMinMax] =
            hillsVarName[sumHillsFofT`Private`getMinMax];
        AppendTo[outputList, "getMinMax"]];
      If[ValueQ[hillsVarName[sumHillsFofT`Private`getGridSize]],
        hillsVarName[getGridSize] =
            hillsVarName[sumHillsFofT`Private`getGridSize];
        AppendTo[outputList, "getGridSize"]];
      If[ValueQ[hillsVarName[sumHillsFofT`Private`getGrid]],
        hillsVarName[getGrid] =
            hillsVarName[sumHillsFofT`Private`getGrid];
        AppendTo[outputList, "getGrid"]];
      If[ValueQ[hillsVarName[sumHillsFofT`Private`getTimes]],
        hillsVarName[getTimes] =
            hillsVarName[sumHillsFofT`Private`getTimes];
        AppendTo[outputList, "getTimes"]];
      If[
        Not[OptionValue[silent]],
        Print["Fixed Values of ", outputList]];
      (* Optional thing to save memory and hopefully improve later \
        processing times *)
      If[OptionValue[doShare], Share[]];];

Options[fixPySummedHills] = {doShare -> True, silent -> False};

fixPySummedHills[hillsVarName_, opts:OptionsPattern[]] :=
    Module[
      {
        dataTP,
        numofCVs,
        minMaxCVs,
        gridSize,
        grid,
        timeChunkwUnits,
        times
      },
      (* One time point of the data *)
      dataTP = hillsVarName[getData][[1]];
      (* Data imported with the Python script will have getData
      assigned but not getMinMax *)
      If[
        Not[
          dataTP[[0]] == List &&
              hillsVarName[getMinMax] == $Failed],
        If[OptionValue[silent],
          Return[],
          Print["Seems to be fixed based on getData and getMinMax"];
          Return[]]];
      numofCVs = hillsVarName[getNumofCVs];
      minMaxCVs = Table[{Min[Data[[All, i]]], Max[dataTP[[All, i]]]},
        {i, numofCVs}];
      (* DeleteDuplicates is probably quite slow. Faster way? *)
      grid = Table[DeleteDuplicates[dataTP[[All, i]]], {i, numofCVs}];
      gridSize = Table[grid[[i, 1]] - grid[[i, 2]], {i, numofCVs}];
      If[Not[OptionValue[silent]],
        Print["I don't yet know how to get time chunk
        with units and times"]];
      hillsVarName[getMinMax] = minMaxCVs;
      hillsVarName[getGridSize] = gridSize;
      hillsVarName[getGrid] = grid;
      If[Not[OptionValue[silent]],
        Print["Set downvalues for getMinMax, getGridSize, and getGrid"]];
      If[OptionValue[doShare], Share[]];
    ];

(* End Private Context *)
End[];

EndPackage[]
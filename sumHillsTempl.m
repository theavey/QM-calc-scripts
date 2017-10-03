#Plot[{varname}] ^:= plotHills[{varname}]
 
#Plot[{varname}, opts$:OptionsPattern[plotHills]] ^:= 
#    plotHills[{varname}, opts$]
 
#{varname}[about] = {about}
 
#{varname}[getData] = {data}
 
#{varname}[getGrid] = (Print["getGrid not currently defined"]; $Failed)
 
#{varname}[getGridSize] = {spacing}
 
#{varname}[getMinMax] = (Print["getMinMax not currently defined"]; $Failed)
 
#{varname}[getNumofCVs] = {numcvs}
 
#{varname}[getTimeChunk] = {stride}
 
#{varname}[getTimeChunkwUnits] = (Print["getTimeChunkwUnits not defined currently"];
     $Failed)
 
#{varname}[getTimes] = (Print["getTimes not currently defined"]; $Failed)
 
plotHills[dataName_, opts:OptionsPattern[plotHills]] := 
    Module[{data, tempopts, timepoint, timeLength, plot}, 
     tempopts = Join[{opts}, Options[plotHills]]; data = dataName[getData]; 
      timeLength = Length[data]; timepoint = 
       If[OptionValue[timePoint] === Automatic, -1, 
        If[Abs[OptionValue[timePoint]] > timeLength, 
         Print["Given timepoint not in data, using last point"]; -1, 
         OptionValue[timePoint]]]; If[OptionValue[manipulate], 
       plot = Manipulate[ListPlot3D[data[[i]], FilterRules[{tempopts}, 
           Options[ListPlot3D]]], {{i, timepoint, "Time Chunk"}, 1, 
          timeLength, 1, Appearance -> "Labeled"}], 
       plot = ListPlot3D[data[[timepoint]], FilterRules[{tempopts}, 
          Options[ListPlot3D]]]]; plot]
 
Options[plotHills] = {manipulate -> True, timePoint -> Automatic, 
     PlotRange -> All, ColorFunction -> "TemperatureMap", 
     AlignmentPoint -> Center, AspectRatio -> Automatic, 
     AutomaticImageSize -> False, Axes -> True, AxesEdge -> Automatic, 
     AxesLabel -> None, AxesOrigin -> Automatic, AxesStyle -> {}, 
     Background -> None, BaselinePosition -> Automatic, BaseStyle -> {}, 
     BoundaryStyle -> GrayLevel[0], Boxed -> True, BoxRatios -> {1, 1, 0.4}, 
     BoxStyle -> {}, ClippingStyle -> Automatic, ColorFunction -> Automatic, 
     ColorFunctionScaling -> True, ColorOutput -> Automatic, 
     ContentSelectable -> Automatic, ControllerLinking -> Automatic, 
     ControllerMethod -> Automatic, ControllerPath -> Automatic, 
     CoordinatesToolOptions -> Automatic, DataRange -> Automatic, 
     DisplayFunction :> $DisplayFunction, Epilog -> {}, FaceGrids -> None, 
     FaceGridsStyle -> {}, Filling -> None, FillingStyle -> Opacity[0.5], 
     FormatType :> TraditionalForm, ImageMargins -> 0., ImagePadding -> All, 
     ImageSize -> Automatic, ImageSizeRaw -> Automatic, 
     InterpolationOrder -> None, LabelStyle -> {}, Lighting -> Automatic, 
     MaxPlotPoints -> Automatic, Mesh -> Automatic, 
     MeshFunctions -> {#1 & , #2 & }, MeshShading -> None, 
     MeshStyle -> Automatic, Method -> Automatic, NormalsFunction -> 
      Automatic, PerformanceGoal :> $PerformanceGoal, PlotLabel -> None, 
     PlotLegends -> None, PlotRange -> {Full, Full, Automatic}, 
     PlotRangePadding -> Automatic, PlotRegion -> Automatic, 
     PlotStyle -> Automatic, PreserveImageOptions -> Automatic, Prolog -> {}, 
     RegionFunction -> (True & ), RotationAction -> "Fit", 
     SphericalRegion -> False, TargetUnits -> Automatic, 
     TextureCoordinateFunction -> Automatic, TextureCoordinateScaling -> 
      Automatic, Ticks -> Automatic, TicksStyle -> {}, 
     TouchscreenAutoZoom -> False, VertexColors -> Automatic, 
     VertexNormals -> Automatic, ViewAngle -> Automatic, 
     ViewCenter -> Automatic, ViewMatrix -> Automatic, 
     ViewPoint -> {1.3, -2.4, 2.}, ViewRange -> All, ViewVector -> Automatic, 
     ViewVertical -> {0, 0, 1}, Alignment -> Automatic, 
     AppearanceElements -> Automatic, AutoAction -> False, 
     AutorunSequencing -> Automatic, BaselinePosition -> Automatic, 
     BaseStyle -> {}, Bookmarks -> {}, ContentSize -> Automatic, 
     ContinuousAction -> Automatic, ControlAlignment -> Automatic, 
     ControllerLinking -> Automatic, ControllerMethod -> Automatic, 
     ControllerPath -> Automatic, ControlPlacement -> Automatic, 
     ControlType -> Automatic, DefaultBaseStyle -> "Manipulate", 
     DefaultLabelStyle -> "ManipulateLabel", Deinitialization :> None, 
     Deployed -> False, Evaluator -> Automatic, Frame -> False, 
     FrameLabel -> None, FrameMargins -> Automatic, ImageMargins -> 0, 
     Initialization :> None, InterpolationOrder -> Automatic, 
     LabelStyle -> {}, LocalizeVariables -> True, Method -> {}, 
     Paneled -> True, PreserveImageOptions -> True, RotateLabel -> False, 
     SaveDefinitions -> False, ShrinkingDelay -> 0, 
     SynchronousInitialization -> True, SynchronousUpdating -> Automatic, 
     TouchscreenAutoZoom -> False, TouchscreenControlPlacement -> Automatic, 
     TrackedSymbols -> Full, UnsavedVariables :> None, 
     UntrackedVariables :> None}
 
plotHillsDiff[dataName_, opts:OptionsPattern[]] := 
    Module[{tempOpts, data, numTimePoints}, 
     tempOpts = Join[{opts}, Options[plotHillsDiff]]; 
      data = dataName[getData]; numTimePoints = Length[data]; 
      Manipulate[ListPlot3D[Transpose[{data[[1,All,1]], data[[1,All,2]], 
          data[[i,All,3]] - data[[i + timeDiff,All,3]]}], 
        FilterRules[{tempOpts}, Options[ListPlot3D]]], 
       {{i, 1, "Ref. Time Point"}, 1, numTimePoints - timeDiff, 1, 
        Appearance -> "Labeled"}, {{timeDiff, 5, "Diff. in Time"}, 1, 
        numTimePoints - i, 1, Appearance -> "Labeled"}]]
 
Options[plotHillsDiff] = {ColorFunction -> "TemperatureMap", 
     AlignmentPoint -> Center, AspectRatio -> Automatic, 
     AutomaticImageSize -> False, Axes -> True, AxesEdge -> Automatic, 
     AxesLabel -> None, AxesOrigin -> Automatic, AxesStyle -> {}, 
     Background -> None, BaselinePosition -> Automatic, BaseStyle -> {}, 
     BoundaryStyle -> GrayLevel[0], Boxed -> True, BoxRatios -> {1, 1, 0.4}, 
     BoxStyle -> {}, ClippingStyle -> Automatic, ColorFunction -> Automatic, 
     ColorFunctionScaling -> True, ColorOutput -> Automatic, 
     ContentSelectable -> Automatic, ControllerLinking -> Automatic, 
     ControllerMethod -> Automatic, ControllerPath -> Automatic, 
     CoordinatesToolOptions -> Automatic, DataRange -> Automatic, 
     DisplayFunction :> $DisplayFunction, Epilog -> {}, FaceGrids -> None, 
     FaceGridsStyle -> {}, Filling -> None, FillingStyle -> Opacity[0.5], 
     FormatType :> TraditionalForm, ImageMargins -> 0., ImagePadding -> All, 
     ImageSize -> Automatic, ImageSizeRaw -> Automatic, 
     InterpolationOrder -> None, LabelStyle -> {}, Lighting -> Automatic, 
     MaxPlotPoints -> Automatic, Mesh -> Automatic, 
     MeshFunctions -> {#1 & , #2 & }, MeshShading -> None, 
     MeshStyle -> Automatic, Method -> Automatic, NormalsFunction -> 
      Automatic, PerformanceGoal :> $PerformanceGoal, PlotLabel -> None, 
     PlotLegends -> None, PlotRange -> {Full, Full, Automatic}, 
     PlotRangePadding -> Automatic, PlotRegion -> Automatic, 
     PlotStyle -> Automatic, PreserveImageOptions -> Automatic, Prolog -> {}, 
     RegionFunction -> (True & ), RotationAction -> "Fit", 
     SphericalRegion -> False, TargetUnits -> Automatic, 
     TextureCoordinateFunction -> Automatic, TextureCoordinateScaling -> 
      Automatic, Ticks -> Automatic, TicksStyle -> {}, 
     TouchscreenAutoZoom -> False, VertexColors -> Automatic, 
     VertexNormals -> Automatic, ViewAngle -> Automatic, 
     ViewCenter -> Automatic, ViewMatrix -> Automatic, 
     ViewPoint -> {1.3, -2.4, 2.}, ViewRange -> All, ViewVector -> Automatic, 
     ViewVertical -> {0, 0, 1}}
 
Attributes[opts$] = {Temporary}
 
Attributes[a$] = {Temporary}
 
Attributes[b$] = {Temporary}
 
plotHillsPoint[dataName_, {x_:Null, y_:Null}, opts:OptionsPattern[]] := 
    Module[{data, lastTimePoint, times, tempopts, minMaxCV1, minMaxCV2, 
      xChecked, yChecked, nearestFunction, nearestFunctionxy, location, 
      locationxy, plotData, plot}, 
     tempopts = Join[{opts}, Options[plotHillsPoint]]; 
      data = dataName[getData]; {minMaxCV1, minMaxCV2} = dataName[getMinMax]; 
      times = dataName[getTimes]; xChecked = 
       If[x === Null ||  !IntervalMemberQ[Interval[minMaxCV1], x], 
        Print["Invalid x coordinate."]; Mean[minMaxCV1], x]; 
      yChecked = If[y === Null ||  !IntervalMemberQ[Interval[minMaxCV2], y], 
        Print["Invalid y coordinate."]; Mean[minMaxCV2], y]; 
      nearestFunction = Nearest[data[[1]][[All,1 ;; 2]] -> Automatic]; 
      nearestFunctionxy = Nearest[data[[1]][[All,1 ;; 2]]]; 
      location = nearestFunction[{xChecked, yChecked}]; 
      locationxy = data[[1]][[location,1 ;; 2]][[1]]; 
      If[OptionValue[dynamic], lastTimePoint = data[[-1]]; 
        DynamicModule[{spotxy = locationxy}, Column[{Dynamic[spotxy], 
           Dynamic[ListLinePlot[Transpose[{times, Flatten[data[[All,
                 nearestFunction[spotxy][[1]],3]]]}], FrameLabel -> 
              {"Time / ps", "Free Energy"}, FilterRules[{tempopts}, 
              Options[ListLinePlot]]]], Show[ListDensityPlot[lastTimePoint, 
             ColorFunction -> "TemperatureMap", FilterRules[{tempopts}, 
              Options[ListDensityPlot]]], Graphics[Locator[
              Dynamic[spotxy]]]]}], UnsavedVariables -> {spotxy}], 
       Print["Taking data at point ", locationxy]; 
        plotData = Transpose[{times, Flatten[data[[All,location,3]]]}]; 
        plot = ListLinePlot[plotData, FrameLabel -> {"Time", "Free Energy"}, 
          FilterRules[{tempopts}, Options[ListLinePlot]]]; Return[plot]]]
 
Options[plotHillsPoint] = {dynamic -> False, PlotRange -> All, Frame -> True, 
     LabelStyle -> GrayLevel[0], ImageSize -> Medium, 
     AlignmentPoint -> Center, AspectRatio -> GoldenRatio^(-1), Axes -> True, 
     AxesLabel -> None, AxesOrigin -> Automatic, AxesStyle -> {}, 
     Background -> None, BaselinePosition -> Automatic, BaseStyle -> {}, 
     ClippingStyle -> None, ColorFunction -> Automatic, 
     ColorFunctionScaling -> True, ColorOutput -> Automatic, 
     ContentSelectable -> Automatic, CoordinatesToolOptions -> Automatic, 
     DataRange -> Automatic, DisplayFunction :> $DisplayFunction, 
     Epilog -> {}, Filling -> None, FillingStyle -> Automatic, 
     FormatType :> TraditionalForm, Frame -> False, FrameLabel -> None, 
     FrameStyle -> {}, FrameTicks -> Automatic, FrameTicksStyle -> {}, 
     GridLines -> None, GridLinesStyle -> {}, ImageMargins -> 0., 
     ImagePadding -> All, ImageSize -> Automatic, ImageSizeRaw -> Automatic, 
     InterpolationOrder -> None, Joined -> True, LabelStyle -> {}, 
     MaxPlotPoints -> Infinity, Mesh -> None, MeshFunctions -> {#1 & }, 
     MeshShading -> None, MeshStyle -> Automatic, Method -> Automatic, 
     PerformanceGoal :> $PerformanceGoal, PlotLabel -> None, 
     PlotLegends -> None, PlotMarkers -> None, PlotRange -> Automatic, 
     PlotRangeClipping -> True, PlotRangePadding -> Automatic, 
     PlotRegion -> Automatic, PlotStyle -> Automatic, 
     PreserveImageOptions -> Automatic, Prolog -> {}, RotateLabel -> True, 
     TargetUnits -> Automatic, Ticks -> Automatic, TicksStyle -> {}, 
     AlignmentPoint -> Center, AspectRatio -> 1, Axes -> False, 
     AxesLabel -> None, AxesOrigin -> Automatic, AxesStyle -> {}, 
     Background -> None, BaselinePosition -> Automatic, BaseStyle -> {}, 
     BoundaryStyle -> None, BoxRatios -> Automatic, ClippingStyle -> None, 
     ColorFunction -> Automatic, ColorFunctionScaling -> True, 
     ColorOutput -> Automatic, ContentSelectable -> Automatic, 
     CoordinatesToolOptions -> Automatic, DataRange -> Automatic, 
     DisplayFunction :> $DisplayFunction, Epilog -> {}, 
     FormatType :> TraditionalForm, Frame -> True, FrameLabel -> None, 
     FrameStyle -> {}, FrameTicks -> Automatic, FrameTicksStyle -> {}, 
     GridLines -> None, GridLinesStyle -> {}, ImageMargins -> 0., 
     ImagePadding -> All, ImageSize -> Automatic, ImageSizeRaw -> Automatic, 
     InterpolationOrder -> None, LabelStyle -> {}, LightingAngle -> None, 
     MaxPlotPoints -> Automatic, Mesh -> None, MeshFunctions -> 
      {#1 & , #2 & }, MeshStyle -> Automatic, Method -> Automatic, 
     PerformanceGoal :> $PerformanceGoal, PlotLabel -> None, 
     PlotLegends -> None, PlotRange -> {Full, Full, Automatic}, 
     PlotRangeClipping -> True, PlotRangePadding -> Automatic, 
     PlotRegion -> Automatic, PreserveImageOptions -> Automatic, 
     Prolog -> {}, RegionFunction -> (True & ), RotateLabel -> True, 
     TargetUnits -> Automatic, Ticks -> Automatic, TicksStyle -> {}, 
     VertexColors -> Automatic}

(* Mathematica Package         *)
(* Created by IntelliJ IDEA    *)

(* :Title: geomInterp2DPackage        *)
(* :Context: geomInterp2DPackage`     *)
(* :Author: Thomas Heavey             *)
(* :Date: 5/28/15              *)

(* :Package Version: 1.0       *)
(* :Mathematica Version:       *)
(* :Copyright: (c) 2015 Thomas Heavey *)
(* :Keywords:                  *)
(* :Discussion:                *)

BeginPackage["geomInterp2DPackage`"]
(* Exported symbols added here with SymbolName::usage *)

makeGeomInterp2D::usage = "makeGeomInterp2D[middle point file name, direction 1 file name,
  direction 2 file name, export base name, [points->sqrt[total points]]]"

Begin["`Private`"] (* Begin Private Context *)

elemMasses=<|"C"->12.0107,"H"->1.00794,"O"->15.9994|>;

debugPrint[args___]:=$debugPrintFunction[args]

$debugPrintFunction=Null&;
(*$debugPrintFunction=Print[##]&*)
(*Uncomment the line above to get more (maybe) useful output during evaluation*)

Options[makeGeomInterp2D] =
    {
      points -> 21,
      trySwitch -> True,
      tryRotate -> True
    }

makeGeomInterp2D[originFile_, xFile_, yFile_, exportName_, opts:OptionsPattern[]]:=
  Module[
    {atomsO, coordsO, atomsX, coordsX, atomsY, coordsY, numAtoms, diffX, diffY, numPoints,
      coordsAll, deltaX, deltaY, numberStrings, fileNames},
    numPoints = OptionValue[points];
    {atomsO, coordsO} = importXYZ[originFile];
    {atomsX, coordsX} = importXYZ[xFile];
    {atomsY, coordsY} = importXYZ[yFile];
    debugPrint[coordsO,coordsX,coordsY];
    If[atomsO==atomsX==atomsY,,Return["Atom lists not the same!"]];
    numAtoms = Length[atomsO];
    (*this will rotate the molecules to be oriented similarly using the origin
     as the reference for both directions.*)
    If[OptionValue[tryRotate],
      debugPrint["trying rotations"];
      {atomsX, coordsX} = minDiffWith3DRot[{atomsO, coordsO}, {atomsX, coordsX}];
      {atomsY, coordsY} = minDiffWith3DRot[{atomsO, coordsO}, {atomsY, coordsY}];
      ];
    (*finally, try to switch atoms in case the ordering is different.
    This may fail (badly) is the molecular geometries are significantly different.
    Again, using the origin as the reference for both directions.*)
    If[OptionValue[trySwitch],
      {atomsX, coordsX} = switchAtomLocation[{atomsO, coordsO}, {atomsX, coordsX}];
      {atomsY, coordsY} = switchAtomLocation[{atomsO, coordsO}, {atomsY, coordsY}];
      ];
    debugPrint[coordsO,coordsX,coordsY];
    diffX = coordsO - coordsX;
    diffY = coordsO - coordsY;
    deltaX = diffX / Floor[(numPoints / 2)];
    deltaY = diffY / Floor[(numPoints / 2)];
    coordsAll = Array[coordsO + #1 * deltaX + #2 * deltaY &,
      {numPoints, numPoints}, {-Floor[numPoints / 2], -Floor[numPoints / 2]}];
    numberStrings = StringPadLeft[Table[ToString[i], {i, numPoints}],
      IntegerLength[numPoints], "0"];
    fileNames = Array[exportName <> numberStrings[[#1]] <> numberStrings[[#2]] <> ".xyz" &,
      {numPoints, numPoints}];
    Array[
      Export[fileNames[[#1, #2]],
        {atomsO, coordsAll[[#1, #2]]},
        {{"VertexTypes", "VertexCoordinates"}}] &,
      {numPoints, numPoints}
    ]
  ]


importXYZ[fileName_]:=centerOfMassCoord[Import[ToString[fileName],{{"VertexTypes", "VertexCoordinates"}}]]

diff[coords1_,coords2_] := Module[{difference},
  debugPrint["called diff"];
  difference = Total[Abs[coords1 - coords2], 2];
  debugPrint[difference];
  difference]

minDiffWith3DRot[xyz1_,xyz2_]:=
(*This function takes two XYZ file inputs and returns the 2nd,
rotated to be similar to the first. It does this by minimizing*)
    Module[{min,rotation,x,y,z,coord1,coord2,elemList, xrefl, yrefl, zrefl},
      elemList=xyz1[[1]];
      coord1=massWeightedCoords[xyz1];
      coord2=massWeightedCoords[xyz2];
      (*Try some reflections first:*)
      xrefl = {-1, 1, 1}*# & /@coord2;
      yrefl = {1, -1, 1}*# & /@coord2;
      zrefl = {1, 1, -1}*# & /@coord2;
      coord2 = If[diff[coord1, coord2] > diff[coord1, xrefl],
        Print["x reflected"];
        xrefl,
        coord2
      ];
      coord2 = If[diff[coord1, coord2] > diff[coord1, yrefl],
        Print["y reflected"];
        yrefl,
        debugPrint["no y reflection"];
        coord2
      ];
      coord2 = If[diff[coord1, coord2] > diff[coord1, zrefl],
        Print["z reflected"];
        zrefl,
        coord2
      ];
      (*Now rotate to minimize difference:*)
      debugPrint[ListPointPlot3D[{coord1,coord2},PlotStyle->PointSize[Medium]]];
      min=Minimize[Total[Abs[coord1-RotationTransform[{{1,0,0},{x,y,z}}]/@coord2],2],{x,y,z}];
      debugPrint[min[[1]]];
      rotation=min[[2]];
      debugPrint[rotation];
      coord2=RotationTransform[{{1,0,0},{x,y,z}}]/@coord2/.rotation;
      debugPrint[ListPointPlot3D[{coord1,coord2},PlotStyle->PointSize[Medium],PlotRange->All]];
      {elemList,unweightedCoords[elemList,coord2]}]


centerOfMassCoord[xyzin_]:=
(*This function takes an XYZ file and will return the xyz file transformed
to center of mass coordinates and along the principle component axes.
Depends on having elemMasses defined as an association of element symbols with masses*)
    Module[{totalMass,elemList,coM,mwCoord,centeredCoord,
      mwcenteredCoord,pcCoords,unweightedCoords,sortOrder},
      elemList=xyzin[[1]];
      totalMass=Total[elemMasses/@elemList];
      mwCoord=massWeightedCoords[xyzin];
      (*Find center of mass*)
      coM=Table[Mean[mwCoord[[All,i]]],{i,3}]/totalMass;
      debugPrint[coM];
      (*Move atoms around center of mass*)
      centeredCoord=Map[#-coM&,xyzin[[2]]];
      mwcenteredCoord=massWeightedCoords[{elemList,centeredCoord}];
      debugPrint[Table[Total[mwcenteredCoord[[All,i]]],{i,3}]];
      debugPrint[ListPointPlot3D[mwcenteredCoord]];
      (*Call helper function that will change axes to principle inertial axes*)
      pcCoords=principalAxesTransform[{elemList,centeredCoord}];
      debugPrint[ListPointPlot3D[pcCoords]];
      (*Sort the elements such that X has most mass, Y second most, Z least mass*)
      mwcenteredCoord=massWeightedCoords[{elemList,pcCoords}];
      sortOrder=Ordering[Total[Abs[mwcenteredCoord],{1}]];
      pcCoords=pcCoords[[All,Reverse[sortOrder]]];
      {elemList,pcCoords}]

principalAxesTransform[xyzin_]:=Module[{elemList,coords,massList,numAtoms,inertiaTensor,transformTensor},
  {elemList,coords}=xyzin;
  massList=elemMasses[#]&/@elemList;
  numAtoms=Length[massList];
  inertiaTensor=Sum[massList[[i]]Array[Norm[coords[[i]]]^2 KroneckerDelta[#1,#2]
      -coords[[i,#1]]*coords[[i,#2]]&,{3,3}],{i,numAtoms}];
  (*Rotation matrix is matrix of eigenvectors, which are returned as row vectors,
hence the Transpose operation*)
  transformTensor=Transpose[Eigenvectors[inertiaTensor]];
  (*Rotate and return the coordinates.*)
  #. transformTensor&/@coords]


switchAtomLocation[refxyz_,testxyz_]:=
(*This essentially assumes the molecules are oriented
very similarly. If they're very different, then the
largest difference might not need to be switched, and that
would kill the loop, while something else really does
need to be moved. Also definitely assumes that at least
the element lists are the same.*)
    Module[{elemList,refCoords,testCoords,mdIndex,mdAtomType,mdAtomCoord,sameAtomIndices,diffs,switchIndex,thisSwitch,lastSwitch,i,imax},
      elemList=refxyz[[1]];
      refCoords=refxyz[[2]];
      testCoords=testxyz[[2]];
      imax=2*Length[elemList];(*arbitrary max number of switches*)
      Do[
      (*Ordering sorts the list and returns the indices of how the list should be sorted. Thus the -1th element of ordering is the index of the largest element. It returns it as a list, so take the first (only) element of it.*)
        mdIndex=Ordering[Total[Abs[refCoords-testCoords],{2}],-1][[1]];
        {mdAtomType,mdAtomCoord}=testxyz\[Transpose][[mdIndex]];
        debugPrint[mdAtomType,mdAtomCoord];
        sameAtomIndices=Flatten[Position[elemList,mdAtomType]];
        debugPrint[sameAtomIndices];
        diffs=Total[Abs[refCoords[[#]]-mdAtomCoord]]&/@sameAtomIndices;
        debugPrint[diffs];
        switchIndex=sameAtomIndices[[Ordering[diffs,1][[1]]]];
        If[switchIndex==mdIndex,
          If[i==1,
            Print["Nothing to switch, returning input XYZ as is."],Print["No more switches to make. Made ",i-1," switches total."]];
            Break[],
          ""];
        testCoords[[{switchIndex,mdIndex}]]=testCoords[[{mdIndex,switchIndex}]];
        thisSwitch=Sort[{switchIndex,mdIndex}];
        If[thisSwitch===lastSwitch,
          Print["Uh oh. Switched atoms ",thisSwitch," again. Breaking."];Break[],
          Print["Switched atoms ",thisSwitch];
          lastSwitch=thisSwitch
        ];
        If[i==imax,
          Print["Hit max number of atom switch iterations, ",i,". Accidental repetition?"],
          ""],
        {i,imax}];
      {elemList,testCoords}
    ]


massWeightedCoords[xyzin_]:=Map[elemMasses[#[[1]]]*#[[2]]&,xyzin\[Transpose]]
unweightedCoords[elemList_,coords_]:=Map[#[[2]]/elemMasses[#[[1]]]&,{elemList,coords}\[Transpose]];



End[] (* End Private Context *)

EndPackage[]
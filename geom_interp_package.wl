(* ::Package:: *)

BeginPackage["geomInterp`"]

makeFullGeomInterp::usage=
	"makeFullGeomInterp[reactantXYZ, transitionstateXYZ, productXYZ, exportname, number of points]
	This function takes given XYZ file names and exports 2\[Times](number of points)+1 files in XYZ format.
	It imports the files, centering at the center of mass, and transforming the axes to the principal
	axes of rotation. It will then rotate both the transition state and product to be similar to the 
	reactant. It will then try to switch atoms of the same type in case they're in different orders.
	If the geometry changes significantly, this may cause some problems. Finally, it will do a linear
	interpolation from the reactant to transition state, then transition state to product, and export
	them as exportname(n).xyz where n is a sequential number."



Begin["`Private`"]

elemMasses=<|"C"->12.0107,"H"->1.00794,"O"->15.9994|>;

debugPrint[args___]:=$debugPrintFunction[args]

$debugPrintFunction=Null&;
(*$debugPrintFunction=Print[##]&*)
(*Uncomment the line above to get more (maybe) useful output during evaluation*)


geomInterp[react_,product_,numPoints_]:=Module[{diff,delta,i,output},
diff=product-react;
delta=diff/numPoints;
i=-1;
output={};
While[i<numPoints,
i++;
AppendTo[output,react+i*delta];
];
output]


importXYZ[fileName_]:=centerOfMassCoord[Import[ToString[fileName],{{"VertexTypes", "VertexCoordinates"}}]]


exportXYZ[dataWithName_]:=Module[{fileName,data},
fileName=dataWithName[[1]];
data=dataWithName[[2;;3]];
Export[ToString[fileName]<>".xyz",data,{{"VertexTypes", "VertexCoordinates"}}]]


exportSetOfXYZ[baseName_,data_]:=Module[{names,numPoints,namedData,numAtoms},
numPoints=Length[data];
numAtoms=Length[data[[All,1]]];
names=makeNames[baseName,numPoints];
namedData=PadLeft[data,{numPoints,3},{names}\[Transpose]];
Map[exportXYZ,namedData]]


makeNames[baseName_,numPoints_]:=Module[{i},
Map[ToString[baseName]<>ToString[#]&,Table[i,{i,numPoints}]]]


makeGeomInterp[reactName_,productName_,numPoints_,exportName_]:=Module[{reactCoord,prodCoord,atomsR,atomsP,interpedData,dataWithAtoms,numAtoms},
{atomsR,reactCoord}=importXYZ[reactName];
{atomsP,prodCoord}=importXYZ[productName];
If[atomsR==atomsP,,Return["Atom lists not the same!"]];
numAtoms=Length[atomsR];
interpedData=Partition[geomInterp[reactCoord,prodCoord,numPoints],1];
dataWithAtoms=PadLeft[interpedData,{numPoints+1,2,numAtoms},atomsR];
exportSetOfXYZ[exportName,dataWithAtoms]]


makeFullGeomInterp[reactName_String,tsName_String,productName_String,exportName_String,numPoints_Integer:10]:=
Module[{reactCoord,tsCoord,prodCoord,atomsR,atomsT,atomsP,interpedDataRT,interpedDataTP,interpedData,dataWithAtoms,numAtoms},
(*importing with the importXYZ function will transform the axis origin to the center of mass and make the axes the principal mass weighted axes.*)
{atomsR,reactCoord}=importXYZ[reactName];
{atomsP,prodCoord}=importXYZ[productName];
{atomsT,tsCoord}=importXYZ[tsName];
If[atomsR==atomsP==atomsT,,Return["Atom lists not the same!"]];
numAtoms=Length[atomsR];
(*this will rotate the molecules to be oriented similarly using the reactant as the reference for both the transition state and product.*)
{atomsP,prodCoord}=minDiffWith3DRot[{atomsR,reactCoord},{atomsP,prodCoord}];
{atomsT,tsCoord}=minDiffWith3DRot[{atomsR,reactCoord},{atomsT,tsCoord}];
(*finally, try to switch atoms in case the ordering is different. This may fail (badly) is the molecular geometries are significantly different. Again, using the reactant as the reference for both the product and transition state.*)
{atomsP,prodCoord}=switchAtomLocation[{atomsR,reactCoord},{atomsP,prodCoord}];
{atomsT,tsCoord}=switchAtomLocation[{atomsR,reactCoord},{atomsT,tsCoord}];
(*With preprocessing done, now call geomInterp to make the linear interpolations.*)
interpedDataRT=Partition[geomInterp[reactCoord,tsCoord,numPoints],1];
interpedDataTP=Partition[geomInterp[tsCoord,prodCoord,numPoints],1];
interpedData=Join[interpedDataRT,interpedDataTP[[2;;All]]];
dataWithAtoms=PadLeft[interpedData,{2*numPoints+1,2,numAtoms},atomsR];
exportSetOfXYZ[exportName,dataWithAtoms]]


minDiffWith3DRot[xyz1_,xyz2_]:=(*This function takes two XYZ file inputs and returns the 2nd, rotated to be similar to the first. It does this by minimizing*)Module[{min,rotation,x,y,z,coord1,coord2,elemList},
elemList=xyz1[[1]];
coord1=massWeightedCoords[xyz1];
coord2=massWeightedCoords[xyz2];
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


massWeightedCoords[xyzin_]:=Map[elemMasses[#[[1]]]*#[[2]]&,xyzin\[Transpose]]
unweightedCoords[elemList_,coords_]:=Map[#[[2]]/elemMasses[#[[1]]]&,{elemList,coords}\[Transpose]];


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


switchAtomLocation[refxyz_,testxyz_]:=(*This essentially assumes the molecules are oriented very similarly. If they're very different, then the largest difference might not need to be switched, and that would kill the loop, while something else really does need to be moved. Also definitely assumes that at least the element lists are the same.*)
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
Print["Nothing to switch, returning input XYZ as is."],Print["No more switches to make. Made ",i," switches total."]];
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

End[]

EndPackage[]





/*********************************************************************
  Single nanopore.
 *********************************************************************/

RevL = 16e-9 ;
radius = 5.4e-9 ;
RevH = 20e-9;
length = 90e-9;


lc = 1.00e-9 ;
lc2 = 2.0e-9;

Point(1) = { 0,0,0,lc};
Point(2) = { RevL,0,0,lc};
Point(3) = { RevL,RevL,0,lc};
Point(4) = { 0,RevL,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Point(5) = { RevL/2,RevL/2,0, 3e-9};
Point(6) = { RevL/2,RevL/2 - radius,0,lc2};
Point(7) = { RevL/2 - radius, RevL/2,0,lc2};
Point(8) = { RevL/2 , RevL/2 + radius ,0,lc2};
Point(9) = { RevL/2 + radius, RevL/2,0,lc2};




Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(2) = { 5,6,7,8};
Line Loop(1) = { 1,2,3,4};
Plane Surface(1) = {1,2};



Plane Surface(2) = {2};
out1[]=Extrude{0,0,-RevH}{Line{1,2,3,4};};
Line Loop(3) = { 9,13,17,21};
Plane Surface(3) = {3};
Surface Loop(1) = {1,2,3,12,20,16,24};
Volume(1) = {1};

out1[]=Extrude{0,0,length}{Surface{2};};
NewSurf1[]=Translate{0,0,length}{ Duplicata{Surface{1};}};
NewSurf2[]=Translate{0,0,length+RevH}{ Duplicata{Surface{20,12,24,16};}};
NewSurf3[]=Translate{0,0,length+2*RevH}{ Duplicata{Surface{3};}};
Surface Loop(10) = {47,46,69,67,62,52,57};
Volume(10) = {10};


//No mesh refinement yet, need to do mesh refinement 

// Define Fields.... ?


// do mesh refinement at the adjuctions of nanopore and reservoir
/*
Field[1] = Attractor;
Field[1].EdgesList = {5,6,7,8,26,27,28,29};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMax = lc;
Field[2].LcMin = lc2;
Field[2].DistMax = 2e-9;
Field[2].DistMin = 1e-9;
*/

Field[1] = Attractor;
Field[1].EdgesList = {5,6,7,8,26,27,28,29,31,32,36,40};
//Field[1].EdgesList = {31,32,36,40};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMax = lc2;
Field[2].LcMin = lc2/5;
Field[2].DistMax = 3e-9;
Field[2].DistMin = 2e-9;


/*Field[1] = Cylinder;
Field[1].Radius =2e-9;
Field[1].VIn = lc2/;
Field[1].VOut = lc/5;
Field[1].ZAxis = 1;
Field[1].XCenter = RevL/2 ;
Field[1].YCenter = RevL/2 ;
Field[1].ZCenter = length/2 ;
*/
// Now need to refine the mesh and the bottom resevoir boundary since
// we need to calculate the flux here
// Use a Box field for this refinement ?

Field[4] = Box;
Field[4].VIn = lc/2;
Field[4].VOut = lc;
Field[4].XMin = -1e-9;
Field[4].XMax = RevL + 1e-9;
Field[4].YMin = -1e-9;
Field[4].YMax = RevL + 1e-9;
Field[4].ZMax = -RevH + 2e-9;
Field[4].ZMin = -RevH - 1e-9;

Field[5] = Box;
Field[5].VIn = lc/3;
Field[5].VOut = lc;
Field[5].XMin = -1e-9;
Field[5].XMax = RevL + 1e-9;
Field[5].YMin = -1e-9;
Field[5].YMax = RevL + 1e-9;
Field[5].ZMax = RevH + length + 1e-9;
Field[5].ZMin = RevH + length - 2e-9;

Field[6] = Min;
Field[6].FieldsList = {2,4,5};
Background Field = 6;






//Compound Surface(5) = {1,2};
//out[]=Extrude{0,0,length}{Surface{2};};
//out1[]=Extrude{0,0,-RevH}{Surface{1};};
//+


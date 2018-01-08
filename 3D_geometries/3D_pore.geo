/*********************************************************************
  Single nanopore. # for KCl ionic conductance validation
 *********************************************************************/

RevL = 40e-9 ;
radius = 5.1e-9 ;
RevH = 100e-9;
length = 34e-9;


lc = 4e-9 ;
lc2 = 5e-9;
layerlc = 5e-9;

Point(1) = { 0,0,0,lc2};
Point(2) = { RevL,0,0,lc2};
Point(3) = { RevL,RevL,0,lc2};
Point(4) = { 0,RevL,0,lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Point(5) = { RevL/2,RevL/2,0, lc2};
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
out1[]=Extrude{0,0,-RevH}{Line{1,2,3,4};Layers{RevH/layerlc};};
Line Loop(3) = { 9,13,17,21};
Plane Surface(3) = {3};
Surface Loop(1) = {1,2,3,12,20,16,24};
Volume(1) = {1};

out1[]=Extrude{0,0,length}{Surface{2};};
NewSurf1[]=Translate{0,0,length}{ Duplicata{Surface{1};}};
out2[]=Extrude{0,0,RevH}{Line{48,49,50,51};Layers{RevH/layerlc};};
Line Loop(4) = { 52,56,60,64};
Plane Surface(555) = {4};
Surface Loop(2) = {46,47,55,59,63,67,555};
Volume(10) = {2};
Line(80) = {5,17};


// Define Fields.... ?

Field[4] = Attractor;
Field[4].EdgesList = {80};

// We can then create a MathEval field with a function that depends on the
// return value of the Attractr Field[4], i.e., depending on the distance to
// point 1 (here using a cubic law, with minumum element size = lc / 100)
Field[5] = MathEval;
Field[5].F = Sprintf("1e-9*(5.1- F4/1e-9)^2 + %g", lc / 20);


Field[6] = Min;
Field[6].FieldsList = {5};
Background Field = 6;




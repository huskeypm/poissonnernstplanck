//Unit Cell A Geometry

// parameter 
spacing=5.2e-9 ;
 
length = 500e-9 ; //length of nanopore 
radius = 18e-9 ; //pore radius
RevL = 300e-9 ; //Reservoir width and length
RevH = 300e-9 ; // Reservoir depth



//parameter for nanopore

lc = 2e-9 ;
lc2 = 2e-9;

Point(1) = { 0, 0, 0, lc };
Point(2) = { 0,RevL,0, lc } ;
Point(3) = { RevH, RevL, 0, lc } ;
Point(4) = { RevH, RevL/2 + radius/2, 0, lc } ;

Point(5) = { RevH + length , RevL/2 + radius/2, 0, lc } ;
Point(6) = { RevH + length , RevL, 0, lc } ;
Point(7) = { 2*RevH + length, RevL, 0, lc } ;
Point(8) = { 2*RevH + length,0, 0, lc } ;
Point(9) = { RevH +length,0, 0, lc } ;
Point(10) = { RevH + length, RevL/2 - radius/2, 0, lc } ;
Point(11) = { RevH, RevL/2 - radius/2, 0, lc } ;
Point(12) = { RevH ,0, 0, lc } ;

// 1/4 pore at


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,1};

Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12};
Plane Surface(1) = {1};



//
//No mesh refinement yet, need to do mesh refinement 

// Define Fields.... ?


// do mesh refinement at the adjuctions of nanopore and reservoir

Field[1] = Attractor;
Field[1].EdgesList = {3,5,11,9,4,10};
//Field[1].EdgesList = {31,32,36,40};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMax = lc2;
Field[2].LcMin = lc2/8;
Field[2].DistMax = 3e-9;
Field[2].DistMin = 2e-9;


// Now need to refine the mesh and the bottom resevoir boundary since
// we need to calculate the flux here
// Use a Box field for this refinement ?

Field[4] = Box;
Field[4].VIn = lc/3;
Field[4].VOut = lc;
Field[4].XMin = -1e-9;
Field[4].XMax = 3e-9;
Field[4].YMin = -1e-9;
Field[4].YMax = RevL + 1e-9;
Field[4].ZMax = 0;
Field[4].ZMin = 0;

Field[5] = Box;
Field[5].VIn = lc/3;
Field[5].VOut = lc;
Field[5].XMin = RevH*2 + length - 3e-9;
Field[5].XMax = RevH*2 + length + 1e-9;
Field[5].YMin = -1e-9;
Field[5].YMax = RevL + 1e-9;
Field[5].ZMax = 0;
Field[5].ZMin = 0;

Field[6] = Min;
Field[6].FieldsList = {2,4,5};
Background Field = 6;

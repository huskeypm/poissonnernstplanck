spacing = 5.2e-9;
length = 90e-9;
radius = 5.4e-9;
RevL = 5.2e-9+2*(5.4e-9);
RevH = 20e-9;

//parameter for nanopore

lc = 1e-9 ;
lc2 = 2e-9;

Point(1) = { 0, 0, 0, lc };
Point(2) = { 0,RevL,0, lc } ;
Point(3) = { RevL, RevL, 0, lc } ;
Point(4) = { RevL, 0, 0, lc } ;


// 1/4 pore at
Point(5) = { radius, 0, 0, lc } ;
Point(6) = { 0 , radius , 0, lc} ;
Circle(103) = {5,1,6};


Point(7) = { radius, RevL, 0, lc } ;
Point(8) = { 0 , RevL - radius , 0, lc} ;
Circle(203) = {7,2,8};

Point(9) = { RevL-radius, RevL, 0, lc } ;
Point(10) = { RevL , RevL - radius , 0, lc} ;
Circle(303) = {9,3,10};


Point(11) = { RevL - radius, 0, 0, lc } ;
Point(12) = { RevL , radius , 0, lc} ;
Circle(403) = {11,4,12};

Line(1) = {2,7};
Line(2) = {7,9};
Line(3) = {9,3};
Line(4) = {3,10};
Line(5) = {10,12};
Line(6) = {12,4};
Line(7) = {4,11};
Line(8) = {11,5};
Line(9) = {5,1};
Line(10) = {1,6};
Line(11) = {6,8};
Line(12) = {8,2};

Line Loop(1) = {103,11,-203,2,303,5,-403,8};
Plane Surface(1) = {1};
Extrude{0,0,-RevH}{Line{1,2,3,4,5,6,7,8,9,10,11,12};};
Line Loop(100) = {404,408,412,416,420,424,428,432,436,440,444,448};
Plane Surface(100) ={100};


Line Loop(2) = {6,7,403};
Line Loop(3) = {3,4,-303};
Line Loop(4) = {12,1,203};
Line Loop(5) = {9,10,-103};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};



Surface Loop(1) = {100,407,411,415,419,423,427,431,435,439,443,447,451,1,2,3,4,5};
Volume(1) = {1};

out1[]=Extrude{0,0,length}{Surface{2,3,4,5};};
NewSurf1[]=Translate{0,0,length}{ Duplicata{Surface{1};}};
NewSurf2[]=Translate{0,0,length+RevH}{ Duplicata{Surface{407,411,415,419,423,427,431,435,439,443,447,451};}};
NewSurf3[]=Translate{0,0,length+2*RevH}{ Duplicata{Surface{100};}};
Surface Loop(4) = {485,502,519,468,520,529,534,539,544,549,554,559,564,569,574,579,584,586};
Volume(10) = {4};

//Line Loop(200) = {453,454,455};
//Line Loop(201) = {470,471,472};
//Line Loop(202) = {487,488,489};
//Line Loop(203) = {504,505,506};
//Plane Surface(200) = {200};
//Plane Surface(201) = {201};
//Plane Surface(202) = {202};
//Plane Surface(203) = {203};
//
//No mesh refinement yet, need to do mesh refinement

// Define Fields.... ?


// do mesh refinement at the adjuctions of nanopore and reservoir

Field[1] = Attractor;
Field[1].EdgesList = {9,10,103,7,6,403,3,4,303,1,12,203,462,457,491,496,508,513,474,479,453,454,455,470,471,472,504,505,506,487,488,489};
//Field[1].EdgesList = {31,32,36,40};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMax = lc2;
Field[2].LcMin = lc2/8;
Field[2].DistMax = 3e-9;
Field[2].DistMin = 2e-9;


/*Field[1] = Cylinder;
Field[1].Radius =2e-9;
Field[1].VIn = lc2/10;
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
Field[4].VIn = lc/3;
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
Field[5].ZMax = length + RevH + 1e-9;
Field[5].ZMin = length + RevH -2e-9;

Field[6] = Min;
Field[6].FieldsList = {2,4,5};
Background Field = 6;
                                                                                    148,7         Bot


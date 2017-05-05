//Unit Cell B
// parameter for reservoir
spacing=5.2e-9 ;
 
length = 90e-9 ; //length of nanopore 
radius = 5.4e-9 ; //pore radius
RevL = spacing + 2*radius ; //Reservoir width and length
RevH = 20e-9 ; // Reservoir depth

lc = 1e-9 ;
lc2 = 1.99e-9;

//
Point(1) = { 0, 0, 0, lc2 };
Point(2) = { 0,RevL,0, lc2 } ;
Point(3) = { RevL, RevL, 0, lc2 } ;
Point(4) = { RevL, 0, 0, lc2 } ;


// 1/2 pore at
Point(5) = { RevL/2, 0, 0, lc2 } ;
Point(51) = { RevL/2 + radius, 0, 0, lc2 } ;
Point(52) = { RevL/2 , radius , 0, lc2} ;
Point(53) = { RevL/2 - radius , 0 , 0, lc2} ;
Circle(501) = {51,5,52};
Circle(502) = {52,5,53};


Point(6) = { 0, RevL/2, 0, lc2 } ;
Point(61) = { 0 , RevL/2 - radius, 0, lc2 } ;
Point(62) = { radius , RevL/2 , 0, lc2} ;
Point(63) = { 0 , RevL/2 + radius , 0, lc2} ;
Circle(601) = {61,6,62};
Circle(602) = {62,6,63};


Point(7) = { RevL/2, RevL, 0, lc2 } ;
Point(71) = { RevL/2 - radius, RevL, 0, lc2 } ;
Point(72) = { RevL/2 , RevL - radius , 0, lc2} ;
Point(73) = { RevL/2 + radius , RevL , 0, lc2} ;
Circle(701) = {71,7,72};
Circle(702) = {72,7,73};


Point(8) = { RevL, RevL/2, 0, lc2 } ;
Point(81) = { RevL, RevL/2 + radius, 0, 0, lc2 } ;
Point(82) = { RevL - radius , RevL/2 , 0, lc2} ;
Point(83) = { RevL, RevL/2 - radius , 0 , 0, lc2} ;
Circle(801) = {81,8,82};
Circle(802) = {82,8,83};

Line(1) = {1,53};
Line(2) = {53,5};
Line(3) = {5,51};
Line(4) = {51,4};
Line(5) = {4,83};
Line(6) = {83,8};
Line(7) = {8,81};
Line(8) = {81,3};
Line(9) = {3,73};
Line(10) = {73,7};
Line(11) = {7,71};
Line(12) = {71,2};
Line(13) = {2,63};
Line(14) = {63,6};
Line(15) = {6,61};
Line(16) = {61,1};


//Line Loop(1) = {103,11,-203,2,303,5,-403,8};
//Plane Surface(1) = {1};
Extrude{0,0,-RevH}{Line{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};};
Line Loop(1) = {4,5,-802,-801,8,9,-702,-701,12,13,-602,-601,16,1,-502,-501};
Plane Surface(1) = {1};

Line Loop(2) = { 801,802,6,7};
Line Loop(3) = { 701,702,10,11};
Line Loop(4) = { 601,602,14,15};
Line Loop(5) = { 501,502,2,3};
Line Loop(100) = { 803,807,811,815,819,823,827,831,835,839,843,847,851,855,859,863};
Plane Surface(100) = {100};

Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

Surface Loop(1) = {100,806,810,814,818,822,826,830,834,838,842,846,850,854,858,862,866,1,2,3,4,5};
Volume(1) = {1};
//Line Loop(100) = {404,408,412,416,420,424,428,432,436,440,444,448};
//Plane Surface(100) ={100};
out1[]=Extrude{0,0,length}{Surface{2,3,4,5};};
NewSurf1[]=Translate{0,0,length}{ Duplicata{Surface{1};}};
NewSurf2[]=Translate{0,0,length+RevH}{ Duplicata{Surface{806,810,814,818,822,826,830,834,838,842,846,850,854,858,862,866};}};
NewSurf3[]=Translate{0,0,length+2*RevH}{ Duplicata{Surface{100};}};
Surface Loop(4) = {910,932,954,888,955,970,975,980,985,990,995,1000,1005,1010,1015,1020,1025,1030,1035,1040,1045,1047};
Volume(10) = {4};


//
//No mesh refinement yet, need to do mesh refinement 

// Define Fields.... ?


// do mesh refinement at the adjuctions of nanopore and reservoir

Field[1] = Attractor;
Field[1].EdgesList = {601,602,501,502,701,702,801,802,912,913,890,891,868,869,934,935,939,940,944,917,918,922,895,896,900,873,874,878};
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

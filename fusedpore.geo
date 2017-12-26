/*********************************************************************
  Fused pore from Ryan
  
 *********************************************************************/


// basics parameter for the fused pore
width = 15e-9 ; //  width of the fused pore
height = 30e-9; // length of the fused pore
lwall = 13.5e-9 ; // wall width at the left side
rwall = 13.5e-9; // wall width at the right side

//---------------
RevD = 40e-9;
RevW = width + lwall + rwall;
RevH = height ;
length = 90e-9; 



lc2 = 5e-9;
layerlc = 5e-9;


Point(1) = { 0, 0, 0, lc2 };
Point(2) = { 0,height,0, lc2 } ;
Point(3) = { width, height, 0, lc2 } ;
Point(4) = { width, 0, 0, lc2 } ;

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};


Point(10)={-lwall, 0,0,lc2};
Point(11)={-lwall, height,0,lc2};
Line(5) = {1,10};
Line(6) = {10,11};
Line(7) = {11,2};

Point(12)={rwall+width, 0,0,lc2};
Point(13)={rwall+width, height,0,lc2};
Line(8) = {4,12};
Line(9) = {12,13};
Line(10) = {13,3};
Line Loop(2) = {1,-7,-6,-5};
Line Loop(3) = {8,9,10,3};

Plane Surface(2) = {2};
Plane Surface(3) = {3};

Point(20)={-lwall, 0,-RevD,lc2};
Point(21)={-lwall, height,-RevD,lc2};
Point(22)={rwall+width, 0,-RevD,lc2};
Point(23)={rwall+width, height,-RevD,lc2};
Line(21) = {20,21};
Line(22) = {21,23};
Line(23) = {23,22};
Line(24) = {22,20};
Line(25) = {20,10};
Line(26) = {21,11};
Line(27) = {22,12};
Line(28) = {23,13};
Line Loop(21) = {21,22,23,24};
Line Loop(22) = {25,6,-26,-21};
Line Loop(23) = {26,7,2,-10,-28,-22};
Line Loop(24) = {28,-9,-27,-23};
Line Loop(25) = {27,-8,4,5,-25,-24};


Plane Surface(21) = {21};
Plane Surface(22) = {22};
Plane Surface(23) = {23};
Plane Surface(24) = {24};
Plane Surface(25) = {25};
Surface Loop(2) = {1,2,3,21,22,23,24,25};
Volume(4) = {2};


out2[]=Extrude{0,0,length}{Surface{1};};
out3[]=Translate{0,0,length}{ Duplicata{Surface{1,2,3};}};


Point(30)={-lwall, 0,RevD+length,lc2};
Point(31)={-lwall, height,RevD+length,lc2};
Point(32)={rwall+width, 0,RevD+length,lc2};
Point(43)={rwall+width, height,RevD+length,lc2};
Line(91) = {30,31};
Line(92) = {31,43};
Line(93) = {43,32};
Line(94) = {32,30};
Line(95) = {30,59};
Line(96) = {31,55};
Line(97) = {43,71};
Line(98) = {32,67};
Line Loop(31) = {91,92,93,94};
Line Loop(32) = {95,-59,-96,-91};
Line Loop(33) = {96,-58,31,-64,-97,-92};
Line Loop(34) = {98,63,-97,93};
Line Loop(35) = {98,-62,33,-60,-95,-94};


Plane Surface(31) = {31};
Plane Surface(32) = {32};
Plane Surface(33) = {33};
Plane Surface(34) = {34};
Plane Surface(35) = {35};
Surface Loop(3) = {61,50,56,31,32,33,34,35};
Volume(7) = {3};
//out4[]=Extrude{0,0,RevD}{Line{107,134,136,105,132,138};Layers{RevD/layerlc};};
//Line(201) = {74,80};
//Line(202) = {76,81};

//Line Loop(7) = {157,-201,-144,-133};
//Line Loop(8) = {160,-202,-149,137};
//Line Loop(9) = {155,-201,143,-139,147,202,159,-151};
//Plane Surface(7) = {7};
//Plane Surface(8) = {8};
//Plane Surface(9) = {9};
//Surface Loop(4) = {7,8,9,130,124,135,154,158,162,142,146,150};
//Volume(6) = {4};

//Field[4] = Attractor;
//Field[4].FacesList = {111,119};

// We can then create a MathEval field with a function that depends on the
// return value of the Attractr Field[4], i.e., depending on the distance to
// point 1 (here using a cubic law, with minumum element size = lc / 100)
//Field[5] = MathEval;
//Field[5].F = Sprintf("1.5e-9*(F4/1e-9)^6 + %g", lc / 18);

Field[4] = Box;
Field[4].VIn = lc2/10;
Field[4].VOut = 2*lc2;
Field[4].XMin = -1e-9;
Field[4].XMax = 2e-9;
Field[4].YMin = -1e-9;
Field[4].YMax = RevH + 1e-9;
Field[4].ZMin = -2e-9;
Field[4].ZMax = length + 2e-9;

Field[5] = Box;
Field[5].VIn = lc2/10;
Field[5].VOut = 2*lc2;
Field[5].XMin = width - 2e-9;
Field[5].XMax = width + 1e-9;
Field[5].YMin = -1e-9;
Field[5].YMax = RevH + 1e-9;
Field[5].ZMin = -2e-9;
Field[5].ZMax = length + 2e-9;



Field[6] = Min;
Field[6].FieldsList = {4,5};
Background Field = 6;




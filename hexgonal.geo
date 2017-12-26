/*********************************************************************
  Hexgonal Unit Cell from Ryan
  
 *********************************************************************/


// basics parameter for the unitcell
radius = 5e-9 ; // radius of pore
hsp = 6e-9 ; // horizontal spacing
vsp = 21e-9; // vertical spacing

//---------------
RevW = hsp + 2*radius; // width of reservoir
RevH = vsp + 2*radius ; // height of reservoir
RevD = 40e-9 ; // depth of reservoir
length = 90e-9; 


lc2 = 5e-9;
layerlc = 5e-9;


Point(1) = { 0, 0, 0, lc2 };
Point(2) = { 0,RevH,0, lc2 } ;
Point(3) = { RevW, RevH, 0, lc2 } ;
Point(4) = { RevW, 0, 0, lc2 } ;


// 1/2 pore at
Point(5) = { RevW/2, 0, 0, lc2 } ;
Point(51) = { RevW/2 + radius, 0, 0, lc2 } ;
Point(52) = { RevW/2 , radius , 0, lc2} ;
Point(53) = { RevW/2 - radius , 0 , 0, lc2} ;
Circle(501) = {51,5,52};
Circle(502) = {52,5,53};


Point(6) = { 0, RevH/2, 0, lc2 } ;
Point(61) = { 0 , RevH/2 - radius, 0, lc2 } ;
Point(62) = { radius , RevH/2 , 0, lc2} ;
Point(63) = { 0 , RevH/2 + radius , 0, lc2} ;
Circle(601) = {61,6,62};
Circle(602) = {62,6,63};


Point(7) = { RevW/2, RevH, 0, lc2 } ;
Point(71) = { RevW/2 - radius, RevH, 0, lc2 } ;
Point(72) = { RevW/2 , RevH - radius , 0, lc2} ;
Point(73) = { RevW/2 + radius , RevH , 0, lc2} ;
Circle(701) = {71,7,72};
Circle(702) = {72,7,73};


Point(8) = { RevW, RevH/2, 0, lc2 } ;
Point(81) = { RevW, RevH/2 + radius, 0, lc2 } ;
Point(82) = { RevW - radius , RevH/2 , 0, lc2} ;
Point(83) = { RevW, RevH/2 - radius  , 0, lc2} ;
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



//Plane Surface(1) = {1};

Line Loop(2) = { 801,802,6,7};
Line Loop(3) = { 701,702,10,11};
Line Loop(4) = { 601,602,14,15};
Line Loop(5) = { 501,502,2,3};

Line Loop(1) = {4,5,-802,-801,8,9,-702,-701,12,13,-602,-601,16,1,-502,-501};
Plane Surface(1) = {1};


Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

out1[]=Extrude{0,0,length}{Surface{2,3,4,5};};

Point(200) = { 0, 0, -RevD, lc2 };
Point(201) = { 0,RevH,-RevD, lc2 } ;
Point(202) = { RevW, RevH, -RevD, lc2 } ;
Point(203) = { RevW, 0, -RevD, lc2 } ;
Line(201) = {200,201};
Line(202) = {201,202};
Line(203) = {202,203};
Line(204) = {203,200};
Line(205) = {200,1};
Line(206) = {201,2};
Line(207) = {202,3};
Line(208) = {203,4};
Line Loop(201) = {201,202,203,204};
Line Loop(202) = {205,-16,-15,-14,-13,-206,-201};
Line Loop(203) = {206,-12,-11,-10,-9,-207,-202};
Line Loop(204) = {207,-8,-7,-6,-5,-208,-203};
Line Loop(205) = {208,-4,-3,-2,-1,-205,-204};
//Line Loop(203) = {201,202,203,-204};
//Line Loop(204) = {201,202,203,-204};
Plane Surface(201) = {201};
Plane Surface(202) = {202};
Plane Surface(203) = {203};
Plane Surface(204) = {204};
Plane Surface(205) = {205};
Surface Loop(2) = {1,2,3,4,5,201,202,203,204,205};
Volume(10) = {2};

NewSurf1[]=Translate{0,0,length}{ Duplicata{Surface{1};}};
Point(300) = { 0, 0, RevD+length, lc2 };
Point(301) = { 0,RevH,RevD+length, lc2 } ;
Point(302) = { RevW, RevH, RevD+length, lc2 } ;
Point(303) = { RevW, 0, RevD+length, lc2 } ;
Line(301) = {300,301};
Line(302) = {301,302};
Line(303) = {302,303};
Line(304) = {303,300};
Line(305) = {300,259};
Line(306) = {301,241};
Line(307) = {302,223};
Line(308) = {303,205};
Line Loop(301) = {301,302,303,304};
Line Loop(302) = {305,-904,-851,-850,-901,-306,-301};
Line Loop(303) = {306,-900,-829,-828,-897,-307,-302};
Line Loop(304) = {307,-896,-807,-806,-893,-308,-303};
Line Loop(305) = {308,-892,-873,-872,-905,-305,-304};

//Line Loop(204) = {201,202,203,-204};
Plane Surface(301) = {301};
Plane Surface(302) = {302};
Plane Surface(303) = {303};
Plane Surface(304) = {304};
Plane Surface(305) = {305};
Surface Loop(3) = {891,846,868,824,890,301,302,303,304,305};
Volume(20) = {3};


Compound Volume(22) = {1,2,3,4,10,20};
//No mesh refinement yet, need to do mesh refinement 

// Define Fields.... ?

Field[4] = Attractor;
//Field[4].EdgesList = {926,948,882,904};
Field[4].EdgesList = {840,818,862,884};


// We can then create a MathEval field with a function that depends on the
// return value of the Attractr Field[4], i.e., depending on the distance to
// point 1 (here using a cubic law, with minumum element size = lc / 100)
Field[5] = MathEval;
Field[5].F = Sprintf("2e-10*(5- F4/1e-9)^4 + %g", 4e-9 / 14);


Field[6] = Min;
Field[6].FieldsList = {5};
Background Field = 6;



//simple 2D 40 x 40nm square to validate against COMSOL

// parameter 
 
length = 40e-9 ; //length of square domain



lc = 0.5e-9 ;
lc2 = 0.4e-9;

Point(1) = { 0, 0, 0, lc };
Point(2) = { 0,length,0, lc } ;
Point(3) = { length,length,0, lc } ;
Point(4) = { length,0,0, lc } ;

// define lines


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};



//
//No mesh refinement yet, need to do mesh refinement 



Field[1] = Attractor;
Field[1].EdgesList = {4};
//Field[1].EdgesList = {31,32,36,40};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMax = lc2;
Field[2].LcMin = lc2/8;
Field[2].DistMax = 7e-9;
Field[2].DistMin = 5e-9;



Field[6] = Min;
Field[6].FieldsList = {2};
Background Field = 6;

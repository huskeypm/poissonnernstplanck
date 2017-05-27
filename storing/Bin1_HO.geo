lc = 50.000000;

Point(1) = {9.734513, 23.008850, 0, lc};
Line(1) = {1,1};
Line Loop(2) = {1};

Point(2) = {7.079646, 24.778761, 0, lc};
Line(3) = {2,2};
Line Loop(4) = {3};

Point(3) = {8.849558, 29.203540, 0, lc};
Line(5) = {3,3};
Line Loop(6) = {5};

Point(4) = {13.274336, 28.318584, 0, lc};
Line(7) = {4,4};
Line Loop(8) = {7};

Point(5) = {14.159292, 24.778761, 0, lc};
Line(9) = {5,5};
Line Loop(10) = {9};

Point(6) = {10.619469, 23.008850, 0, lc};
Line(11) = {6,6};
Line Loop(12) = {11};
// Container
Point(9) = {0.000000, 0.000000, 0, lc};
Point(10) = {0.000000, 35.398230, 0, lc};
Point(11) = {35.398230, 35.398230, 0, lc};
Point(12) = {35.398230, 0.000000, 0, lc};
Line(15)={9,10};
Line(16)={10,11};
Line(17)={11,12};
Line(18)={12,9};
Line Loop(19) ={15,16,17,18};
Plane Surface(20) ={2,4,6,8,10,12,19};

// Create an 'attractor' around these lines to refine mesh in this area
// Need to play around with the DistMin/DistMax values to get something 
// reasonable 
Field[4] = Attractor;
Field[4].EdgesList = {1,3,5,7,9,11};
Field[5] = Threshold;
Field[5].IField = 4;
Field[5].LcMin = lc / 5;
Field[5].LcMax = lc;
Field[5].DistMin = 30;
Field[5].DistMax = 45;

// Using a cubic law, with minumum element size = lhigh
// In practice, this doesnt work too well 
//Field[5] = MathEval;
//lrate =0.001;
//Field[5].F = Sprintf("(F4*%g)^3 + %g", 0.02,lc/5);


// Use the minimum of all the fields as the background mesh field
Field[7] = Min;
Field[7].FieldsList = {5};
Background Field = 7;

Physical Line(111) = {1,3,5,7,9,11};
Physical Surface(3)={20};
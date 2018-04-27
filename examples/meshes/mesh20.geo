lc = 0.8e-9;

length = 100e-9;
width = 50e-9;
gap = 20e-9;

Point(1) = {0, 0, 0, lc};
Point(2) = {0.25*length, 0, 0, lc};
Point(3) = {0.25*length, 0.5*(width - gap), 0, lc};
Point(4) = {0.75*length, 0.5*(width - gap), 0, lc};
Point(5) = {0.75*length, 0, 0, lc};
Point(6) = {length, 0, 0, lc};
Point(7) = {length, width, 0, lc};
Point(8) = {0.75*length, width, 0, lc};
Point(9) = {0.75*length, 0.5*(width + gap), 0, lc};
Point(10) = {0.25*length, 0.5*(width + gap), 0, lc};
Point(11) = {0.25*length, width, 0, lc};
Point(12) = {0, width, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(1) = {1};


Field[1] = Attractor;
Field[1].EdgesList = {3, 9};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMax = lc;
Field[2].LcMin = lc/8;
Field[2].DistMax = 1e-9;
Field[2].DistMin = 5e-10;

Field[6] = Min;
Field[6].FieldsList = {2};
Background Field = 6;

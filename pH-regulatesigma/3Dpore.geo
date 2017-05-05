//build a single 3D nanopore 

// parameter for reservoir
RevL = 20e-9 ;
RevH = 20e-9 ;



//parameter for nanopore

length = 30e-9 ;
radius = 4e-9 ;

lc1 =1e-9 ;
lc2 = 0.2e-9 ;


Point(1) = { 0, 0, 0, lc1 };
Point(2) = { 0,RevL,0, lc1 } ;
Point(3) = { RevL, RevL, 0, lc1 } ;
Point(4) = { RevL, 0, 0, lc1 } ;

Point(5) = { RevL/2, RevL/2, 0, 1 } ;
Point(6) = { RevL/2 , RevL/2 - radius , 0, lc2} ;
Point(7) = { RevL/2 - radius , RevL/2, 0, lc2 } ;
Point(8) = { RevL/2 , RevL/2 + radius, 0, lc2 } ;
Point(9) = { RevL/2 + radius , RevL/2, 0, lc2} ;



Line(1) = { 1,2};
Line(2) = { 2,3} ;
Line(3) = { 3,4} ;
Line(4) = { 4,1} ;

Circle(5)={6,5,7} ;
Circle(6)={7,5,8} ;
Circle(7)={8,5,9} ;
Circle(8)={9,5,6} ;


Line Loop(1) = {5,6,7,8} ;
Line Loop(2) = {1,2,3,4} ;

Plane Surface(1) = {1} ;
Plane Surface(2) = {2} ;
//Physical Surface(2) = {1,2};

Extrude{0,0,length}{Surface{1};Layers{40};} ;
Extrude{0,0,-RevH}{Surface{2};Layers{40};} ;
NewSurf[]=Translate{0,0,length}{ Duplicata{Surface{2};}} ;
Extrude{0,0,RevH}{Surface{NewSurf[0]};Layers{20};} ;


length = 80e-9 ;
interdis = 20e-9 ;
radius = 10e-9 ;

lc1 =2e-9 ;
lc2 = 1e-10 ;


Point(1) = { 0, 0, 0, lc1 };
Point(2) = { 0,length,0, lc1 } ;
Point(3) = { length, length, 0, lc1 } ;
Point(4) = { length, 0, 0, lc1 } ;

Point(5) = { length/2 - interdis/2 - radius , length/2, 0, 1 } ;
Point(6) = { length/2 - interdis/2 - radius , length/2 - radius , 0, lc2} ;
Point(7) = { length/2 - interdis/2 - 2*radius , length/2, 0, lc2 } ;
Point(8) = { length/2 - interdis/2 - radius , length/2 + radius, 0, lc2 } ;
Point(9) = { length/2 - interdis/2 , length/2, 0, lc2} ;


Point(10) = { length/2 + interdis/2 + radius , length/2, 0, lc2 } ;
Point(11) = { length/2 + interdis/2 + radius , length/2 - radius, 0, lc2 } ;
Point(12) = { length/2 + interdis/2  , length/2, 0, lc2 } ;
Point(13) = { length/2 + interdis/2 + radius , length/2 + radius, 0, lc2} ;
Point(14) = { length/2 + interdis/2 + 2*radius , length/2, 0, lc2 } ;

Line(1) = { 1,2};
Line(2) = { 2,3} ;
Line(3) = { 3,4} ;
Line(4) = { 4,1} ;

Circle(5)={6,5,7} ;
Circle(6)={7,5,8} ;
Circle(7)={8,5,9} ;
Circle(8)={9,5,6} ;

Circle(9) = {11,10,12} ;
Circle(10) = {12,10,13} ;
Circle(11) = {13,10,14} ;
Circle(12) = {14,10,11} ;

Line Loop(1) = {5,6,7,8} ;
Line Loop(2) = {1,2,3,4} ;
Line Loop(3) = {9,10,11,12} ;

Plane Surface(1) = {1} ;
Plane Surface(2) = {2,1,3} ;

Physical Surface(1) = {2};

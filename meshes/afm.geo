/*
Gmsh model of a nanocyrstalline silicon film.

2007 Ondrej Certik
*/

h=1.0;
r=2.0;

Point(1) = {0, h, 0, 0.1};
Point(2) = {0, 0, 0, 0.1};
Point(3) = {r, 0, 0, 0.1};
Point(4) = {r, h, 0, 0.1};
Point(5) = {0.5*r/2, h, 0, 0.1};
Point(6) = {0.4*r/2, 0.8*h, 0, 0.1};
Point(7) = {0.3*r/2, 0.6*h, 0, 0.1};
Point(8) = {0.2*r/2, 0.4*h, 0, 0.1};
Point(9) = {0.1*r/2, 0.2*h, 0, 0.1};
Point(10) = {0, 0.1*h, 0, 0.1};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,1};
Line(5) = {1,10};
Line(6) = {10,2};
Line(7) = {10,9};
Line(8) = {9,8};
Line(9) = {8,7};
Line(10) = {7,6};
Line(11) = {6,5};
Line Loop(12) = {3,-11,-10,-9,-8,-7,6,1,2};
Plane Surface(13) = {12};
Line Loop(14) = {4,5,7,8,9,10,11};
Plane Surface(15) = {14};
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{13,15};
}
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{57,89};
}
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{131,163};
}
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{237,205};
}
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{311,269};
}
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{385,353};
}
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{417,459};
}
Extrude {{0,1,0}, {-0.9+0.9,0,0}, Pi/4} {
  Surface{491,533};
}
Surface Loop(606) = {177,103,29,578,505,431,325,283,-268,181,107,33,-564,-490,-416,329,333,-264,185,111,37,-560,-486,-412,-408,337,-260,189,115,41,-556,-482,-478,-404,341,-256,193,119,45,-552,-548,-474,-400,344,-252,196,122,48,310,204,130,56,605,532,458,352,348,306,200,126,52,601,528,454};
Volume(607) = {606};
Surface Loop(608) = {329,-268,181,107,33,-564,-490,-416,-396,-470,-544,-68,-142,-216,-248,-364,-412,333,-264,185,111,37,-560,-486,-482,-408,337,-260,189,115,41,-556,-552,-478,-404,341,-256,193,119,45,48,122,196,-252,344,-400,-474,-548};
Volume(609) = {608};

//zbytek
Physical Volume (100)={607};
//zrno
Physical Volume (101)={609};

//dole - substrat
Physical Surface(2) = {528,454,348,306,200,126,52,601};
//nahore klobouk zrna
Physical Surface(1) = {216,248,142,68,544,470,396,364};
//nahore - zbytek kolem zrna
Physical Surface(3) = {103,177,283,325,431,29,578,505};
//po stranach - obvod
Physical Surface(4) = {605,532,458,352,310,204,130,56};
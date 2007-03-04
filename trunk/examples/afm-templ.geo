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

Line(534) = {303,4};
Line(535) = {3,332};

//rt=0.5;
rt=%f;

x0=1.1;
y0=rt;
eps=0.3;

Point(333) = {x0,h,y0,0.1};
Point(334) = {x0+eps,h,y0,0.1};
Point(335) = {x0-eps,h,y0,0.1};
Point(336) = {x0,h,y0+eps,0.1};
Point(337) = {x0,h,y0-eps,0.1};
Circle(536) = {334,333,337};
Circle(537) = {337,333,335};
Circle(538) = {335,333,336};
Circle(539) = {336,333,334};
Line Loop(540) = {537,538,539,536};
Plane Surface(541) = {540};
Extrude {{0,1,0}, {0,0,0}, Pi/4} {
  Surface{491};
}
Line Loop(573) = {3,-551,-493,534};
Plane Surface(574) = {573,540};
Line Loop(575) = {500,-535,-1};
Plane Surface(576) = {575};
Line Loop(577) = {501,534,-2,535};
Plane Surface(578) = {577};
Surface Loop(579) = {578,-533,-574,13,-572,-568,-564,-560,-556,576,-541};
Volume(580) = {579};

//substrate:
Physical Surface(1) = {126,200,306,348,454,52,576,528};

//tip:
Physical Surface(2) = {541};

//amorphous:
Physical Volume(100) = {1,3,5,8,9,12,14,580};

//grain:
Physical Volume(101) = {2,4,6,7,10,11,13,15};

//sides:
Physical Surface(3) = {578,532,458,352,310,204,130,56};

//top, without tip:
Physical Surface(4) = {431,325,283,177,505,103,29,574,396,364,248,216,470,552,68,142};

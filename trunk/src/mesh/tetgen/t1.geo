R=10; //radius
n=20; //number of nodes on the radius.

lt=R/n;
phi=Pi/2;

Point(1) = {0,0,0,lt};

Point(2) = {R,0,0,lt};
Line(1) = {1,2};

Extrude Line {1, {0.0,0.0,1.0}, {0.0,0.0,0.0}, phi};
Extrude Line {2, {0.0,0.0,1.0}, {0.0,0.0,0.0}, phi};

Extrude Surface {4, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {7, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};

Extrude Surface {19, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {31, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {43, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {55, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {67, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};
Extrude Surface {79, {1.0,0.0,0.0}, {0.0,0.0,0.0}, phi};

//Transfinite Line {1,2,5,-11,-35,-59} = n Using Progression 1.05;

//Physical Volume(100) = {1,2,3,4,5,6,7,8};
Physical Surface(1) = {100,29,77,53,39,87,63,15};

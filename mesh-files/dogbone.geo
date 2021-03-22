Mesh.MshFileVersion = 2.2;

l = 0.5;
h = 0.5;
r = 0.25;

E = 0.05;
e = 0.01;

Point(1) = { l, -h, 0, E};
Point(2) = { l,  h, 0, E};
Point(3) = { r,  h, 0, E};
Point(4) = { 0,  r, 0, e};
Point(5) = {-r,  h, 0, E};
Point(6) = {-l,  h, 0, E};
Point(7) = {-l, -h, 0, E};
Point(8) = {-r, -h, 0, E};
Point(9) = { 0, -r, 0, e};
Point(10) = {r, -h, 0, E};

Point(11) = {0,  h, 0, e};
Point(12) = {0, -h, 0, e};

Line(1) = {1, 2};
Line(2) = {2, 3};
Circle(3) = {4, 11, 3};
Circle(4) = {5, 11, 4};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Circle(8) = {9, 12, 8};
Circle(9) = {10, 12, 9};
Line(10) = {10, 1};
Line(11) = {4, 9};

Line Loop(1) = {1, 2, -3, -4, 5, 6, 7, -8, -9, 10};
Plane Surface(1) = {1};
Line{11} In Surface {1};

Physical Surface("all") = {1};
Physical Line("left") = {6};
Physical Line("right") = {1};

//Recombine Surface "*";

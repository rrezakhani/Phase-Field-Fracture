Mesh.MshFileVersion = 2.2;

Nb_elems = 1;

x_size = 1;
y_size = 1;

h_e = x_size/Nb_elems;

Point(1) = {0,           0, 0, h_e};
Point(2) = {x_size,      0, 0, h_e};
Point(3) = {x_size, y_size, 0, h_e};
Point(4) = {0,      y_size, 0, h_e};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Line("top") = {3};
Physical Line("bottom") = {1};
Physical Line("left") = {4};
Physical Line("right") = {2};

Physical Surface("block") = {1};

Transfinite Surface "*";
Recombine Surface "*";

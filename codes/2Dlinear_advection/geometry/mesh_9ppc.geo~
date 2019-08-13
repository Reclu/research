lc = 0.1;
length = 4.;
height = 4.;

MPoints=12;
dx=4./(MPoints-1);
dy=4./(MPoints-1);

NNodes=7;

Point(1) = {-dx/2., -dy/2., 0, lc};
Point(2) = {length+dx/2., -dy/2., 0, lc};
Point(3) = {length+dx/2., height+dy/2., 0, lc};
Point(4) = {-dx/2.,height+dy/2. , 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(6) = {1, 2, 3,4};
Plane Surface(6) = {6};

// Numer of nodes on the line
Transfinite Line {1,2,3,4} = 7;

Transfinite Surface "*";
Recombine Surface "*";

lc = 0.1;
dx=0.125;
dy=0.125;
length = 4.;
height = 3.;

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
Transfinite Line {1} = 34;
Transfinite Line {3} = 34;
Transfinite Line {2} = 26;
Transfinite Line {4} = 26;

Transfinite Surface "*";
Recombine Surface "*";

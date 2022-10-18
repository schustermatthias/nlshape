delta = 0.1;
cl__1 = 1;
lc = 1;
lc_fine = lc/2;
lc_point= lc/(15);

Point(999) = {0.5, 0.5, 0.0, lc_point};


//----------------------------------------
// Omega
Point(1) = {0, 0, 0, 1};
Point(2) = {0, 1, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {1, 0, 0, 1};

// for circle
Point(5) = {0.25, 0.5, 0, 1};
Point(6) = {0.5, 0.5, 0, 1};

// Omega_I
Point(8) = {1. + delta, -delta, -0, lc};
Point(9) = {-delta, -delta, -0, lc};
Point(10) = {1 + delta, 1 + delta, -0, lc};
Point(11) = {-delta, 1 + delta, -0, lc};

Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};

Circle(5) = {5, 6, 5};

Line(6) = {11, 9};
Line(7) = {9, 8};
Line(8) = {8, 10};
Line(9) = {10, 11};

Line Loop(15) = {6, 7, 8, 9, -4, -3, -2, -1};
Plane Surface(15) = {15};

Line Loop(17) = {1, 2, 3, 4, -5};
Plane Surface(17) = {17};

Line Loop(18) = {5};
Plane Surface(18) = {18};

Physical Line(9) = {1, 2, 3, 4};
Physical Line(12) = {5};
Physical Line(13) = {6, 7, 8, 9};

Physical Surface(3) = {15};
Physical Surface(2) = {17};
Physical Surface(1) = {18};


// For coarsening mesh around midpoint
Point {999} In Surface {18};


// INTERFACE
Field[1] = Distance;
Field[1].EdgesList = {5};
Field[1].NNodesByEdge = 5000;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_fine;// element size inside DistMin
Field[2].LcMax = lc;  // element size outside DistMax
Field[2].DistMin = 0.08;
Field[2].DistMax = 0.1;


// Define minimum of threshold and function field
Field[5] = Min;
Field[5].FieldsList = {2};


// Use the min as the background field
Background Field = 5;


delta = 0.1;
cl__1 = 1;
lc = 1;
lc2 = 1;
lc3 = 1;
lc_fine = lc/2;
lc_point= lc/(15);


// SHAPE
Point(12) = {0.27, 0.5, 0, lc};
Point(13) = {0.73, 0.5, 0, lc};
Point(15) = {0.5, 0.27, 0, lc};
Point(16) = {0.5, 0.73, 0, lc};

// Point for coarsening (in the center of the shape)
Point(999) = {0.5, 0.5, 0.0, lc_point};


//------------------------------------------------------------------------------
// OMEGA
Point(1) = {0, 0, 0,  lc};
Point(2) = {0, 1, 0,  lc};
Point(3) = {1, 1, 0,  lc};
Point(4) = {1, 0, 0,  lc};

// Omega_I
Point(8) = {1. + delta, -delta, -0, lc};
Point(9) = {-delta, -delta, -0, lc};
Point(10) = {1 + delta, 1 + delta, -0, lc};
Point(11) = {-delta, 1 + delta, -0, lc};



Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line(6) = {11, 9};
Line(7) = {9, 8};
Line(8) = {8, 10};
Line(9) = {10, 11};

Spline(10) = {12, 15, 13, 16, 12};

Line Loop(14) = {10};
Plane Surface(15) = {14};
Line Loop(16) = {4, 1, 2, 3};
Plane Surface(17) = {14, 16};
Line Loop(18) = {9, 6, 7, 8};
Plane Surface(19) = {16, 18};

//=============== LABELING ===============//
// Interface
Physical Line(12) = {10};
Physical Line(9) = {1, 2, 3, 4};
Physical Line(13) = {6, 7, 8, 9};

// Omega_(...)
Physical Surface(1) = {15};
Physical Surface(2) = {17};
Physical Surface(3) = {19};
 


// SHAPE ATTRACTORS

// For coarsening mesh around midpoint
Point{999} In Surface {15};

// INTERFACE
Field[1] = Distance;
Field[1].EdgesList = {10};
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




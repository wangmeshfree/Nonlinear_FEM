SetFactory("OpenCASCADE");

// ------------------------------------------------------------
// Parameters: L is the beam length and R is the beam height.
// The beam is centred about y = 0: -R/2 <= y <= R/2.
// ------------------------------------------------------------
L = 10;
R = 2;
mesh_size = 0.25;

// ------------------------------------------------------------
// Rectangle corners
// ------------------------------------------------------------
Point(1) = {0, -R/2, 0, mesh_size};
Point(2) = {L, -R/2, 0, mesh_size};
Point(3) = {L,  R/2, 0, mesh_size};
Point(4) = {0,  R/2, 0, mesh_size};

// Boundary lines: bottom, right, top, and fixed left edge.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups are read by read_lsdyna_k() in the notebook.
Physical Surface("BEAM") = {1};
Physical Line("FREE_BOTTOM") = {1};
Physical Line("TRACTION_RIGHT") = {2};
Physical Line("FREE_TOP") = {3};
Physical Line("FIXED_LEFT") = {4};

// Generate a quad-dominant 2D mesh for the Q4 element code.
Mesh.CharacteristicLengthMin = mesh_size;
Mesh.CharacteristicLengthMax = mesh_size;
Mesh.RecombineAll = 1;
Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 1;
Mesh 2;

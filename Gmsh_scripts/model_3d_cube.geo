SetFactory("OpenCASCADE");

// ------------------------------------------------------------
// Cube of size Lx x Ly x Lz meshed with n elements per direction
// (structured 8-node hexahedra for the H8 code in src/FEM3D).
// ------------------------------------------------------------
Lx = 1;
Ly = 1;
Lz = 1;
n  = 4;

Box(1) = {0, 0, 0, Lx, Ly, Lz};

// OpenCASCADE box surface tags: 1:x=0, 2:x=Lx, 3:y=0, 4:y=Ly, 5:z=0, 6:z=Lz
Physical Volume("CUBE") = {1};
Physical Surface("XMIN") = {1};
Physical Surface("XMAX") = {2};


// Structured hex mesh: n elements per edge.
Transfinite Curve {:} = n + 1;
Transfinite Surface {:};
Recombine Surface {:};
Transfinite Volume {:};

Mesh.ElementOrder = 1;
Mesh 3;

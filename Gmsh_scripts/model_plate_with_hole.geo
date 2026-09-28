SetFactory("OpenCASCADE");

// --------------------
// Parameters
// --------------------
L = 5;
R = 1;

// --------------------
// Four blocks
// --------------------
Rectangle(1) = {0,   0,   0, L/2, L/2, 0};
Rectangle(2) = {L/2, 0,   0, L/2, L/2, 0};
Rectangle(3) = {L/2, L/2, 0, L/2, L/2, 0};
Rectangle(4) = {0,   L/2, 0, L/2, L/2, 0};

// --------------------
// Hole cutter
// --------------------
Disk(5) = {0, 0, 0, R, R};

// --------------------
// Cut hole from lower-left block
// --------------------
s1[] = BooleanDifference{
    Surface{1}; Delete;
}{
    Surface{5}; Delete;
};

// --------------------
// Make the four regions conforming
// --------------------
all[] = BooleanFragments{
    Surface{s1[], 2, 3, 4}; Delete;
}{};

// --------------------
// Mesh size
// --------------------
Mesh.CharacteristicLengthMin = 0.3;
Mesh.CharacteristicLengthMax = 0.3;

// --------------------
// Quad-dominant mesh
// --------------------
Mesh.RecombineAll = 1;
Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 1;

// --------------------
// Physical group
// --------------------
Physical Surface("PLATE") = {all[]};
Physical Line("FIXED_LEFT") = {2,13 };
Physical Line("FIXED_BOTTOM") = {5,6};
Physical Line("TRACTION_RIGHT") = {7,9};
Physical Line("TRACTION_TOP") = {10,12};

Mesh 2;
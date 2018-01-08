

Include "Sphere.geo";
Include "Cube.geo";


////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
// Geometry dimensions
///////////////////////////////////////////////////////////////

outer_r = 1.0; // outer r
inner_r = 0.25; // inner r
n1 = 3;
n2 = 3;

////////////////////////////////////////////////////////////////
// Type of elements
////////////////////////////////////////////////////////////////

NonLinearQ = 0;
IsTetraQ =  1;
IsPrismQ = 0;

//Call MakeSphere;

Call MakeCube;

Coherence Mesh;

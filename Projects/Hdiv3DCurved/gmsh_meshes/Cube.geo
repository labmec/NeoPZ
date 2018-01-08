
// ---- Gmsh Macro ----
// ---- Box egion  ----
// Created 05/01/2018 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Macro MakeCube

If (NonLinearQ == 1)
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
EndIf

l = 2.0;
x_length = 1.0;
y_length = 1.0;
z_length = 1.0;


p1 = newp; Point(p1) = { 0.0, 0.0, 0.0, l};
p2 = newp; Point(p2) = { x_length, 0.0, 0.0, l};
p3 = newp; Point(p3) = { x_length,  y_length, 0.0, l};
p4 = newp; Point(p4) = { 0.0,  y_length, 0.0, l};

p5 = newp; Point(p5) = {0.0, 0.0, z_length, l};
p6 = newp; Point(p6) = { x_length, 0.0, z_length, l};
p7 = newp; Point(p7) = { x_length,  y_length, z_length, l};
p8 = newp; Point(p8) = { 0.0,  y_length, z_length, l};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

l5 = newl; Line(l5) = {p5,p6};
l6 = newl; Line(l6) = {p6,p7};
l7 = newl; Line(l7) = {p7,p8};
l8 = newl; Line(l8) = {p8,p5};

l9  = newl; Line(l9)  = {p5,p1};
l10 = newl; Line(l10) = {p6,p2};
l11 = newl; Line(l11) = {p7,p3};
l12 = newl; Line(l12) = {p8,p4};

ll1  = newll; Line Loop(ll1) = {l1,  l2,   l3, l4}; // Bottom
ll2  = newll; Line Loop(ll2) = {l5,  l6,   l7, l8}; // Top
ll3  = newll; Line Loop(ll3) = {l1, -l10, -l5, l9}; // South
ll4  = newll; Line Loop(ll4) = {l2, -l11, -l6, l10}; // East
ll5  = newll; Line Loop(ll5) = {l3, -l12, -l7, l11}; // North
ll6  = newll; Line Loop(ll6) = {l4, -l9,  -l8, l12}; // West

s1  = news; Plane Surface(s1) = {ll1}; // Bottom unstructured region
s2  = news; Plane Surface(s2) = {ll2}; // Top unstructured region
s3  = news; Plane Surface(s3) = {ll3}; // South unstructured region
s4  = news; Plane Surface(s4) = {ll4}; // East unstructured region
s5  = news; Plane Surface(s5) = {ll5}; // North unstructured region
s6  = news; Plane Surface(s6) = {ll6}; // West unstructured region

res_B[] = {s1};
res_T[] = {s2};
res_S[] = {s3};
res_E[] = {s4};
res_N[] = {s5};
res_W[] = {s6};

lateral_boundaries[] = {res_S[],res_E[],res_N[],res_W[]};
top_bott_boundaries[] = {res_B[],res_T[]};


// reservoir
sl1 = newsl; Surface Loop(sl1) = {lateral_boundaries[],top_bott_boundaries[]};
v1  = newv; Volume(v1) = {sl1} ;
cube[] = {v1};



If(IsTetraQ)


Else

Transfinite Surface "*";
Transfinite Volume "*";

EndIf


// Tagging physical entities
Physical Volume("volume") = {cube[]};
Physical Surface("inner") = {lateral_boundaries[]};
Physical Surface("outer") = {top_bott_boundaries[]};


Return

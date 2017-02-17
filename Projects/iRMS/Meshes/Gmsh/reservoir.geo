// ---- Cylinder wellbore Region Gmsh scritp ----
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

// Settings
ExpertMode = 1;


Include "CadReservoir.geo";
Include "BoxReservoir.geo";
Include "BoxSideBurden.geo";
Include "BuildOmegas.geo";
Include "drill_well.geo";
Include "drill_well_box.geo";
Include "PhysicalEntities.geo";

well_index = 0;
well_lids = {};

well_p_bores = {};
well_p_regions = {};
well_p_v_regions = {};

well_i_bores = {};
well_i_regions = {};
well_i_v_regions = {};


geomechanicQ = 0;
dimension = 2;
nolinearQ = 0;
CADReservoirQ = 0;

xzQ = 0;
hexahedronsWQ = 0;
hexahedronsRQ = 0;
hexahedronsSBQ = 0;

If (nolinearQ == 1)
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
EndIf

If (hexahedronsWQ == 1 || hexahedronsRQ == 1 || hexahedronsSBQ == 1)
Mesh.Algorithm3D = 6 ;
Else
Mesh.Algorithm3D = 1 ;
EndIf


// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 1;
cl2 = 0.1;
cl3 = 10.0;
cl4 = 250.0;
cl5 = 5000.0;

////////////////////////////////////////////////////////////////////////////
// reservoir region geometry
////////////////////////////////////////////////////////////////////////////

// reservoir box dimensions
x_length = 1000.0;
y_length = 1000.0;
z_length = 100.0;

////////////////////////////////////////////////////////////////////////////
// side-burden region geometry
////////////////////////////////////////////////////////////////////////////

// side-burden box dimensions
sb_x_length = 10000.0;
sb_y_length = 10000.0;
sb_z_length = 4000.0;

////////////////////////////////////////////////////////////////////////////
// reservoir rock
////////////////////////////////////////////////////////////////////////////
If(CADReservoirQ == 1)
Call ReservoirCAD;
Else
Call ReservoirBox;
EndIf

////////////////////////////////////////////////////////////////////////////
// side-burden rock
////////////////////////////////////////////////////////////////////////////
If (geomechanicQ == 1)
Call SideBurdenBox;
EndIf


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Drilling wells 
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
// well geometry and settings
////////////////////////////////////////////////////////////////////////////

rw = 0.5;
wl = 5.0;

wbr = 10.0;
ela = 20.0;
rw_cell= 1.0;
wr_cell= 20.0;

// Orientation and length
alfa = 0.0*Pi/2.0;
beta = 0.0;

////////////////////////////////////////////////////////////////////////////
// well 1 
////////////////////////////////////////////////////////////////////////////

// well location
//wcx = 0.0;
//wcy = 0.0;
//wcz = 0.0;

wcx = -200.0;
wcy = 100.0;
wcz = 150.0;

wcx = -138;
wcy = 100.0;
wcz = 0.0;

IsInjectorQ = 0;
//Call DrillWell;


////////////////////////////////////////////////////////////////////////////
// well 2 
////////////////////////////////////////////////////////////////////////////

// well location
//wcx = 400.0;
//wcy = 400.0;
//wcz = 0.0;

wcx = 150.0;
wcy = 0.0;
wcz = 40.0;
IsInjectorQ = 1;
//Call DrillWell;


////////////////////////////////////////////////////////////////////////////
// well 3 
////////////////////////////////////////////////////////////////////////////

// well location
//wcx = -400.0;
//wcy = -400.0;
//wcz = 0.0;

wcx = -650.0;
wcy = 500.0;
wcz = 20.0;
IsInjectorQ = 1;
//Call DrillWell;

////////////////////////////////////////////////////////////////////////////
// well 4 
////////////////////////////////////////////////////////////////////////////

// well location
wcx = -650.0;
wcy = -500.0;
wcz = 20.0;
IsInjectorQ = 1;
//Call DrillWell;


////////////////////////////////////////////////////////////////////////////
// well 5 
////////////////////////////////////////////////////////////////////////////

// well location
wcx = -150.0;
wcy = -300.0;
wcz = 120.0;
IsInjectorQ = 0;
//Call DrillWell;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Computational domain, boundaries:: Tagging physical entities
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

Call DefineOmegas;
Call DrawBoundaries;

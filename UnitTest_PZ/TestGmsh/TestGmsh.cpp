/**
 * @file
 * @brief Contains Unit Tests for reading .msh files generated from gmsh
 */

#include "TPZGmshReader.h"
#include "pzcheckgeom.h"

#include <catch2/catch_test_macros.hpp>


TEST_CASE("Periodic meshes", "[gmsh]"){

  /*simple periodic mesh on a unit square,
    testing multiple periodicity*/
  SECTION("Multiple periodicities 2D"){
    TPZGmshReader meshReader;
    TPZAutoPointer<TPZGeoMesh> gmesh =
      meshReader.GeometricGmshMesh("multiple_periodic_2d.msh");
    REQUIRE(true);


    TPZCheckGeom check(gmesh.operator->());
    REQUIRE(check.CheckIds()==0);
  }

  /**
     the mesh is a rectangle, periodic on the x direction
     given the larger element sizes, we have non-periodic elements
     with vertices on both periodic edges

     i've found a bug in a similar situation, so it is important
     to test (the problematic element was not periodic, so it was hard to find it)
   */
  SECTION("Periodic layer"){
    TPZGmshReader meshReader;
    TPZAutoPointer<TPZGeoMesh> gmesh =
      meshReader.GeometricGmshMesh("periodic_layer.msh");
    REQUIRE(true);


    TPZCheckGeom check(gmesh.operator->());
    REQUIRE(check.CheckIds()==0);
  }
}

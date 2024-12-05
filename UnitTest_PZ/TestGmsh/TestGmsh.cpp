/**
 * @file
 * @brief Contains Unit Tests for reading .msh files generated from gmsh
 */

#include "TPZGmshReader.h"
#include "TPZGeoMeshTools.h"
#include "pzcheckgeom.h"

#include <catch2/catch_test_macros.hpp>
#include <numeric> //for iota

///check node ordering of every periodic element in mesh
void TestPeriodicElements(TPZGeoMesh *gmesh, const SPZPeriodicData& periodic_data);

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

    auto data = meshReader.GetPeriodicData();
    TestPeriodicElements(gmesh.operator->(), *data);
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

    auto data = meshReader.GetPeriodicData();
    TestPeriodicElements(gmesh.operator->(), *data);
  }
}


void TestPeriodicElements(TPZGeoMesh *gmesh, const SPZPeriodicData& periodic_data){

  auto &data = periodic_data;

  TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> el_map_vec;
  TPZGeoMeshTools::FindPeriodicElements(gmesh, data.dep_mat_ids, data.indep_mat_ids, data.nodes_map,
                                        el_map_vec);

  for(auto el_map : el_map_vec){
    for(auto [dep_idx, indep_idx] : *el_map){
      auto dep_el = gmesh->Element(dep_idx);
      auto indep_el = gmesh->Element(indep_idx);
      REQUIRE(dep_el->NCornerNodes() == indep_el->NCornerNodes());
      const auto nnodes = dep_el->NCornerNodes();
      constexpr int max_nnodes{8};
      TPZManVector<int64_t,max_nnodes> dep_ids(nnodes,0),indep_ids(nnodes,0);;
      for(auto i=0; i<nnodes; i++) {
        dep_ids[i] = dep_el->NodePtr(i)->Id();
        indep_ids[i] = indep_el->NodePtr(i)->Id();
      }

      //now we compare their ordering
      TPZManVector<int64_t,max_nnodes> dep_order(nnodes,0), indep_order(nnodes,0);
      std::iota(dep_order.begin(), dep_order.end(), 0);
      std::iota(indep_order.begin(), indep_order.end(), 0);

      // sort indexes based on comparing values in v
      // using std::stable_sort instead of std::sort
      // to avoid unnecessary index re-orderings
      // when v contains elements of equal values 
      std::sort(dep_order.begin(), dep_order.end(),
                [&dep_ids](const int64_t i1, const int64_t i2)
                {return dep_ids[i1] < dep_ids[i2];});
      std::sort(indep_order.begin(), indep_order.end(),
                [&indep_ids](const int64_t i1, const int64_t i2)
                {return indep_ids[i1] < indep_ids[i2];});

      for(auto i=0; i<nnodes; i++) {
        REQUIRE(dep_order[i]==indep_order[i]);
      }
    }
  }
}
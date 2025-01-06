/**
 * @file
 * @brief Contains Unit Tests for reading .msh files generated from gmsh
 */

#include "TPZGmshReader.h"
#include "TPZGeoMeshTools.h"
#include "pzcheckgeom.h"
#include "pzcmesh.h"
#include "TPZVTKGenerator.h"
#include <TPZVTKGeoMesh.h>
#include "Projection/TPZL2Projection.h"


#include <catch2/catch_test_macros.hpp>
#include <numeric> //for iota
///generates vtk files for visualising the solution on periodic meshes, also prints cmesh in .txt format
void PlotPeriodicSol(TPZGeoMesh *gmesh, const TPZVec<int> &vol_ids,
                     const SPZPeriodicData& data,
                     const std::string &name);
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
    //useful for debugging
    if constexpr (1){
      auto mesh_dim = meshReader.Dimension();
      auto matinfo_vec = meshReader.GetDimNamePhysical();
      TPZVec<int> vol_ids;
      for(auto [name,id] : matinfo_vec[mesh_dim]){
        vol_ids.push_back(id);
      }

      PlotPeriodicSol(gmesh.operator->(), vol_ids, data, "multiple_periodic");
    }
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

    //useful for debugging
    if constexpr (1){
      auto mesh_dim = meshReader.Dimension();
      auto matinfo_vec = meshReader.GetDimNamePhysical();
      TPZVec<int> vol_ids;
      for(auto [name,id] : matinfo_vec[mesh_dim]){
        vol_ids.push_back(id);
      }

      PlotPeriodicSol(gmesh.operator->(), vol_ids, data, "periodic_layer");
    }
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

void SetPeriodic(TPZCompMesh * cmesh, const std::map<int64_t,int64_t> &periodic_els)
{

  
  auto gmesh = cmesh->Reference();
  gmesh->ResetReference();
  cmesh->LoadReferences();
  const bool complex_mesh = cmesh->GetSolType() == ESolType::EComplex;

  for(auto [dep, indep] : periodic_els){
    //geometric elements
    auto *dep_gel = gmesh->Element(dep);
    const auto *indep_gel = gmesh->Element(indep);
    //computational element
    auto *indep_cel = indep_gel->Reference();
    auto *dep_cel = dep_gel->Reference();
    //not necessarily all elements belong to the comp mesh
    if(!indep_cel && !dep_cel){continue;}
    if(!indep_cel || !dep_cel){
      PZError<<__PRETTY_FUNCTION__
             <<"\nError: found only one of periodic els!\n";
      if(!indep_cel){
        PZError<<"Could not find comp el associated with geo el "<<indep<<std::endl;
      }
      if(!dep_cel){
        PZError<<"Could not find comp el associated with geo el "<<dep<<std::endl;
      }
      DebugStop();
    }
    //number of connects
    const auto n_dep_con = dep_cel->NConnects();
    const auto n_indep_con = indep_cel->NConnects();
    //just to be sure
    if(n_dep_con != n_indep_con){
      PZError<<__PRETTY_FUNCTION__
             <<"\nindep cel "<<indep_cel->Index()
             <<" has "<<n_indep_con<<" connects"
             <<"\ndep cel "<<dep_cel->Index()
             <<" has "<<n_dep_con<<" connects"<<std::endl;
      DebugStop();
    }

    //now we create dependencies between connects
    for(auto ic = 0; ic < n_indep_con; ic++){
      const auto indep_ci = indep_cel->ConnectIndex(ic);
      const auto dep_ci = dep_cel->ConnectIndex(ic);

      auto &dep_con = dep_cel->Connect(ic);
      //TODO: think of a more robust way on how to proceed in this scenario
      if(dep_con.HasDependency()){continue;}
      const auto ndof = dep_con.NDof(*cmesh);
      if(ndof==0) {continue;}
      constexpr int64_t ipos{0};
      constexpr int64_t jpos{0};


      if(complex_mesh){
        TPZFNMatrix<400,CSTATE> mat(ndof,ndof);
        mat.Identity();
        dep_con.AddDependency(dep_ci, indep_ci, mat, ipos,jpos,ndof,ndof);
      }else{
        TPZFNMatrix<400,STATE> mat(ndof,ndof);
        mat.Identity();
        dep_con.AddDependency(dep_ci, indep_ci, mat, ipos,jpos,ndof,ndof);
      }
    }//con 
  }//els
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
}

void PlotPeriodicSol(TPZGeoMesh *gmesh, const TPZVec<int> &vol_ids,
                     const SPZPeriodicData& data,
                     const std::string &name){

  TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> el_map_vec;
  TPZGeoMeshTools::FindPeriodicElements(gmesh, data.dep_mat_ids, data.indep_mat_ids, data.nodes_map,
                                        el_map_vec);
  //we will join all maps in one
  std::map<int64_t, int64_t> el_map;

  //we only care about elements of dimension dim-1
  const auto mesh_dim = gmesh->Dimension();

  
  for(const auto & el_map_orig : el_map_vec){
    for(auto [dep_idx, indep_idx] : *el_map_orig){
      if(gmesh->Element(dep_idx)->Dimension() != mesh_dim - 1){continue;}
      if(el_map.count(dep_idx) != 0){
          DebugStop();
      }
      el_map[dep_idx] = indep_idx;
      std::cout<<"element "<<dep_idx<<" depends on "<<indep_idx<<std::endl;
      
    }
  }
  TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
  cmesh->SetAllCreateFunctionsContinuous();
  TPZMaterial *mat{nullptr};
  for(auto id : vol_ids){
    mat = new TPZL2Projection<STATE> (id, mesh_dim);
    cmesh->InsertMaterialObject(mat);
  }
  //create boundary materials
  TPZFMatrix<STATE> v1;
  TPZVec<STATE> v2;
  for(auto id : data.dep_mat_ids){
    auto vmat = dynamic_cast<TPZMaterialT<STATE>*>(mat);
    auto *bcmat = vmat->CreateBC(vmat,id,0,v1,v2);
    cmesh->InsertMaterialObject(bcmat);
  }
  for(auto id : data.indep_mat_ids){
    auto vmat = dynamic_cast<TPZMaterialT<STATE>*>(mat);
    auto *bcmat = vmat->CreateBC(vmat,id,0,v1,v2);
    cmesh->InsertMaterialObject(bcmat);
  }
  //we need to check higher order functions for orientation
  cmesh->SetDefaultOrder(3);
  cmesh->AutoBuild();
  SetPeriodic(cmesh.operator->(), el_map);

  
  //we want to debug, right?
  {
    std::ofstream cmesh_file("cmesh_"+name+".txt");
    cmesh->Print(cmesh_file);
    std::ofstream gmesh_file_txt("gmesh_"+name+".txt");
    gmesh->Print(gmesh_file_txt);
    std::ofstream gmesh_file_vtk("gmesh_"+name+".vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, gmesh_file_vtk, true);
  }

  //so we are able to see the bubbles properly
  constexpr auto vtk_res{2};
  auto vtk =
    TPZVTKGenerator(cmesh, {"Solution"}, "sol_"+name, vtk_res);
  
  TPZFMatrix<STATE> sol = cmesh->Solution();

  for(auto [dep_idx, indep_idx] : el_map){
  //now we get the independent element
    auto gel = gmesh->Element(indep_idx);
    if(!gel){DebugStop();}
    auto cel  = gel->Reference();
    if(!cel){DebugStop();}
    const auto nc = cel->NConnects();
    int64_t pos{-1};
    auto &block = cmesh->Block();

    std::cout<<"setting solution for element "<<indep_idx<<std::endl;
    for(auto ic = 0; ic < nc; ic++){
      auto con_idx = cel->ConnectIndex(ic);
      auto &con = cmesh->ConnectVec()[con_idx];
      auto dep = con.FirstDepend();
      int64_t seqnum{-1};
      if(dep){
        auto indep_idx = dep->fDepConnectIndex;
        std::cout<<"\tconnect "<<con_idx
                 <<" depends on connect "<<indep_idx<<std::endl;
        auto &indep_con = cmesh->ConnectVec()[indep_idx];
        seqnum = indep_con.SequenceNumber();
      }else{
        std::cout<<"\tconnect "<<con_idx<<std::endl;
        seqnum = con.SequenceNumber();
      }
      const auto ndof = con.NDof();
      for(auto idof = 0; idof < ndof; idof++){
        pos = block.Position(seqnum) + idof;
        sol.Zero();
        sol.PutVal(pos,0,1);
        cmesh->LoadSolution(sol);
        vtk.Do();
      }
    }
  }

  //so it wont break
  for(auto &con : cmesh->ConnectVec()){
    con.RemoveDepend();
  }
  
}
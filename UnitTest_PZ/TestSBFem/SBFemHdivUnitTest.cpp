/**
 * @file MultiThreadUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of multi-thread computations
 *
 */
#include <iostream>
#include <chrono>

#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZBuildSBFemMultiphysics.h"
#include "TPZSpStructMatrix.h"
#include "TPZExactFunction.h"
#include "TPZDarcySBFemHdiv.h"
#include "TPZGeoMeshTools.h"
#include "TPZBndCond.h"
#ifdef USING_MKL
#include "TPZSSpStructMatrix.h"
#else
#include "pzskylstrmatrix.h"
#endif
#include "TPZAnalyticSolution.h"

#include <catch2/catch_test_macros.hpp>


namespace SBFemTest {

  void SBFemHdivDarcy2D(const int nThreads);
  
  TPZGeoMesh * CreateGMesh(const int nelx);

  TPZCompMesh * CreateCMeshPressure(TPZAutoPointer<TPZGeoMesh> & gmesh);
  TPZCompMesh * CreateCMeshFlux(TPZAutoPointer<TPZGeoMesh> & gmesh);
  TPZMultiphysicsCompMesh * CreateCMeshMultiphysics(TPZAutoPointer<TPZGeoMesh> & gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf);

  void Analysis(TPZLinearAnalysis & an, const int nThreads, TPZManVector<REAL> & errorVec);
  
  // Material Data
  constexpr int EGroup{1};
  constexpr int ESkeleton{6};
  constexpr int EMat1{0};
  constexpr int EBc1{2};

  // Solution data
  constexpr int pOrder{3};
  TLaplaceExample1 LaplaceExact;
}

TEST_CASE("SBFEM Hdiv convergence test","[sbfemhdiv][analysis]")
{
  std::cout << "#######################\nTesting SBFEM-Hdiv 2D approximation\n";
  SBFemTest::SBFemHdivDarcy2D(4);
}

TPZGeoMesh * SBFemTest::CreateGMesh(const int nelx)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh();

    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);

    // Setting boundary elements
    {
        gengrid.Read(gmesh, SBFemTest::EGroup);
        gengrid.SetBC(gmesh, 4, SBFemTest::EBc1);
        gengrid.SetBC(gmesh, 5, SBFemTest::EBc1);
        gengrid.SetBC(gmesh, 6, SBFemTest::EBc1);
        gengrid.SetBC(gmesh, 7, SBFemTest::EBc1);
    }
    gmesh->BuildConnectivity();

    return gmesh;
}

TPZCompMesh * SBFemTest::CreateCMeshPressure(TPZAutoPointer<TPZGeoMesh> & gmesh)
{
  auto *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);
  cmesh->SetDimModel(2);

  // Volumetric material
  auto mat = new TPZNullMaterial(SBFemTest::EMat1, 2, 1);
  cmesh->InsertMaterialObject(mat);

  // Boundary conditions
  TPZFMatrix<STATE> val1(2,2,0.);
  TPZManVector<STATE> val2(2,0.);
  {
      auto bcond = mat->CreateBC(mat, SBFemTest::EBc1, 0, val1, val2);
      bcond->SetForcingFunctionBC(LaplaceExact.ExactSolution(),7);
      cmesh->InsertMaterialObject(bcond);
  }
  {
      auto bcond = mat->CreateBC(mat, SBFemTest::ESkeleton, 1, val1, val2);
      cmesh->InsertMaterialObject(bcond);
  }

  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  set<int> matid = {SBFemTest::EMat1, SBFemTest::EBc1, SBFemTest::ESkeleton};
  cmesh->AutoBuild(matid);

  for(auto newnod : cmesh->ConnectVec())
  {
      newnod.SetLagrangeMultiplier(1);
  }

  return cmesh;
}

TPZCompMesh * SBFemTest::CreateCMeshFlux(TPZAutoPointer<TPZGeoMesh> & gmesh)
{
  auto dim = 2; auto nstate = 1;
  auto cmeshcollapsed = new TPZCompMesh(gmesh);
  cmeshcollapsed->SetDefaultOrder(SBFemTest::pOrder);
  cmeshcollapsed->SetDimModel(dim);
  cmeshcollapsed->CleanUp();
  
  // Volumetric material
  auto mat = new TPZNullMaterial(SBFemTest::EMat1, dim, nstate);
  cmeshcollapsed->InsertMaterialObject(mat);

  cmeshcollapsed->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
  cmeshcollapsed->AutoBuild();
  
  return cmeshcollapsed;
}

TPZMultiphysicsCompMesh * SBFemTest::CreateCMeshMultiphysics(TPZAutoPointer<TPZGeoMesh> & gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf)
{
  int dim = gmesh->Dimension();
  int nstate = 1;

  auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
  cmesh->SetDefaultOrder(SBFemTest::pOrder);
  cmesh->SetDimModel(dim);
  cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

  auto mat = new TPZHybridPoissonCollapsed(SBFemTest::EMat1,2);
  mat->TPZIsotropicPermeability::SetConstantPermeability(1.);
  cmesh->InsertMaterialObject(mat);

  mat->SetBigNumber(1e16);

  TPZFMatrix<STATE> val1(2,2,0.);
  TPZManVector<STATE> val2(2,0.);
  {
      auto bcond = mat->CreateBC(mat, SBFemTest::EBc1, 0, val1, val2);
      bcond->SetForcingFunctionBC(LaplaceExact.ExactSolution(),7);
      cmesh->InsertMaterialObject(bcond);
  }
  {
      auto bcond = mat->CreateBC(mat, ESkeleton, 1, val1, val2);
      cmesh->InsertMaterialObject(bcond);
  }

  cout << "Creating multiphysics mesh\n";
  TPZManVector< TPZCompMesh *, 2> meshvector(2);
  meshvector[0] = cmeshf;
  meshvector[1] = cmeshp;

  TPZManVector<MSpaceConfig> active(2,ETestTrial);
  cmesh->BuildMultiphysicsSpace(active, meshvector);
  cmesh->LoadReferences();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}

void SBFemTest::Analysis(TPZLinearAnalysis & an, const int nThreads, TPZManVector<REAL> &errorVec)
{
#ifdef USING_MKL
  TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matssp(an.Mesh());
  matssp.SetNumThreads(nThreads);
  an.SetStructuralMatrix(matssp);
#else
    TPZSkylineStructMatrix<> sklstr(an.Mesh());
    an.SetStructuralMatrix(sklstr);
#endif

  TPZStepSolver <STATE> *direct = new TPZStepSolver<STATE>;
  direct->SetDirect(ELDLt);
  an.SetSolver(*direct);
  delete direct;
  direct = 0;

  an.Assemble();
  an.Solve();
  an.SetThreadsForError(nThreads);

  an.PostProcessError(errorVec, false);
};

void SBFemTest::SBFemHdivDarcy2D(const int nThreads){
  LaplaceExact.fExact = TLaplaceExample1::EHarmonicPoly;

  constexpr int nelx = 2;

  constexpr int nDiv{1};
  constexpr int dim{2};

  TPZAutoPointer<TPZGeoMesh> gMesh = CreateGMesh(nelx);
  std::map<int,int> matmap;
  matmap[SBFemTest::EGroup] = SBFemTest::EMat1;
  TPZBuildSBFemMultiphysics build(gMesh, SBFemTest::ESkeleton, matmap);
  build.StandardConfiguration();
  build.DivideSkeleton(0);

  auto *cmeshp = CreateCMeshPressure(gMesh);
  cmeshp->LoadReferences();
  auto *cmeshf = CreateCMeshFlux(gMesh);
  auto *cmeshm = CreateCMeshMultiphysics(gMesh, cmeshp, cmeshf);
  

  build.BuildMultiphysicsCompMesh(*cmeshm);

  TPZLinearAnalysis an(cmeshm, RenumType::ENone);
  an.SetExact(LaplaceExact.ExactSolution());

  TPZManVector<REAL> errorVecPar;
  auto start = std::chrono::system_clock::now();
  Analysis(an, nThreads, errorVecPar);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedParallel = end - start;
  std::cout << "\tParallel time: " << elapsedParallel.count() << "s\n";
  std::cout << "\tParallel errors: " << errorVecPar << std::endl;

  bool pass = true;
  for (auto i = 0; i < 3; i++) {
    if (fabs(errorVecPar[i]) > 1e-10) {
      pass = false;
    }
  }
  REQUIRE(pass);

  delete cmeshm;
}

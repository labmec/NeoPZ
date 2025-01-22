/**
 * @file TestIterative.cpp
 * @brief Define a Unit Test using Catch2 for usage of preconditioners in NeoPZ
 *
 */

#include "pzcmesh.h"
#include "TPZGeoMeshTools.h"
#include "TPZLinearAnalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZIdentitySolver.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZSSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZSimpleTimer.h"
#include <iostream>


#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>


//useful for debugging
// #define TEST_1D

template<>
struct Catch::StringMaker<long double> {
    static std::string convert(long double ref);
    static int precision;
};

int Catch::StringMaker<long double>::precision = 10;

std::string Catch::StringMaker<long double>::convert(long double value) {
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << value;
    return out.str();
}

namespace testiterative{
  //! Creates computational mesh for 2d poisson problem
  TPZCompMesh *CreateCMesh();
  //! Creates computational mesh for 1d poisson problem
  TPZCompMesh *CreateCMesh1D();
  //! Creates vector of active equations by removing all dirichlet bc eqs
  void
  FilterBoundaryEquations(TPZCompMesh * cmesh,
                          TPZVec<int64_t> &activeEquations);
  constexpr int niter_wo_eqfilt{40};
  constexpr int niter{40};
  constexpr STATE tol{std::numeric_limits<STATE>::epsilon() * 10000};
}

TEST_CASE("Poisson equation","[iterative_testss]") {
  auto oldPrecision = Catch::StringMaker<STATE>::precision;
  Catch::StringMaker<STATE>::precision = std::numeric_limits<STATE>::max_digits10;
  

  TPZAutoPointer<TPZCompMesh> cmesh = testiterative::CreateCMesh();
  TPZLinearAnalysis an(cmesh,RenumType::ENone);
#if defined (PZ_USING_MKL) || (PZ_USING_EIGEN)
  TPZSSpStructMatrix<STATE> str(cmesh);
#else
  TPZSkylineStructMatrix<STATE>str(cmesh);
#endif

  an.SetStructuralMatrix(str);
  
  {
    TPZStepSolver<STATE> solv;
    solv.SetDirect(ECholesky);
    an.SetSolver(solv);
  }
  an.Assemble();

  auto *solver = dynamic_cast<TPZStepSolver<STATE>*>(an.Solver());

  solver->Matrix()->SetSymmetry(SymProp::Herm);
  solver->Matrix()->SetDefPositive(true);


  constexpr int64_t fromCurrent{0};
  SECTION("Identity matrix as precond"){
    std::cout<<"testing identity matrix as a precond"<<std::endl;
    bool success = false;
    TPZIdentitySolver<STATE> pre;
    solver->SetGMRES(testiterative::niter_wo_eqfilt,
                     testiterative::niter_wo_eqfilt, pre,
                     testiterative::tol, fromCurrent);
    an.Solve();
    success = true;
    REQUIRE(success);
  }
  SECTION("Exact inverse as precond"){
    std::cout<<"testing exact inverse as precond"<<std::endl;
    //we copy the original matrix
    TPZAutoPointer<TPZMatrix<STATE>> mtrxcp = solver->Matrix()->Clone();
    //we use a direct solver as a precond
    TPZStepSolver<STATE> pre(mtrxcp);
    pre.SetDirect(ELDLt);
    solver->SetGMRES(testiterative::niter,
                     testiterative::niter, pre,
                     testiterative::tol, fromCurrent);
    an.Solve();
    const auto niter_res = solver->NumIterations();
    REQUIRE(niter_res == 1);
    const auto tol_res = solver->GetTolerance();
    REQUIRE(tol_res <= testiterative::tol);
  }

  auto pre_type = GENERATE(Precond::Jacobi,
                           Precond::BlockJacobi,
                           Precond::Element,
                           Precond::NodeCentered);

  const bool overlap = GENERATE(true,false);
  const auto precondname = std::string(Precond::Name(pre_type)) +
                                       (overlap ? " with overlap" : " without overlap");
  SECTION(precondname){
    std::cout<<"testing "<<precondname<<std::endl;
    TPZAutoPointer<TPZMatrixSolver<STATE>> pre =
      an.BuildPreconditioner<STATE>(pre_type, overlap);
    const bool created_precond = pre;
    REQUIRE(created_precond);
    solver->SetGMRES(testiterative::niter_wo_eqfilt,
                     testiterative::niter_wo_eqfilt,
                     *pre, testiterative::tol, fromCurrent);
    an.Solve();
    const auto tol_res = solver->GetTolerance();
    REQUIRE(tol_res <= testiterative::tol);
  }
  Catch::StringMaker<STATE>::precision = oldPrecision;
}

TEST_CASE("Equation filter and precond","[iterative_testss]") {
  auto oldPrecision = Catch::StringMaker<STATE>::precision;
  Catch::StringMaker<STATE>::precision = std::numeric_limits<STATE>::max_digits10;
  TPZAutoPointer<TPZCompMesh> cmesh = testiterative::CreateCMesh();
  TPZLinearAnalysis an(cmesh,RenumType::ENone);
  {
    TPZStepSolver<STATE> solv;
    solv.SetDirect(ECholesky);
    an.SetSolver(solv);
  }
#if defined (PZ_USING_MKL) || (PZ_USING_EIGEN)
  TPZSSpStructMatrix<STATE> str(cmesh);
#else
  TPZSkylineStructMatrix<STATE>str(cmesh);
#endif
  TPZVec<int64_t> active_eqs;
  testiterative::FilterBoundaryEquations(cmesh.operator->(), active_eqs);
  str.EquationFilter().SetActiveEquations(active_eqs);
  an.SetStructuralMatrix(str);
  an.Assemble();

  auto *solver = dynamic_cast<TPZStepSolver<STATE>*>(an.Solver());

  solver->Matrix()->SetSymmetry(SymProp::Herm);
  solver->Matrix()->SetDefPositive(true);


  constexpr int64_t fromCurrent{0};
  SECTION("Exact inverse as precond"){
    std::cout<<"testing exact inverse as precond"<<std::endl;
    //we copy the original matrix
    TPZAutoPointer<TPZMatrix<STATE>> mtrxcp = solver->Matrix()->Clone();
    //we use a direct solver as a precond
    TPZStepSolver<STATE> pre(mtrxcp);
    pre.SetDirect(ELDLt);
    solver->SetGMRES(testiterative::niter, testiterative::niter, pre,
                     testiterative::tol, fromCurrent);
    an.Solve();
    const auto niter_res = solver->NumIterations();
    REQUIRE(niter_res == 1);
    const auto tol_res = solver->GetTolerance();
    REQUIRE(tol_res <= testiterative::tol);
  }

  auto pre_type = GENERATE(Precond::Jacobi,
                           Precond::BlockJacobi,
                           Precond::Element,
                           Precond::NodeCentered);
  
  const bool overlap = GENERATE(true,false);
  const auto precondname = std::string(Precond::Name(pre_type)) +
                                       (overlap ? " with overlap" : " without overlap");
  SECTION(precondname){
    std::cout<<"testing "<<precondname<<std::endl;
    TPZSimpleTimer timer(precondname,true);
    TPZAutoPointer<TPZMatrixSolver<STATE>> pre =
      an.BuildPreconditioner<STATE>(pre_type, overlap);
    const bool created_precond = pre;
    REQUIRE(created_precond);
    solver->SetGMRES(testiterative::niter, testiterative::niter,
                     *pre, testiterative::tol, fromCurrent);
    an.Solve();
    const auto tol_res = solver->GetTolerance();
    REQUIRE(tol_res <= testiterative::tol);
  }
  Catch::StringMaker<STATE>::precision = oldPrecision;
}


#ifdef TEST_1D
/*
  The following test case is mostly for debugging.
  No need to run it every time
 */
TEST_CASE("1D H1 problem","[iterative_tests]") {
  auto oldPrecision = Catch::StringMaker<STATE>::precision;
  Catch::StringMaker<STATE>::precision = std::numeric_limits<STATE>::max_digits10;
  TPZAutoPointer<TPZCompMesh> cmesh = testiterative::CreateCMesh1D();
  constexpr bool optimizeBandwidth{false};
  TPZLinearAnalysis an(cmesh,optimizeBandwidth);
  {
    TPZStepSolver<STATE> solv;
    solv.SetDirect(ECholesky);
    an.SetSolver(solv);
  }
#if defined (PZ_USING_MKL) || (PZ_USING_EIGEN)
  TPZSSpStructMatrix<STATE> str(cmesh);
#else
  TPZSkylineStructMatrix<STATE>str(cmesh);
#endif
  TPZVec<int64_t> active_eqs;
  testiterative::FilterBoundaryEquations(cmesh.operator->(), active_eqs);
  str.EquationFilter().SetActiveEquations(active_eqs);
  an.SetStructuralMatrix(str);
  an.Assemble();

  auto *solver = dynamic_cast<TPZStepSolver<STATE>*>(an.Solver());

  solver->Matrix()->SetSymmetry(SymProp::Herm);
  solver->Matrix()->SetDefPositive(true);


  constexpr int64_t fromCurrent{0};

  auto pre_type = GENERATE(Precond::Jacobi,
                           Precond::BlockJacobi,
                           Precond::Element,
                           Precond::NodeCentered);
  
  const bool overlap = GENERATE(true,false);
  const auto precondname = std::string(Precond::Name(pre_type)) +
                                       (overlap ? " with overlap" : " without overlap");
  SECTION(precondname){
    std::cout<<"1d: testing "<<precondname<<std::endl;
    TPZAutoPointer<TPZMatrixSolver<STATE>> pre =
      an.BuildPreconditioner<STATE>(pre_type, overlap);
    const bool created_precond = pre;
    REQUIRE(created_precond);
    solver->SetGMRES(testiterative::niter, testiterative::niter,
                     *pre, testiterative::tol, fromCurrent);
    an.Solve();
    const auto tol_res = solver->GetTolerance();
    REQUIRE(tol_res <= testiterative::tol);
  }
  Catch::StringMaker<STATE>::precision = oldPrecision;
}

#endif

namespace testiterative{

  TPZCompMesh *CreateCMesh()
  {
    constexpr int rhsPOrder{3};
    const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
      const REAL &x = loc[0];
      const REAL &y = loc[1];
      u[0] = 2*y*y+2*x*x-4;
      u[0] *= -1;
    };
    
    //dimension of the problem
    constexpr int dim{2};
    //n divisions in x direction
    constexpr int nDivX{3};
    //n divisions in y direction
    constexpr int nDivY{3};
  
    //TPZManVector<Type,N> is a vector container with static + dynamic storage. one can also use TPZVec<Type> for dynamic storage
    TPZManVector<int,2> nDivs={nDivX,nDivY};

    //all geometric coordinates in NeoPZ are in the 3D space

    //lower left corner of the domain
    TPZManVector<REAL,3> minX={-1,-1,0};
    //upper right corner of the domain
    TPZManVector<REAL,3> maxX={ 1, 1,0};

    /*vector for storing materials(regions) identifiers
     * for the TPZGeoMeshTools::CreateGeoMeshOnGrid function.
     * In this example, we have different materials in each
     * section of the boundary*/
    TPZManVector<int,5> matIdVec={1,-1,-2,-3,-4};
    //whether to create boundary elements
    constexpr bool genBoundEls{true};
    //type of elements
    constexpr MMeshType meshType{MMeshType::ETriangular};


    //polynomial order of elements
    constexpr int pOrder{3};
    //defining the geometry of the problem
    //TPZAutoPointer is a smart pointer from the NeoPZ library
    TPZAutoPointer<TPZGeoMesh> gmesh =
      TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,meshType,genBoundEls);

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    auto *mat = new TPZMatPoisson<STATE>(matIdVec[0],dim);
    //TPZMatLaplacian solves div(k grad(u)) = -f
    mat->SetForcingFunction(rhs,rhsPOrder);
    cmesh->InsertMaterialObject(mat);

    //now we insert the boundary conditions
    for(auto i = 0; i < 4; i++)
    {
      //TPZFMatrix<T> implements a full matrix of type T
      /* val1 and val2 are used for calculating the boundary
       * conditions. val1 goes in the matrix and val2 in the rhs.
       * for dirichlet boundary conditions, only the value of 
       * val2 is used.*/
      TPZFMatrix<STATE> val1(1,1,0.);
      TPZManVector<STATE,1> val2={0};
      //dirichlet=0,neumann=1,robin=2
      constexpr int boundType{0};
      //TPZBndCond is a material type for boundary conditions
      TPZBndCond * bnd = mat->CreateBC(mat,matIdVec[i+1],boundType,val1,val2);
      cmesh->InsertMaterialObject(bnd);
    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();
    return cmesh;
  }

  TPZCompMesh *CreateCMesh1D()
  {
    /*
      this rhs corresponds to the solution (x^2-1)^2
     */
    constexpr int rhsPOrder{2};
    const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
      const REAL &x = loc[0];
      const REAL &y = loc[1];
      u[0] = 4-12*x*x;
    };
    
    //dimension of the problem
    constexpr int dim{1};
    //number of elements
    constexpr int nEls{5};
    //all geometric coordinates in NeoPZ are in the 3D space

    //left boundary
    constexpr REAL minX{-1};
    //right boundary
    constexpr REAL maxX{ 1};

    /*vector for storing materials(regions) identifiers
     * for the TPZGeoMeshTools::CreateGeoMeshOnGrid function.
     * In this example, we have different materials in each
     * section of the boundary*/
    TPZManVector<int,5> matIdVec={1,-1,-1};
    //whether to create boundary elements
    constexpr bool genBoundEls{true};
    //type of elements
    constexpr MMeshType meshType{MMeshType::ETriangular};


    //polynomial order of elements
    constexpr int pOrder{3};
    //defining the geometry of the problem
    //TPZAutoPointer is a smart pointer from the NeoPZ library
    TPZAutoPointer<TPZGeoMesh> gmesh =
      TPZGeoMeshTools::CreateGeoMesh1D(minX, maxX, nEls, matIdVec, genBoundEls);

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    auto *mat = new TPZMatPoisson<STATE>(matIdVec[0],dim);
    //TPZMatLaplacian solves div(k grad(u)) = -f
    mat->SetForcingFunction(rhs,rhsPOrder);
    cmesh->InsertMaterialObject(mat);

    //now we insert the boundary conditions
    {
      //TPZFMatrix<T> implements a full matrix of type T
      /* val1 and val2 are used for calculating the boundary
       * conditions. val1 goes in the matrix and val2 in the rhs.
       * for dirichlet boundary conditions, only the value of 
       * val2 is used.*/
      TPZFMatrix<STATE> val1(1,1,0.);
      TPZManVector<STATE,1> val2={0};
      //dirichlet=0,neumann=1,robin=2
      constexpr int boundType{0};
      //TPZBndCond is a material type for boundary conditions
      TPZBndCond * bnd = mat->CreateBC(mat,matIdVec[1],boundType,val1,val2);
      cmesh->InsertMaterialObject(bnd);

    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();
    return cmesh;
  }
  
  void FilterBoundaryEquations(TPZCompMesh *cmesh,
                               TPZVec<int64_t> &activeEquations) {
    TPZManVector<int64_t, 1000> allConnects;
    std::set<int64_t> boundConnects;

    for (int iel = 0; iel < cmesh->NElements(); iel++) {
      TPZCompEl *cel = cmesh->ElementVec()[iel];
      if (cel == nullptr) {
        continue;
      }
      if (cel->Reference() == nullptr) {//there is no associated geometric el
        continue;
      }
      TPZBndCond *mat = dynamic_cast<TPZBndCond *>(
        cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
      if (mat && mat->Type() == 0) {//check for dirichlet bcs
        std::set<int64_t> boundConnectsEl;
        cel->BuildConnectList(boundConnectsEl);
        for(auto val : boundConnectsEl){
          if (boundConnects.find(val) == boundConnects.end()) {
            boundConnects.insert(val);
          }
        }
      }
    }

    //certainly we have less equations than this, but we will avoid repeated resizes
    activeEquations.Resize(cmesh->NEquations());
    int neq = 0;
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
      if (boundConnects.find(iCon) == boundConnects.end()) {
        TPZConnect &con = cmesh->ConnectVec()[iCon];
        const auto hasdep = con.HasDependency();
        const auto seqnum = con.SequenceNumber();
        const auto pos = cmesh->Block().Position(seqnum);
        const auto blocksize = cmesh->Block().Size(seqnum);
      
        if(hasdep || seqnum < 0 || !blocksize) { continue; }
        const auto vs = neq;
        for (auto ieq = 0; ieq < blocksize; ieq++) {
          activeEquations[vs + ieq] = pos + ieq;
        }
        neq += blocksize;
      }
    }
    activeEquations.Resize(neq);
  }
}

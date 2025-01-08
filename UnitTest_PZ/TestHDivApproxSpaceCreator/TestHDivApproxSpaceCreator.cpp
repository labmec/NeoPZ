// ----- PZ includes -----
#include <TPZGenGrid3D.h>
#include <TPZGenGrid2D.h>
#include <pzgmesh.h>
#include <TPZVTKGeoMesh.h>

#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <pzfstrmatrix.h>
#include "TPZEnumApproxFamily.h"

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <Elasticity/TPZMixedElasticityND.h>
#include "TPZHDivApproxCreator.h"
#include "TPZVTKGenerator.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzcheckgeom.h"
#include "pzcheckmesh.h"
#include "TPZAnalyticSolution.h"
#include "TPZMHMHDivApproxCreator.h"

#include <pzlog.h>
#include "TPZSimpleTimer.h"

// ----- Unit test includes -----
//#define USE_MAIN

#ifndef USE_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#endif

using namespace std;

/// Creates a simple mesh used for testing
TPZGeoMesh *Create2DGeoMesh(ProblemType& pType, MMeshType &mType);
TPZGeoMesh *Create3DGeoMesh(ProblemType& pType, MMeshType &mType);

void InsertMaterials(TPZHDivApproxCreator &approxCreator, ProblemType &ptype);

void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder,
                                bool isRigidBodySpaces, MMeshType mType, int extrapOrder,
                                bool isCondensed, HybridizationType hType, bool isRef, bool isMHM);

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, HDivFamily hdivfam);

void CheckError(TPZMultiphysicsCompMesh *cmesh, TPZVec<REAL> &error, ProblemType pType);

void Refinement(TPZGeoMesh *gmesh);

void SolveSystem(TPZMultiphysicsCompMesh* cmesh, const bool isTestKnownSol);

void TestKnownSol(TPZLinearAnalysis& an, const REAL cteSol, TPZMultiphysicsCompMesh* mpcmesh);

void PostProcessVTK(TPZMultiphysicsCompMesh* cmesh, ProblemType probType);

void CheckNEqCondensedProb(TPZMultiphysicsCompMesh* mpcmesh,
                           TPZHDivApproxCreator& hdivcreator,
                           MMeshType& elType);

enum MaterialIds {EDomain,EBCDirichlet,EBCNeumann,EBCDisplacementLeft,EBCDisplacementRight};

constexpr const char* HDivFamilyToChar(HDivFamily hdivfam) {
    switch (hdivfam){
        case HDivFamily::EHDivStandard: return "EHDivStandard";
        case HDivFamily::EHDivConstant: return "EHDivConstant";
        case HDivFamily::EHDivKernel: return "EHDivKernel";
        default: std::invalid_argument("Unimplemented item");
    }
    return "Unimplemented item";//silences compiler warning on gcc
}
constexpr const char* HybridizationTypeToChar(HybridizationType hType) {
    switch (hType){
        case HybridizationType::ENone: return "ENone";
        case HybridizationType::EStandard: return "EStandard";
        case HybridizationType::EStandardSquared: return "EStandardSquared";
        case HybridizationType::ESemi: return "ESemi";
        default: std::invalid_argument("Unimplemented item");
    }
    return "Unimplemented item";//silences compiler warning on gcc
}
constexpr const char* ProblemTypeToChar(ProblemType ptype) {
    switch (ptype){
        case ProblemType::EDarcy: return "EDarcy";
        case ProblemType::EElastic: return "EElastic";
        case ProblemType::EStokes: return "EStokes";
        case ProblemType::ENone: return "ENone";
        default: std::invalid_argument("Unimplemented item");
    }
    return "Unimplemented item";//silences compiler warning on gcc
}

auto exactSolDarcy = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = 1.;
    gradU.Zero();
};

auto exactSolElastic = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&disp,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    disp.Fill(0.);
    disp[0] = (1. + x)/2;
    disp[1] = 0.;
    
    gradU.Zero();
    gradU(0,0) = 0.5;
};


#ifndef USE_MAIN
TEST_CASE("Approx Space Creator", "[hdiv_space_creator_test][!mayfail]") {
    std::cout << "Testing HDiv Approx Space Creator \n";
    
    // HDivFamily sType = GENERATE(HDivFamily::EHDivConstant,HDivFamily::EHDivStandard);
    HDivFamily sType = GENERATE(HDivFamily::EHDivConstant,HDivFamily::EHDivStandard);
    ProblemType pType = GENERATE(ProblemType::EDarcy,ProblemType::EElastic);
    // ProblemType pType = GENERATE(ProblemType::EElastic);
    int pOrder = GENERATE(1);
    bool isRBSpaces = GENERATE(false);
    MMeshType mType = GENERATE(MMeshType::EQuadrilateral,MMeshType::ETriangular,MMeshType::EHexahedral,MMeshType::ETetrahedral);
    int extraporder = GENERATE(0);
//    bool isCondensed = GENERATE(true);
    bool isCondensed = GENERATE(false,true);
//    HybridizationType hType = GENERATE(HybridizationType::ENone);
//    HybridizationType hType = GENERATE(HybridizationType::ESemi);
    HybridizationType hType = GENERATE(HybridizationType::EStandard,HybridizationType::ESemi,HybridizationType::ENone);
    // HybridizationType hType = GENERATE(HybridizationType::EStandard,HybridizationType::ENone);
    bool isRef = GENERATE(true,false);
    // bool isRef = GENERATE(true);
    bool isMHM = GENERATE(false);

    if (!(isRef && hType == HybridizationType::ESemi) && !(isCondensed && sType == HDivFamily::EHDivConstant) &&
    !(sType == HDivFamily::EHDivConstant && pType == ProblemType::EElastic && mType == MMeshType::ETetrahedral) &&
    !(sType != HDivFamily::EHDivConstant && hType == HybridizationType::ESemi))
    {
        TestHdivApproxSpaceCreator(sType, pType, pOrder, isRBSpaces, mType, extraporder, isCondensed, hType, isRef, isMHM);
    }
    std::cout << "Finish test HDiv Approx Space Creator \n";
}
#else
int main(){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    //    HDivFamily sType = HDivFamily::EHDivStandard;
    //    HDivFamily sType = HDivFamily::EHDivKernel;
    HDivFamily sType = HDivFamily::EHDivConstant;
    
       ProblemType pType = ProblemType::EElastic;
    // ProblemType pType = ProblemType::EDarcy;
    
    const int pord = 1;
    const bool isRBSpaces = true;
    
      MMeshType mType = MMeshType::EQuadrilateral;
    // MMeshType mType = MMeshType::ETriangular;
    //    MMeshType mType = MMeshType::EHexahedral;
    // MMeshType mType = MMeshType::ETetrahedral;
    
    int extraporder = 0;
    bool isCondensed = true;
    // bool isCondensed = false;
//    HybridizationType hType = HybridizationType::EStandard;
    // HybridizationType hType = HybridizationType::ESemi;
    HybridizationType hType = HybridizationType::ENone;
    
    // this will create a mesh with hanging nodes
    bool isRef = true;
    // bool isRef = false;
    
    bool isMHM = false;
    
    TestHdivApproxSpaceCreator(sType,pType,pord,isRBSpaces,mType,extraporder,isCondensed,hType,isRef,isMHM);
    
    return 0;
}
#endif


TPZGeoMesh *Create2DGeoMesh(ProblemType& pType, MMeshType &mType) {
    
    // ----- Create Geo Mesh -----
    const TPZManVector<REAL,2> minX = {-1.,-1.};
    const TPZManVector<REAL,2> maxX = {1.,1.};
    const TPZManVector<int,2> nelDiv = {2,1};
    TPZGenGrid2D gen2d(nelDiv,minX,maxX);
    gen2d.SetElementType(mType);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen2d.Read(gmesh,EDomain);

    if(pType == ProblemType::EDarcy){
        gen2d.SetBC(gmesh, 4, EBCDirichlet);
        gen2d.SetBC(gmesh, 5, EBCDirichlet);
        gen2d.SetBC(gmesh, 6, EBCDirichlet);
        gen2d.SetBC(gmesh, 7, EBCDirichlet);
    }
    else if (pType == ProblemType::EElastic){
        const bool isrbmProb = false;
        if(isrbmProb){
            gen2d.SetBC(gmesh, 4, EBCDirichlet);
            gen2d.SetBC(gmesh, 5, EBCDirichlet);
            gen2d.SetBC(gmesh, 6, EBCDirichlet);
            gen2d.SetBC(gmesh, 7, EBCDirichlet);
        }
        else{
            gen2d.SetBC(gmesh, 4, EBCNeumann);
            gen2d.SetBC(gmesh, 5, EBCDisplacementRight);
            gen2d.SetBC(gmesh, 6, EBCNeumann);
            gen2d.SetBC(gmesh, 7, EBCDisplacementLeft);
        }
    }

    return gmesh;
}

TPZGeoMesh *Create3DGeoMesh(ProblemType& pType, MMeshType &mType) {
    
    // ----- Create Geo Mesh -----
    const TPZManVector<REAL,3> minX = {-1.,-1.,-1.};
    const TPZManVector<REAL,3> maxX = {1.,1.,1.};
    const TPZManVector<int,3> nelDiv = {2,1,1};
    TPZGenGrid3D gen3d(minX,maxX,nelDiv,mType);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gmesh = gen3d.BuildVolumetricElements(EDomain);

    if(pType == ProblemType::EDarcy){
        gmesh = gen3d.BuildBoundaryElements(EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet);
    }
    else if (pType == ProblemType::EElastic){
        const bool isrbmProb = false;
        if(isrbmProb){
            gmesh = gen3d.BuildBoundaryElements(EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet);
        }
        else {
            gmesh = gen3d.BuildBoundaryElements(EBCNeumann,EBCDisplacementLeft,EBCNeumann,EBCDisplacementRight,EBCNeumann,EBCNeumann);
        }
    }
    
    return gmesh;
}

auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    // REAL a1 = 1./4;
    // REAL alpha = M_PI/2;
    // u[0] = x*a1*cos(x*alpha)*cosh(y*alpha) + y*a1*sin(x*alpha)*sinh(y*alpha) + x*x - y*y;
    // gradU(0,0) = -a1*(cosh(alpha*y)*(cos(alpha*x) - alpha*x*sin(alpha*x)) + alpha*y*cos(alpha*x)*sinh(alpha*y));
    // gradU(1,0) = -a1*(alpha*y*cosh(alpha*y)*sin(alpha*x) + (alpha*x*cos(alpha*x) + sin(alpha*x))*sinh(alpha*y));

    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = (3.*x*x*y - y*y*y);
    gradU(1,0) = (x*x*x - 3.*y*y*x);

};

void InsertMaterials(TPZHDivApproxCreator &approxCreator, ProblemType &ptype){
    
    if (!approxCreator.GeoMesh()) {
        cout << "\nError! Please set the geomesh before inserting materials" << endl;
        DebugStop();
    }
    const int dim = approxCreator.GeoMesh()->Dimension();

    approxCreator.ProbType() = ptype;

    TPZMaterial *mat = nullptr;
    TPZMixedDarcyFlow* matdarcy = nullptr;
    TPZMixedElasticityND* matelas = nullptr;
    if (ptype == ProblemType::EDarcy) {
        matdarcy = new TPZMixedDarcyFlow(EDomain,dim);
        matdarcy->SetConstantPermeability(1.);
        matdarcy->SetExactSol(exactSolDarcy,1);
        mat = matdarcy;
    }
    else if (ptype == ProblemType::EElastic){
        REAL E = 1., nu = 0., fx = 0., fy = 0.;
        const int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress
        matelas = new TPZMixedElasticityND(EDomain, E, nu, fx, fy, plain, dim);
        matelas->SetExactSol(exactSolElastic,1);
        mat = matelas;
    }

    approxCreator.InsertMaterialObject(mat);

    // ========> Boundary Conditions
    // -----------------------------
    
    TPZBndCondT<STATE> *BCond1 = nullptr, *BCond2 = nullptr, *BCond3 = nullptr, *BCond4 = nullptr;
    const int dirType = 0, neuType = 1;

    if (ptype == ProblemType::EDarcy) {
        
        // matdarcy->SetExactSol(exactSol,4);

        TPZFMatrix<STATE> val1(1,1,0.);
        TPZManVector<STATE> val2(1,1.);
        BCond1 = matdarcy->CreateBC(matdarcy, EBCDirichlet, dirType, val1, val2);
        // BCond1->SetForcingFunctionBC(exactSol,4);
        val2[0] = 0.;
        BCond2 = matdarcy->CreateBC(matdarcy, EBCNeumann, neuType, val1, val2);
        // BCond2->SetForcingFunctionBC(exactSol,4);
    }
    else if (ptype == ProblemType::EElastic){
        TPZFMatrix<STATE> val1(dim,dim,0.);
        TPZManVector<STATE> val2(dim,0.);

        val2[0] = 0.;
        BCond1 = matelas->CreateBC(matelas, EBCNeumann, neuType, val1, val2);
        val2[0] = 1.;
        BCond2 = matelas->CreateBC(matelas, EBCDisplacementRight, dirType, val1, val2);
        val2[0] = 0.;
        BCond3 = matelas->CreateBC(matelas, EBCDisplacementLeft, dirType, val1, val2);
        val2[0] = 1.;
        val2[1] = 1.;
        if(dim == 3) val2[2] = 1.;
        BCond4 = matelas->CreateBC(matelas, EBCDirichlet, dirType, val1, val2);
    }
    else{
        DebugStop(); // yet not supported material
    }
    
    if(BCond1) approxCreator.InsertMaterialObject(BCond1);
    if(BCond2) approxCreator.InsertMaterialObject(BCond2);
    if(BCond3) approxCreator.InsertMaterialObject(BCond3);
    if(BCond4) approxCreator.InsertMaterialObject(BCond4);
}


void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder,
                                bool isRigidBodySpaces, MMeshType mType, int extrapOrder,
                                bool isCondensed, HybridizationType hType, bool isRef, bool isMHM){

    
    TPZSimpleTimer totaltime;
    // ==========> Initial headers <==========
    // =======================================
    static int globcount = 0;
    cout << "\n------------------ Starting test " << globcount++ << " ------------------" << endl;
    cout << "HdivFam = " << HDivFamilyToChar(hdivFam) <<
    "\nProblemType = " << ProblemTypeToChar(probType) <<
    "\npOrder = " << pOrder << "\nisRBSpaces = " << std::boolalpha << isRigidBodySpaces << 
    "\nMeshType = " << mType << "\nExtra POrder = " << extrapOrder <<
    "\nisCondensed = " << std::boolalpha << isCondensed <<
    "\nHybridization type = " << HybridizationTypeToChar(hType) <<
    "\nisRef = " << std::boolalpha << isRef << endl << endl;
    
    // TODO: WARNING!!!! Things to be fixed and for now we are skipping
     if(isCondensed && hdivFam == HDivFamily::EHDivConstant){
         std::cout << "\n\t======> WARNING! Condensing hdivconst spaces leads to singular K00. This happens because of the rotation space that has linear and constant functions."
         " One option is to separate the linear and constant functions is two different spaces. Please implement before using." << std::endl;
         cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
         DebugStop();
         return;
     }
    if(isRef && hType == HybridizationType::ESemi){
        cout << "\n\t======> This test is not working because the global matrix has wrong size!!\n" << endl;
        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
        DebugStop();
        return;
    }
    if(hdivFam == HDivFamily::EHDivConstant && probType == ProblemType::EElastic && mType == MMeshType::ETetrahedral){
        cout << "\n\t======> We don't know why this configuration does not work!!\n" << endl;
        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
        DebugStop();
        return;
    }
    if(hdivFam == HDivFamily::EHDivConstant && extrapOrder > 0){
        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
        DebugStop();
        return;
    }
    if (hdivFam == HDivFamily::EHDivKernel && isRigidBodySpaces) {
        std::cout << "ERROR! Hdiv kernel currently does not support rigid body spaces \n";
        DebugStop();
        return;
    }
    if (hdivFam != HDivFamily::EHDivConstant && hType == HybridizationType::ESemi) {
        std::cout << "ERROR! The only HDiv space with available Semi hybridization is HDivConstant \n";
        DebugStop();
        return;
    }
    
    // ==========> Creating GeoMesh <==========
    // ========================================
    TPZGeoMesh *gmesh;
    if (mType == MMeshType::EQuadrilateral || mType == MMeshType::ETriangular){
        gmesh = Create2DGeoMesh(probType,mType);
    } else if (mType == MMeshType::EHexahedral || mType == MMeshType::ETetrahedral){
        gmesh = Create3DGeoMesh(probType,mType);
    } else {
        DebugStop();
    }
    
    if(isRef) Refinement(gmesh);
    
    std::ofstream out("GeoMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

    // ==========> Creating Multiphysics mesh <==========
    // ==================================================
    TPZHDivApproxCreator* hdivCreator = nullptr;
    if(isMHM){
        hdivCreator = new TPZMHMHDivApproxCreator(gmesh);
    }
    else{
        hdivCreator = new TPZHDivApproxCreator(gmesh);
    }
    hdivCreator->HdivFamily() = hdivFam;
    hdivCreator->ProbType() = probType;
    hdivCreator->IsRigidBodySpaces() = isRigidBodySpaces;
    hdivCreator->SetDefaultOrder(pOrder);
    hdivCreator->SetExtraInternalOrder(extrapOrder);
    hdivCreator->SetShouldCondense(isCondensed);
    hdivCreator->HybridType() = hType;
//    hdivCreator->SetHybridizeBoundary();
    InsertMaterials(*hdivCreator,probType);
    TPZMultiphysicsCompMesh *cmesh = hdivCreator->CreateApproximationSpace();
    std::ofstream outtxt("geomeshmodified.txt");
//    gmesh->Print(outtxt);

    // ==========> Check number of equations for condensed problems <==========
    // ========================================================================
    if(isCondensed && hType != HybridizationType::ENone) {
        CheckNEqCondensedProb(cmesh,*hdivCreator,mType);
    }
    
#ifdef PZDEBUG
//    hdivCreator->PrintMeshElementsConnectInfo(cmesh);
    hdivCreator->PrintAllMeshes(cmesh);
#endif
    
    // ==========> Solving problem <==========
    // =======================================
    const bool isTestKnownSol = false;
    SolveSystem(cmesh,isTestKnownSol);
    if(isTestKnownSol) return; // do not perform checks if just testing known sol
    
    // ==========> Post processing <==========
    // =======================================
#ifdef USE_MAIN
    PostProcessVTK(cmesh,probType);
#endif

    // ==========> Unit test checks <==========
    // ========================================
    // Checks if the integral over the domain is a known value (most of the cases a constant value so just the volume of the domain)
    CheckIntegralOverDomain(cmesh,probType,hdivFam);
    // Checks if error with respect to exact solution is close to 0
    TPZManVector<REAL,5> error;
    CheckError(cmesh,error,probType);
    
    cout << "\n------------------ Test ended without crashing ------------------" << endl;
    std::cout << "==> Total time: " << totaltime.ReturnTimeDouble()/1000. << " seconds" << std::endl << endl;
}

void SolveSystem(TPZMultiphysicsCompMesh* cmesh, const bool isTestKnownSol) {
#ifdef USE_MAIN
    constexpr int nThreads{0};
#else
    constexpr int nThreads{12};
#endif
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matsp(cmesh);
    matsp.SetNumThreads(nThreads);

    std::cout << "\n=====> Number of equations = " << cmesh->NEquations() << std::endl << std::endl;
#ifdef USE_MAIN
    TPZLinearAnalysis an(cmesh,RenumType::ENone);
#else
    TPZLinearAnalysis an(cmesh);
#endif

    an.SetStructuralMatrix(matsp);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
        
    if(isTestKnownSol){
        TestKnownSol(an,1.,cmesh);
    }
    else{
        an.Run();
    }
}

void PostProcessVTK(TPZMultiphysicsCompMesh* cmesh, ProblemType probType) {
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
    const std::string plotfile = "PostProcess"; //sem o .vtk no final
    constexpr int vtkRes{0};

    TPZManVector<std::string,2> fields = {"Flux","Pressure"};
    if(probType == ProblemType::EElastic){
        fields[0] = "SigmaX";
        fields[1] = "Displacement";
    }
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

    vtk.Do();

}

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, HDivFamily hdivfam){

    int dim = cmesh -> Dimension();
    TPZVec<std::string> fields(2);
    switch (probType){
    case ProblemType::EDarcy:
        fields[0] = "Flux";
        fields[1] = "Pressure";
        break;
    case ProblemType::EElastic:
        fields[0] = "Stress";
        fields[1] = "Displacement";
        break;

    default:
        DebugStop();
        break;
    }
   
    std::set<int> matids = {EDomain};
    cmesh->Reference()->ResetReference();
    TPZVec<STATE> vecint = cmesh->Integrate(fields[0], matids);

    std::cout << "\n--------------- Integral of Flux --------------" <<  std::endl;
    std::cout << "Number of components = " << vecint.size() <<  std::endl;
    for (int i = 0; i < vecint.size(); i++)
    {
        std::cout << "Integral(" << i << ") = "  << vecint[i] << std::endl;
#ifndef USE_MAIN
        if (probType == ProblemType::EDarcy){
            REQUIRE(fabs(vecint[i]) < 1.e-10);
        }        
#endif
    }
    std::cout << std::endl;
#ifndef USE_MAIN
    if (probType == ProblemType::EElastic){
        if (dim == 2) REQUIRE(fabs(vecint[0]) == Catch::Approx( 2.0 ));
        if (dim == 3) REQUIRE(fabs(vecint[0]) == Catch::Approx( 4.0 ));
        for (int i = 1; i < vecint.size(); i++) REQUIRE(fabs(vecint[i]) < 1.e-10);
    }   
#endif
    if (hdivfam != HDivFamily::EHDivKernel){
        TPZVec<STATE> vecintp = cmesh->Integrate(fields[1], matids);
        std::cout << "\n--------------- Integral of State Var --------------" <<  std::endl;
        std::cout << "Number of components = " << vecintp.size() <<  std::endl;
        for (int i = 0; i < vecintp.size(); i++)
        {
            std::cout << "Integral(" << i << ") = "  << vecintp[i] << std::endl;
#ifndef USE_MAIN
            if (probType == ProblemType::EDarcy){
                if (dim == 2) REQUIRE(fabs(vecintp[i]) == Catch::Approx( 4.0 ));
                if (dim == 3) REQUIRE(fabs(vecintp[i]) == Catch::Approx( 8.0 ));
            }
#endif
        }
        std::cout << std::endl;
#ifndef USE_MAIN
        if (probType == ProblemType::EElastic){
            if (dim == 2) REQUIRE(fabs(vecintp[0]) == Catch::Approx( 2.0 ));
            if (dim == 3) REQUIRE(fabs(vecintp[0]) == Catch::Approx( 4.0 ));
            REQUIRE(fabs(vecintp[1]) < 1.e-10);
            REQUIRE(fabs(vecintp[2]) < 1.e-10);
        }
#endif
    }
}


void Refinement(TPZGeoMesh *gmesh){
    // children[0]->Divide(children);
    const int nel = gmesh->NElements();
    // for (int i = 0; i < nel; i++){
        int dim = gmesh->Dimension();
        TPZManVector<TPZGeoEl*,10> children;
        if (gmesh->ElementVec()[0]->Dimension() == dim) gmesh->ElementVec()[0]->Divide(children);
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->Dimension() != dim-1) continue;
            if(gel->HasSubElement()) continue;
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if(neighbour.HasSubElement()) {
                TPZManVector<TPZGeoEl*,10> children;
                gel->Divide(children);
            }
        }
    // }
    
    
    
//    children[0]->Divide(children);
}

void CheckError(TPZMultiphysicsCompMesh *cmesh, TPZVec<REAL> &error, ProblemType pType){

    cmesh->LoadReferences();
    cmesh->EvaluateError(false,error);
    const int dim = cmesh->Dimension();

    const bool isHDivConst = cmesh->MeshVector()[0]->ApproxSpace().HDivFam() == HDivFamily::EHDivConstant;
    std::cout << "Error = \n";
    for (int i = 0 ; i < error.size() ; i++) {
        std::cout << error[i] <<"\n";
#ifndef USE_MAIN
        if( pType == ProblemType::EElastic && i == 6 ) {
            // In Elastic mat i == 6 means energy error of the exact solution
            if(dim == 2){
                REQUIRE(error[i] == Catch::Approx(1.));
            }
            else if(dim == 3) {
                REQUIRE(error[i] == Catch::Approx(M_SQRT2));
            }
            continue;
        }
        
        // if elastic we dont check error in displacement for HDivConstant spaces
        if(pType == ProblemType::EElastic && i == 3 && isHDivConst){
            continue;
        }
        REQUIRE(fabs(error[i]) < 1.e-10);
#endif
    }
    std::cout << std::endl;
}

void TestKnownSol(TPZLinearAnalysis& an, const REAL cteSol, TPZMultiphysicsCompMesh* mpcmesh) {
    
    // Assemble matrix and rhs
    an.Assemble();
    
    // Fill the connects related with displacement with solution cteSol
    int64_t nc = mpcmesh->NConnects();
    TPZFMatrix<STATE> &sol = mpcmesh->Solution();
    for (int64_t ic = 0; ic<nc; ic++) {
        TPZConnect &c = mpcmesh->ConnectVec()[ic];
        int64_t seqnum = c.SequenceNumber();
        if(seqnum < 0) continue;
        unsigned char lagrange = c.LagrangeMultiplier();
        STATE fill = 0.;
        if(lagrange == 1)
        {
            fill = cteSol;
        }
        int ndof = c.NShape() * c.NState();
        for (int idf = 0; idf < ndof ; idf++) {
            int64_t index = mpcmesh->Block().Index(seqnum, idf);
            sol(index,0) = fill;
        }
    }
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mpcmesh->MeshVector(), mpcmesh);
    
    // Get Matrix
    TPZMatrix<STATE>* mat = an.MatrixSolver<STATE>().Matrix().operator->();
    
    // Multiply matrix by known solution vector
    const int neq = mpcmesh->NEquations();
    TPZFMatrix<STATE> res(neq,1,0.);
    mat->Multiply(mpcmesh->Solution(), res);
    res = res - an.Rhs();
    std::ofstream out("problematicElsGlob.txt");
    an.PrintVectorByElement(out, res, 1.e-6);
        
    for(int i = 0 ; i < mpcmesh->MeshVector().size() ; i++){
        std::string name = "Mesh_" + to_string(i) + ".txt";
        std::ofstream outmesh(name);
        mpcmesh->MeshVector()[i]->Print(outmesh);
    }
}

void CheckNEqCondensedProb(TPZMultiphysicsCompMesh* mpcmesh,
                           TPZHDivApproxCreator& hdivcreator,
                           MMeshType& elType) {
    
    std::cout << "\n--------------- Checking number of equations --------------" <<  std::endl;
    TPZCompMesh* pmesh = mpcmesh->MeshVector()[1];
    const int lagmatid = hdivcreator.HybridData().fLagrangeMatId;
    const int dim = mpcmesh->Dimension();
    
    int nlag = 0;
    for(auto cel : pmesh->ElementVec()) {
        if(!cel) continue;
        TPZGeoEl* gel = cel->Reference();
        if(!gel) DebugStop();
        const int gelmatid = gel->MaterialId();
        if (gelmatid == lagmatid) {
            ++nlag;
        }
    }
    
    const int nEquations = mpcmesh->NEquations();
    const int pOrder = mpcmesh->GetDefaultOrder();
    
    int nstate = 0;
    if (hdivcreator.ProbType()==ProblemType::EDarcy){
        nstate = 1;
    } else if (hdivcreator.ProbType()==ProblemType::EElastic) {
        nstate = dim;
    } else {
        DebugStop();
    }
    
    int expNEquations = nlag * (pOrder+1) * nstate;
    if(elType == MMeshType::ETetrahedral){
        expNEquations = nstate * nlag * (pOrder+1) * (pOrder+2) / 2;
    }
    if(elType == MMeshType::EHexahedral) {
        expNEquations = nstate * nlag * (pOrder+1) * (pOrder+1);
    }
    
    std::cout << "Expected equations: " << expNEquations << std::endl;
    std::cout << "Mesh equations: " << nEquations << std::endl;

#ifndef USE_MAIN
    REQUIRE(nEquations == expNEquations);
#endif
}

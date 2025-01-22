//
//  TestHangingNode.cpp
//  PZ
//
//  Created by Jeferson Fernandes and Nathan Shauer on 25/08/22.
//  Copyright 2022 UNICAMP. All rights reserved.
//

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"
#include "pzvec_extras.h"

using namespace pzshape;

#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h>
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include "Projection/TPZL2Projection.h"
#include "Projection/TPZL2ProjectionHDiv.h"
#include "Projection/TPZL2ProjectionHCurl.h"
#include "Projection/TPZL2ProjectionCS.h"
#include <TPZGeoMeshTools.h>
#include "TPZEnumApproxFamily.h"
#include "TPZLinearAnalysis.h"
#include "TPZVTKGenerator.h"
#include "TPZNullMaterial.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzlog.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckrestraint.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPattern.h"
#include "pzreftetrahedra.h"
#include "pzreftriangle.h"
#include "pzrefquad.h"
#include "TPZRefCube.h"

#undef USE_MAIN

//#define EXPORTVTK
//#define DEBUGTEST

#ifndef USE_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#endif

using namespace pzrefine;

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] meshType element type to be created.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template<class tshape>
TPZGeoMesh* CreateGeoMesh(const TPZVec<int> &nDivs, const int volId, const int bcId);


enum SpaceType {EH1, EHDivConstant, EHDivKernel, EHDivStandard, EHCurl, EHCurlNoGrads};

/*
    Test KernelHdiv problem
*/
template<class tshape>
void TestConstrainedSpace(const int &xdiv, const int &pOrder, SpaceType stype);
MMeshType GetMeshType(const MElementType elType);
void CheckSideOrientation(TPZGeoElSide &side1, TPZGeoElSide &side2);
void CheckElementInterfaces(TPZCompMesh *cmesh);

void Refinement(TPZGeoMesh *gmesh, SpaceType stype);

TPZCompMesh * CreateCMeshH1(TPZGeoMesh* gmesh, int pOrder, const int volId);
TPZCompMesh * CreateCMeshHDiv(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam);
TPZCompMesh * CreateCMeshHDiv2(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam);
TPZCompMesh * CreateCMeshHCurl(TPZGeoMesh* gmesh, int pOrder, const int volId, HCurlFamily hcurlfam);

constexpr const char* STypeToChar(SpaceType stype) {
  switch (stype){
    case SpaceType::EHDivStandard: return "EHDivStandard";
    case SpaceType::EHDivConstant: return "EHDivConstant";
    case SpaceType::EHDivKernel: return "EHDivKernel";
    case SpaceType::EH1: return "EH1";
    case SpaceType::EHCurl: return "EHCurl";
    case SpaceType::EHCurlNoGrads: return "EHCurlNoGrads";
    default: std::invalid_argument("Unimplemented item");
  }
  DebugStop();
  return "Unimplemented item";
}

// Test Hanging Nodes: FOR DEBUGGING PURPOSES
#ifndef USE_MAIN
TEST_CASE("Constrained Space", "[constrained_space_test]") {
    std::cout << "Testing Hanging Nodes \n";
    
    const int xdiv = GENERATE(2);
    const int pOrder = GENERATE(1,2);
    SpaceType sType = GENERATE(EHDivConstant,EHDivStandard,
                               EHCurl,
                               EH1);
    
    TestConstrainedSpace<pzshape::TPZShapeTriang>(xdiv,pOrder,sType);
    TestConstrainedSpace<pzshape::TPZShapeQuad>(xdiv,pOrder,sType);
    TestConstrainedSpace<pzshape::TPZShapeTetra>(xdiv,pOrder,sType);
    //NEEDS FIX
    if(sType != EHCurl){TestConstrainedSpace<pzshape::TPZShapeCube>(xdiv, pOrder, sType);}
    std::cout << "Finish test Constrained Space \n";
}
#else
int main(){

    const int xdiv = 2;
    const int pOrder = 1;
    SpaceType sType = EHCurl;
    // SpaceType sType = EHDivStandard;
    // SpaceType sType = EHCurlNoGrads;
    // SpaceType sType = EHDivConstant;

    TestConstrainedSpace<pzshape::TPZShapeCube>(xdiv,pOrder,sType);

    return 0;
}
#endif

MMeshType GetMeshType(const MElementType elType){
    MMeshType meshType;

    switch (elType)
    {
    case ETriangle:
        meshType = MMeshType::ETriangular;
        break;
    case EQuadrilateral:
        meshType = MMeshType::EQuadrilateral;
        break;
    case ETetraedro:
        meshType = MMeshType::ETetrahedral;
        break;
    case ECube:
        meshType = MMeshType::EHexahedral;
        break;
        case EPrisma:
        meshType = MMeshType::EPrismatic;
        break;
    default:
        DebugStop();
    }

    return meshType;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

template<class tshape>
TPZGeoMesh* CreateGeoMesh(const TPZVec<int> &nDivs, const int volId, const int bcId)
{
    MMeshType meshType = GetMeshType(tshape::Type());

    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*tshape::Dimension+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(tshape::Dimension, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);

    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = 1.;
};

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

template<class tshape>
void TestConstrainedSpace(const int &xdiv, const int &pOrder, SpaceType stype){

    MMeshType type = GetMeshType(tshape::Type());
    static int globcount = 0;
    std::cout << "\n------------------ Starting test " << globcount++ << " ------------------" << std::endl;
    std::cout << "SpaceType = " << STypeToChar(stype) <<
    "\nxdiv = " << xdiv <<
    "\nTopology = " << type <<
    "\npOrder = " << pOrder << std::endl << std::endl;

    constexpr int volId{1};
    constexpr int bcId{2};
    int DIM = tshape::Dimension;

    // TPZRefPatternDataBase ref;
    // ref.InitializeUniformRefPattern(tshape::Type());

    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,1};
    if (DIM == 3) nDivs = {xdiv,1,1};
    
    auto gmesh = CreateGeoMesh<tshape>(nDivs, volId, bcId);

    Refinement(gmesh, stype);

#ifdef EXPORTVTK
    //Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
#endif

    TPZCompMesh *cmesh;

    switch (stype)
    {
    case EH1:
        cmesh = CreateCMeshH1(gmesh, pOrder, volId);
        break;

    case EHDivStandard:
        cmesh = CreateCMeshHDiv2(gmesh,pOrder,volId,HDivFamily::EHDivStandard);
        break;

    case EHDivConstant:
        cmesh = CreateCMeshHDiv2(gmesh,pOrder,volId,HDivFamily::EHDivConstant);
        break;

    case EHDivKernel:
        cmesh = CreateCMeshHDiv2(gmesh,pOrder,volId,HDivFamily::EHDivKernel);
        break;
            
    case EHCurl:
        cmesh = CreateCMeshHCurl(gmesh,pOrder,volId,HCurlFamily::EHCurlStandard);
        break;

    case EHCurlNoGrads:
        cmesh = CreateCMeshHCurl(gmesh,pOrder,volId,HCurlFamily::EHCurlNoGrads);
        break;

    default:
        break;
    };
    
    CheckElementInterfaces(cmesh);

    TPZLinearAnalysis an(cmesh);

    an.Assemble();
    an.Solve();

    TPZVec<std::string> fields = {"Solution"};

    std::set<int> matids;
    matids.insert(volId);
    TPZVec<STATE> vecint = cmesh->Integrate(fields[0], matids);

    std::cout << "\n--------------- Integral of Solution --------------" <<  std::endl;
    std::cout << "Number of components = " << vecint.size() <<  std::endl;
    for (int i = 0; i < vecint.size(); i++)
    {
        std::cout << "Integral(" << i << ") = "  << vecint[i] << std::endl;
    }
    std::cout << std::endl;

#ifdef EXPORTVTK
    const std::string plotfile = "solution";//sem o .vtk no final
    constexpr int vtkRes{0};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.Do();
#endif

#ifndef USE_MAIN
    REAL tol = 1.e-10;
    for (int i = 0; i < vecint.size(); i++) REQUIRE(fabs(vecint[i]-1.)<tol);
#endif
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void Refinement(TPZGeoMesh *gmesh, SpaceType stype){

    // children[0]->Divide(children);

    if (stype != EHCurlNoGrads) {
        TPZManVector<TPZGeoEl*,10> children;
        gmesh->ElementVec()[0]->Divide(children);
        
        children[0]->Divide(children); 
    }

}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * CreateCMeshH1(TPZGeoMesh* gmesh, int pOrder, const int volId){

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(pOrder);

    TPZL2Projection<> *mat = new TPZL2Projection<>(volId,gmesh->Dimension());
    mat->SetForcingFunction(forcefunction,4);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * CreateCMeshHDiv(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam){
    //Create flux mesh
    TPZCompMesh *cmeshflux = new TPZCompMesh(gmesh);
    cmeshflux->SetDimModel(gmesh->Dimension());
    cmeshflux->SetDefaultOrder(pOrder);
    
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(volId);
    cmeshflux->InsertMaterialObject(mat);

    cmeshflux->ApproxSpace().SetHDivFamily(hdivfam);
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(gmesh->Dimension());
    cmeshflux->AutoBuild();

    //Create pressure cmesh
    TPZCompMesh *cmeshpres = new TPZCompMesh(gmesh);
    cmeshpres->SetDimModel(gmesh->Dimension());

    TPZNullMaterial<> *mat2 = new TPZNullMaterial<>(volId);
    cmeshpres->InsertMaterialObject(mat2);

    switch (hdivfam)
    {
    case HDivFamily::EHDivStandard:
        cmeshpres->SetDefaultOrder(pOrder);
        cmeshpres->SetAllCreateFunctionsContinuous();
        cmeshpres->ApproxSpace().CreateDisconnectedElements(true);
        break;
    case HDivFamily::EHDivConstant:
        cmeshpres->SetDefaultOrder(0);
        cmeshpres->SetAllCreateFunctionsDiscontinuous();
        break;
    
    default:
        DebugStop();//HDivKernel nao tem malha de pressao
    }
    
    cmeshpres->AutoBuild();
    int ncon = cmeshpres->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmeshpres->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(1);
    }

    TPZManVector<TPZCompMesh *,2> meshvector = {cmeshflux,cmeshpres};

    //Create Multiphysics cmesh
    TPZMultiphysicsCompMesh* cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(gmesh->Dimension());

    auto *matL2 = new TPZL2ProjectionCS<>(volId,gmesh->Dimension());
    matL2->SetForcingFunction(forcefunction,4);
    cmesh->InsertMaterialObject(matL2);

    TPZManVector<MSpaceConfig> active(meshvector.size(),ETestTrial);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * CreateCMeshHDiv2(TPZGeoMesh* gmesh, int pOrder, const int volId, HDivFamily hdivfam){
    //Create flux mesh
    TPZCompMesh *cmeshflux = new TPZCompMesh(gmesh);
    cmeshflux->SetDimModel(gmesh->Dimension());
    cmeshflux->SetDefaultOrder(pOrder);
    
    TPZL2ProjectionHDiv<> *mat = new TPZL2ProjectionHDiv<>(volId,gmesh->Dimension());
    mat->SetForcingFunction(forcefunction,4);
    cmeshflux->InsertMaterialObject(mat);

    cmeshflux->ApproxSpace().SetHDivFamily(hdivfam);
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHDiv(gmesh->Dimension());
    cmeshflux->AutoBuild();    

    return cmeshflux;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * CreateCMeshHCurl(TPZGeoMesh* gmesh, int pOrder, const int volId, HCurlFamily hcurlfam){
    //Create flux mesh
    TPZCompMesh *cmeshflux = new TPZCompMesh(gmesh);
    cmeshflux->SetDimModel(gmesh->Dimension());
    cmeshflux->SetDefaultOrder(pOrder);
    
    TPZL2ProjectionHCurl<> *mat = new TPZL2ProjectionHCurl<>(volId,gmesh->Dimension());
    mat->SetForcingFunction(forcefunction,4);
    cmeshflux->InsertMaterialObject(mat);

    cmeshflux->ApproxSpace().SetHCurlFamily(hcurlfam);
    cmeshflux->ApproxSpace().SetAllCreateFunctionsHCurl(gmesh->Dimension());
    cmeshflux->AutoBuild();

    return cmeshflux;
}


void CheckElementInterfaces(TPZCompMesh *cmesh){

    TPZGeoMesh *gmesh = cmesh->Reference();
    int fDim = cmesh->Dimension();
    
    //Prints computational mesh properties
#ifdef  DEBUGTEST
    std::string txt = "CMesh.txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);
#endif
    for (auto gel : gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->Dimension() != fDim) continue;

        //Loop over all neighbours
        int nSides = gel->NSides();
        int nEdges = gel->NSides(fDim-1);
        int nCorner = gel->NCornerNodes();
        for (int iSide = 0; iSide < nEdges; iSide++){
            TPZGeoElSide neigh = gel->Neighbour(nSides-1-nEdges+iSide);
            TPZGeoElSide gelside(gel,nSides-1-nEdges+iSide);

            if (neigh.Element()->Dimension() != fDim) continue;

            // std::cout << "Side Nodes = " ;
            // for (int i = 0; i < gelside.NSideNodes(); i++){
            //     std::cout << gelside.SideNodeIndex(i) << " ";
            // }
            // std::cout << std::endl;

            // std::cout << "Neigh Nodes = " ;
            // for (int i = 0; i < neigh.NSideNodes(); i++){
            //     std::cout << neigh.SideNodeIndex(i) << " ";
            // }
            // std::cout << std::endl;
            // std::cout << std::endl;

            CheckSideOrientation(gelside,neigh);
            

            // TPZCheckRestraint check = new TPZCheckRestraint(neigh.Reference(),gelside.Reference());

            // std::cout << "Side = " << iSide << " " << gelside << std::endl;
            // std::cout << "Neigh = " << iSide << " " << neigh << std::endl;
        }
        

        
    }

}


void CheckSideOrientation(TPZGeoElSide &thisGeoSide, TPZGeoElSide &largeGeoSide){

    TPZManVector<int64_t> thisSideNodes(thisGeoSide.NSideNodes(),0);
    TPZManVector<int64_t> largeSideNodes(largeGeoSide.NSideNodes(),0);
    
    // thisGeoSide.SideNodeIndex(i)
    for (int i = 0; i < thisGeoSide.NSideNodes(); i++) thisSideNodes[i] = thisGeoSide.SideNodeIndex(i); 
    for (int i = 0; i < largeGeoSide.NSideNodes(); i++) largeSideNodes[i] = largeGeoSide.SideNodeIndex(i); 

    // REAL orient = 1.;
    // if (largeGeoSide.NSideNodes() != thisGeoSide.NSideNodes()) DebugStop();
    // for (int i = 0; i < thisGeoSide.NSideNodes(); i++)
    // {
    //     if (thisSideNodes[i] != largeSideNodes[i]){
    //         orient = -1;
    //         break;
    //     }
        
    // }
    // std::cout << "thisGeoSide Nodes = " ;
    // for (int i = 0; i < thisGeoSide.NSideNodes(); i++){
    //     std::cout << thisGeoSide.SideNodeIndex(i) << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "largeGeoSide Nodes = " ;
    // for (int i = 0; i < largeGeoSide.NSideNodes(); i++){
    //     std::cout << largeGeoSide.SideNodeIndex(i) << " ";
    // }
    // std::cout << std::endl;

    TPZVec<REAL> normalLarge,normalThis;
    auto neighTransf = thisGeoSide.NeighbourSideTransform(largeGeoSide);
    TPZManVector<REAL, 3> neighXi(largeGeoSide.Dimension(), 0);
    TPZManVector<REAL,3> xiSide(thisGeoSide.Dimension(),0);
    neighTransf.Apply(xiSide, neighXi);
    largeGeoSide.Normal(neighXi,normalLarge);
    thisGeoSide.Normal(xiSide,normalThis);
    REAL det;
    TPZFMatrix<REAL> inv;
    neighTransf.Mult().DeterminantInverse(det,inv);

    // std::cout << "Normal Large = " << normalLarge << std::endl;
    // std::cout << "Normal This = " << normalThis << std::endl;
    // std::cout << "Determinant = " << det << std::endl << std::endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

#undef EXPORTVTK
#undef DEBUGTEST
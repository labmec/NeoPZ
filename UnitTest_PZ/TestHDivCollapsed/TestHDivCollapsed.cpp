// ----- PZ includes -----
#include <TPZGenGrid3D.h>
#include <TPZGenGrid2D.h>
#include <pzgmesh.h>
#include <TPZVTKGeoMesh.h>
#include <pzgeoquad.h>
#include <tpzgeoelrefpattern.h>
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzshapequad.h>
#include <pzshapelinear.h>
#include <TPZCompElHDivCollapsed.h>
#include <TPZHybridizeHDiv.h>

#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <pzfstrmatrix.h>

#include <TPZNullMaterial.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <DarcyFlow/TPZMixedDarcyFractureFlow.h>

#include <pzlog.h>

// ----- Unit test includes -----
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
// ----- Run tests with or without main -----
//#define RUNWITHMAIN

// ----- Functions -----
void TestHdivCollapsed(const bool& is3D, const bool& isRefMesh, const bool& isLinPVar, const bool& isFracIntersect);
TPZGeoMesh *Create3DGeoMesh(const bool& isRefMesh);
TPZGeoMesh *Create3DGeoMeshIntersectFrac(const bool& isRefMesh);
TPZGeoMesh *Create2DGeoMesh(const bool& isRefMesh, const bool& isFracIntersect);
TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
void CreateFractureHDivCollapsedEl(TPZCompMesh* cmesh);
void SplitConnectsAtInterface(TPZCompElSide& compside);
TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *>& meshvector,
                                          TPZGeoMesh * gmesh, const bool& isLinPVar, TPZHybridizeHDiv* hybridizer);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *>& meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh,
                             const std::string &plotfile);
const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname);
void HybridizeIntersections(TPZVec<TPZCompMesh *>& meshvec_Hybrid, TPZHybridizeHDiv *hybridizer);
void CreateIntersectionInterfaceElements(TPZMultiphysicsCompMesh* cmesh, TPZHybridizeHDiv *hybridizer);

using namespace std;

// Materials

enum EMatid {ENone,
    EPressureFracBnd,
    EDomain,
    EFaceBCPressure, /*Pressure at domain boundary*/
    ENoFlux, /*No flux in domain boundary*/
    ENoFluxFracBorder, /*No flux at fracture boundary*/
    EFracture,
    EIntersection, /*Used mostly for finding where hybridization should occur*/
    EKd2 /*Pressure loss at intersection*/
};

#ifndef RUNWITHMAIN

// ----- Test cases -----
// ---- Test 0 ----
TEST_CASE("2D_1_frac_element","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = false;
    const bool isLinP = false;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 1 ----
TEST_CASE("2D_uniformly_refined_mesh","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = true;
    const bool isLinP = false;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 2 ----
TEST_CASE("3D_1_frac_element","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = false;
    const bool isLinP = false;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 3 ----
TEST_CASE("3D_uniformly_refined_mesh","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = true;
    const bool isLinP = false;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 4 ----
TEST_CASE("2D_1_frac_element_linP","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = false;
    const bool isLinP = true;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 5 ----
TEST_CASE("2D_uniformly_refined_mesh_linP","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = true;
    const bool isLinP = true;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 6 ----
TEST_CASE("3D_1_frac_element_linP","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = false;
    const bool isLinP = true;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 7 ----
TEST_CASE("3D_uniformly_refined_mesh_linP","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = true;
    const bool isLinP = true;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 8 ----
TEST_CASE("2D_2_frac_intersect","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = false;
    const bool isLinP = false;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 9 ----
TEST_CASE("2D_2_frac_intersect_uniformly_refined","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = true;
    const bool isLinP = false;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 10 ----
TEST_CASE("3D_2_frac_intersect","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = false;
    const bool isLinP = false;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 11 ----
TEST_CASE("3D_2_frac_intersect_uniformly_refined","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = true;
    const bool isLinP = false;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 12 ----
TEST_CASE("2D_2_frac_intersect_linP","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = false;
    const bool isLinP = true;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 13 ----
TEST_CASE("2D_2_frac_intersect_uniformly_refined_linP","[hdivcollapsed]"){
    const bool is3D = false;
    const bool isRefMesh = true;
    const bool isLinP = true;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 14 ----
TEST_CASE("3D_2_frac_intersect_linP","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = false;
    const bool isLinP = true;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}
// ---- Test 15 ----
TEST_CASE("3D_2_frac_intersect_uniformly_refined_linP","[hdivcollapsed]"){
    const bool is3D = true;
    const bool isRefMesh = true;
    const bool isLinP = true;
    const bool isFracIntersect = true;
    TestHdivCollapsed(is3D,isRefMesh,isLinP,isFracIntersect);
}

#else

// Left in case needs some serious debugging. Catch2 does not stop in debugstops in xcode
int main(int argc, char* argv[]){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    const bool isRefMesh = false;
    const bool is3D = true;
    const bool isLinPVar = true;
    const bool isFracIntersect = false;
    TestHdivCollapsed(is3D,isRefMesh,isLinPVar,isFracIntersect);

    return 0;
}

#endif

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
void TestHdivCollapsed(const bool& is3D, const bool& isRefMesh, const bool& isLinPVar, const bool& isFracIntersect){
    
    // ----- porder and header -----
    constexpr int pOrder{1};
    std::stringstream header;
    const int dim = is3D ? 3 : 2;
    header << "Dim: " << dim;
    if(isRefMesh) header << "  Ref";
    if(isFracIntersect) header << "  Intersect";
    if(isLinPVar) header << "  LinP";
    std::cout << "\n ============ Running Problem " << header.str() << " ============\n" << std::endl;
    
    TPZGeoMesh* gmesh = nullptr;
    if (is3D){
        if (isFracIntersect) {
            gmesh = Create3DGeoMeshIntersectFrac(isRefMesh);
        }
        else {
            gmesh = Create3DGeoMesh(isRefMesh);
        }
    }
    else{
        gmesh = Create2DGeoMesh(isRefMesh,isFracIntersect);
    }
        
    
    std::stringstream gfileroot,gfilename;
    gfileroot << "GeoMesh" << gmesh->Dimension() << "D_";
    if(isRefMesh) gfileroot << "Ref";
    if(isFracIntersect) gfileroot << "Intersect";
    gfilename << gfileroot.str() <<  ".vtk";
    // ----- Print gmesh -----
    std::ofstream outgmesh(gfilename.str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outgmesh);
    
    // ----- Create Flux mesh -----
    TPZCompMesh * cmeshflux = FluxCMesh(dim,pOrder,gmesh);
    {
        std::stringstream fluxmeshname;
        fluxmeshname << "Flux" << gfileroot.str() << ".txt";
        std::ofstream out(fluxmeshname.str());
        cmeshflux->Print(out);
    }
    // ----- Pressure mesh -----
    TPZCompMesh * cmeshpressure = PressureCMesh(dim,pOrder,gmesh);
    {
        std::stringstream fluxmeshname;
        fluxmeshname << "Pressure" << gfileroot.str() << ".txt";
        std::ofstream out(fluxmeshname.str());
        cmeshpressure->Print(out);
    }

    // ----- Multiphysics mesh -----
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshflux;
    meshvector[1] = cmeshpressure;

    // ----- Hybridize fracture intersection if any -----
    TPZHybridizeHDiv* hybridizer = nullptr;
    if (isFracIntersect) {
        hybridizer = new TPZHybridizeHDiv(meshvector);
        HybridizeIntersections(meshvector,hybridizer);
        
        std::stringstream fluxmeshname;
        fluxmeshname << "FluxHybridized" << gfileroot.str() << ".txt";
        std::ofstream out(fluxmeshname.str());
        cmeshflux->Print(out);
    }
    
    TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,meshvector,gmesh,isLinPVar,hybridizer);
    std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;

    // ----- Solve system -----
    TPZLinearAnalysis an(cmesh,RenumType::ENone);
    SolveProblemDirect(an,cmesh);

    // ----- Print results -----
    std::stringstream pfilename;
    pfilename << "Result" << gmesh->Dimension();
    if(isRefMesh) pfilename << "Ref";
    if(isFracIntersect) pfilename << "Intersect";
    if(isLinPVar) pfilename << "LinP";
    std::stringstream pfilenamefrac;
    pfilenamefrac << pfilename.str();
    pfilename << ".vtk";
    pfilenamefrac << "frac.vtk";

    PrintResultsMultiphysic(dim,meshvector,an,cmesh,pfilename.str());
    PrintResultsMultiphysic(dim-1,meshvector,an,cmesh,pfilenamefrac.str());

    // ----- Compute integral of pressure and flux over domain and compare with analytical solution -----
    const std::string pvarname = "Pressure";
    const STATE integratedpressure = ComputeIntegralOverDomain(cmesh,pvarname);
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;
    
    const std::string qvarname = "Flux";
    STATE integratedflux = ComputeIntegralOverDomain(cmesh,qvarname);
    if (fabs(integratedflux) < 1.e-14 ) integratedflux = 0.; // to make Catch::Approx(0.) work
    std::cout << "\nintegral of flux  = " << integratedflux << std::endl;
    
    // ----- Comparing with analytical solution -----
    // For 3d:
    // Domain volume is 2*2*2=8. If p cte: 1*8 = 8. If p varies linearly from 2 to 0: ((2-0)/2) * 8 = 8
    // For 2d:
    // Domain volume is 2*2=4. If p cte: 1*4 = 4. If p varies linearly from 2 to 0: ((2-0)/2) * 4 = 8
#ifndef RUNWITHMAIN
    if (is3D) {
        REQUIRE( integratedpressure == Catch::Approx( 8.0 ) ); // Approx is from catch2 lib
        if (isLinPVar)
            REQUIRE( integratedflux == Catch::Approx( 8./3. ) ); // Approx is from catch2 lib
        else
            REQUIRE( integratedflux == Catch::Approx( 0.) ); // Approx is from catch2 lib
    }
    else{
        REQUIRE( integratedpressure == Catch::Approx( 4.0 ) ); // Approx is from catch2 lib
        if (isLinPVar)
            REQUIRE( integratedflux == Catch::Approx( 4./3. ) ); // Approx is from catch2 lib
        else
            REQUIRE( integratedflux == Catch::Approx( 0.) ); // Approx is from catch2 lib
    }
#endif
    
    delete gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    // ===> Creating compmesh
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    // ===> Domain material
    const int nstate = 1;
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain,dim,nstate);
    cmesh->InsertMaterialObject(mat);

    // ===> Boundary conditions
    TPZNullMaterial<> *matbc = new TPZNullMaterial<>(EFaceBCPressure,dim-1,nstate);
    cmesh->InsertMaterialObject(matbc);
    TPZNullMaterial<> *matbcnf = new TPZNullMaterial<>(ENoFlux,dim-1,nstate);
    cmesh->InsertMaterialObject(matbcnf);
    
    TPZNullMaterial<> *matbcnfpt = new TPZNullMaterial<>(ENoFluxFracBorder,dim-2,nstate);
    cmesh->InsertMaterialObject(matbcnfpt);

    
    
    
    // ===> Fracture material
    TPZNullMaterial<> *matfrac = new TPZNullMaterial<>(EFracture,dim-1,nstate);
    cmesh->InsertMaterialObject(matfrac);
    
    // ===> Intersection material for transverse flux
    TPZNullMaterial<> *matintersect = new TPZNullMaterial<>(EKd2,dim-2,nstate);
    cmesh->InsertMaterialObject(matintersect);
    
    // ===> Fracture boundary conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,1.);
    TPZNullMaterial<> *matbcfrac = new TPZNullMaterial<>(EPressureFracBnd,dim-2,nstate);
    cmesh->InsertMaterialObject(matbcfrac);
    
    // ===> Creating space for 3D elements
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ApproxSpace().CreateDisconnectedElements(false); // we need to disconnect by hand at the fracture location later
    cmesh->SetDefaultOrder(pOrder);
    std::set<int> buildmatids;
    buildmatids.insert(EDomain);
    cmesh->AutoBuild(buildmatids);
    
    // Print flux mesh
    std::ofstream myfile("FluxMeshOnlyVol.txt");
    cmesh->Print(myfile);
    
    // ===> Creating fracture element
    // split the volumetric elements in the flux mesh
    CreateFractureHDivCollapsedEl(cmesh);
    cmesh->CleanUpUnconnectedNodes();
    
    // ===> Create BCs for fracture
    buildmatids.clear();
    buildmatids.insert(EPressureFracBnd);
    buildmatids.insert(ENoFluxFracBorder);
    gmesh->ResetReference();
    for (auto cel : cmesh->ElementVec()) {
        if (!cel) {
            continue;
        }
        if (cel->Reference()->MaterialId() == EFracture) {
            cel->LoadElementReference();
        }
    }
    cmesh->SetDimModel(dim-1);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild(buildmatids);
    // adjust the side orientation of the boundary sides
    for (auto cel : cmesh->ElementVec()) {
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(buildmatids.find(matid) == buildmatids.end()) continue;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        intel->SetSideOrient(gel->NSides()-1, 1);
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour(gelside.Neighbour());
        while(neighbour != gelside)
        {
            TPZGeoEl *neighgel = neighbour.Element();
            TPZCompEl *neighcel = neighgel->Reference();
            if(neighcel && neighgel->Dimension() == gel->Dimension()+1)
            {
                TPZInterpolationSpace *neighintel = dynamic_cast<TPZInterpolationSpace *>(neighcel);
                neighintel->SetSideOrient(neighbour.Side(), 1);
            }
            neighbour = neighbour.Neighbour();
        }
    }

    
    
    // ===> Create BCs for 3D domain
    buildmatids.clear();
    buildmatids.insert(EFaceBCPressure);
    buildmatids.insert(ENoFlux);
    gmesh->ResetReference();
    for (auto cel : cmesh->ElementVec()) {
        if (!cel) {
            continue;
        }
        if (cel->Reference()->MaterialId() == EDomain) {
            cel->LoadElementReference();
        }
    }
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild(buildmatids);
    cmesh->CleanUpUnconnectedNodes();
    
    // ===> Fix blocks
    cmesh->SetDimModel(dim);
    cmesh->InitializeBlock();

    // Print flux mesh
//    std::ofstream myfile2("FluxMeshEnd.txt");
//    cmesh->Print(myfile2);

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateFractureHDivCollapsedEl(TPZCompMesh* cmesh) {
    
    // ===> Create the fracture hdivcollapsed elements
    // These should be unconnected with the 3D elements so we reset references in the gmesh
    // They are, however, connected (with connects) between themselves
    TPZGeoMesh* gmesh = cmesh->Reference();
    gmesh->ResetReference(); // So it does not try to use connects of neighbor elements
    const int gmeshdim = gmesh->Dimension();
    cmesh->SetDimModel(gmeshdim-1);
    for(auto gel : gmesh->ElementVec()) {
        const int matid = gel->MaterialId();
        if (matid != EFracture) continue;
        const int hassubel = gel->HasSubElement();
        if (hassubel) { // the mesh can be uniformly refined
            continue;
        }
        TPZInterpolationSpace* hdivcollapsed = nullptr;
        if (gmeshdim == 2){
            hdivcollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeLinear>(*cmesh,gel);
        }
        else {
            hdivcollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>(*cmesh,gel);
        }
        
    }
    
    // ===> Set HdivCollapsed elements connection with adjacent 3d elements
    // This is done through top and bottom connects of the hdivcollapsed el.
    // Since bottom connect's default sideorient is negative one (-1),
    // we should connect it with the 3D element side that has sideorient positive one (+1)
    cmesh->LoadReferences(); // So the hdivcollapsed element can see the 3D the flux mesh adjacent elements to it
    
    // if there are superposed fracture elements then first connects need to be created to link the superposed
    // fracture elements
    // then the extremes need to link to the 3D elements
    for(auto cel : cmesh->ElementVec()) {
        if (!cel) continue;
        TPZGeoEl* gel = cel->Reference();
        const int matid = gel->MaterialId();
        if (matid != EFracture) continue;
        int64_t index;
        TPZInterpolationSpace* hdivcollapsed = nullptr;
        if (gmeshdim == 2) {
            hdivcollapsed = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeLinear>* >(cel);
        }
        else{
//            TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>* hdivcollapsed = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>* >(cel);
            
            // this code needs to be generalised to consider triangular fracture elements
            hdivcollapsed = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>* >(cel);
        }
        
        if (!hdivcollapsed) {
            DebugStop();
        }
        
        int nconnects = hdivcollapsed->NConnects();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neigh = gelside.Neighbour();
        int icon = 0;
        // loop over all neighbours of the fracture element
        while(neigh != gelside){
            TPZGeoEl* gelnei = neigh.Element();
            if (!gelnei || gelnei->Dimension() != gmeshdim || gelnei->MaterialId() != EDomain) {
                neigh++;
                continue;
            }
            // we found a volumetric element
            TPZCompElSide compside = TPZCompElSide(gelnei->Reference(), neigh.Side());
            if (icon == 0){
                // split the connect shared by two adjacent 3d elements if it is the first time
                SplitConnectsAtInterface(compside);
            }
            
            const int conindex = compside.ConnectIndex();
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(compside.Element());
            const int sideorient = intel->GetSideOrient(compside.Side());
            if (sideorient == 1) { // this is a bottom connection (bottom is -1)
                hdivcollapsed->SetConnectIndex(nconnects, conindex);
                icon++;
            }
            else if (sideorient == -1){
                hdivcollapsed->SetConnectIndex(nconnects+1, conindex);
                icon++;
            }
            else{
                DebugStop();
            }
            neigh++;
        }
        if (icon != 2) {
            DebugStop();
        }
    }
    
    cmesh->ExpandSolution();
//    std::ofstream out("cmeshaftercollapsed.txt");
//    cmesh->Print(out);
    
    cmesh->SetDimModel(gmeshdim);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void SplitConnectsAtInterface(TPZCompElSide& compside) {
    
    TPZCompMesh* fluxmesh = compside.Element()->Mesh();
    const int gmeshdim = fluxmesh->Reference()->Dimension();
    
    // ===> Find 3D element neighbor
    TPZCompElSide compsideright;
    TPZGeoElSide gleft(compside.Reference());
    TPZGeoElSide neigh = gleft.Neighbour();
    int icon = 0;
    while(neigh != gleft){
        TPZGeoEl* gelnei = neigh.Element();
        if (!gelnei || gelnei->Dimension() != gmeshdim || gelnei->MaterialId() != EDomain) {
            neigh++;
            continue;
        }
        compsideright = TPZCompElSide(gelnei->Reference(), neigh.Side());
        break;
    }
    if (!compsideright.Element()) {
        DebugStop(); // Could not find neighbor element!
    }
    
    compside.SplitConnect(compsideright);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    // ===> Creating compmesh
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);

    // ===> Domain material
    const int nstate = 1;
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain,dim,nstate);
    cmesh->InsertMaterialObject(mat);

    // ===> Fracture material
    TPZNullMaterial<> *matfrac = new TPZNullMaterial<>(EFracture,dim-1,nstate);
    cmesh->InsertMaterialObject(matfrac);
    
    // ===> P order distinction
    if(pOrder == 0)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
    // ===> Creating space for 3D elements
    std::set<int> buildmatids;
    buildmatids.insert(EDomain);
    buildmatids.insert(EFracture);
    cmesh->AutoBuild(buildmatids);
        
    // ===> Set lagrange multiplier for order of assembly
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    // Print pressure mesh
//    std::ofstream myfile("PressureMesh.txt");
//    cmesh->Print(myfile);

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
auto exactSolFrac = [](const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& gradU){
    const auto& x = loc[0];
    const auto& y = loc[1];
    const auto& z = loc[2];
    const REAL gap = 1./3.;
    if(y < -1.e-2){
        u[0] = 2.0 - (y+1)*gap;
    }
    else if(y > 1.e-2){
        u[0] = gap - y*gap;
    }
    else{
        u[0] = 1. - y;
    }
    // Not using derivatives
};
auto exactSolLinP = [](const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& gradU){
    const auto& x = loc[0];
    const auto& y = loc[1];
    const auto& z = loc[2];
    const REAL gap = 1./3.;
    if(y < 0){
        u[0] = 2.0 - (y+1)*gap;
    }
    else{
        u[0] = gap - y*gap;
    }
    // Not using derivatives
};
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *>& meshvector,TPZGeoMesh * gmesh,
                                          const bool& isLinPVar, TPZHybridizeHDiv* hybridizer)
{
    // ===> Creating mp cmesh
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    // ===> Domain mat
    auto mat = new TPZMixedDarcyFlow(EDomain, dim);
    mat->SetConstantPermeability(1.);
    cmesh->InsertMaterialObject(mat);
    
    // ===> Fracture mat
    auto matfrac = new TPZMixedDarcyFractureFlow(EFracture, dim-1);
    matfrac->SetConstantPermeability(1.0);
    cmesh->InsertMaterialObject(matfrac);
            
    // ===> Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,1.);
    // domain bcs
    auto * BCond0 = mat->CreateBC(mat, EFaceBCPressure, 0, val1, val2);
    if (isLinPVar)
        BCond0->SetForcingFunctionBC(exactSolLinP,2);
    cmesh->InsertMaterialObject(BCond0);
    
    TPZManVector<STATE> val2n(1,0.);
    auto * BCondNoFlux = mat->CreateBC(mat,ENoFlux, 1, val1, val2n);
    cmesh->InsertMaterialObject(BCondNoFlux);
    auto * BCondNoFluxPt = mat->CreateBC(mat,ENoFluxFracBorder, 1, val1, val2n);
    cmesh->InsertMaterialObject(BCondNoFluxPt);
    
    // ===> Intersection material for transverse flux
    TPZFMatrix<STATE> val1inter(1,1,0.5);
    TPZManVector<STATE> val2inter(1,0.);
    auto * BCondintersect = mat->CreateBC(mat, EKd2, 2, val1inter, val2inter);
    cmesh->InsertMaterialObject(BCondintersect);
    
    // frac bcs
    auto * BCond1 = mat->CreateBC(mat, EPressureFracBnd, 0, val1, val2);
//    auto * BCond1 = mat->CreateBC(mat, EPressureFracBnd, 1, val1, val2n);
    if (isLinPVar)
        BCond1->SetForcingFunctionBC(exactSolFrac,2);
    cmesh->InsertMaterialObject(BCond1);
    
    // ===> Materials for hybridizing intersection between fractures
    if (hybridizer) {
        hybridizer->InsertPeriferalMaterialObjects(cmesh);
    }
    
    // ===> Build space
    TPZManVector<MSpaceConfig> active(2,ETestTrial);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    
    if (hybridizer){
        CreateIntersectionInterfaceElements(cmesh,hybridizer);
    }

//    Prints Multiphysics mesh
//    std::ofstream myfile("MultiPhysicsMesh.txt");
//    cmesh->Print(myfile);

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    constexpr int nThreads{0};
//    TPZFStructMatrix<STATE> matskl(cmesh);
//    TPZSkylineStructMatrix<STATE> matskl(cmesh);
#if defined (USING_MKL) || (USING_EIGEN)
    TPZSpStructMatrix<STATE> matskl(cmesh);
#else
    TPZSkylineNSymStructMatrix<STATE> matskl(cmesh);
#endif
//    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    an.Assemble();
    
//  {
//    std::ofstream outmat("mat.nb");
//    TPZMatrixSolver<STATE>* matsol = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
//    matsol->Matrix()->Print("singmat=",outmat,EMathematicaInput);
//    std::ofstream outrhs("rhs.nb");
//      TPZFMatrix<STATE> rhs = an.Rhs();
//    rhs.Print("rhs=",outrhs,EMathematicaInput);
//  }
    
    ///solves the system
    an.Solve();
//    {
//        std::ofstream outsol("sol.nb");
//        TPZFMatrix<STATE> sol = an.Solution();
//        sol.Print("sol = ",outsol,EMathematicaInput);
//    }
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *>& meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh,
                             const std::string &plotfile)
{
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    TPZManVector<std::string,2> scalnames(1), vecnames(1);
    
    scalnames[0] = "Pressure";
    vecnames[0]= "Flux";
    
    int div = 0;
//    std::string plotfile = "PostProcess1Frac.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *Create3DGeoMesh(const bool& isRefMesh) {
    
    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,2,1};
    const MMeshType elType = MMeshType::EHexahedral;
    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    TPZGeoMesh *gmesh = gen3d.BuildVolumetricElements(EDomain);
    
    // ----- Fracture element and bcs -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {2,3};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressureFracBnd,*gmesh,index);
    nodesId = {3,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressureFracBnd,*gmesh,index);
    nodesId = {9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressureFracBnd,*gmesh,index);
    nodesId = {8,2};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressureFracBnd,*gmesh,index);
    TPZManVector<int64_t,4> nodesIdVec = {2,3,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,EFracture,*gmesh,index);
    gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure);
    gmesh->BuildConnectivity();
    
    if (isRefMesh) {
        gRefDBase.InitializeUniformRefPattern(ECube);
        gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
        gRefDBase.InitializeUniformRefPattern(EOned);
        for (auto gel : gmesh->ElementVec()){
            TPZManVector<TPZGeoEl*,10> children;
            gel->Divide(children);
        }
        gmesh->BuildConnectivity();
    }
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *Create3DGeoMeshIntersectFrac(const bool& isRefMesh) {
    
    TPZGeoMesh* gmesh = nullptr;
    
    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,2,2};
    const MMeshType elType = MMeshType::EHexahedral;

    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EDomain);

    // ----- Fracture element -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {6,7};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressureFracBnd,*gmesh,index);
    nodesId = {7,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {9,11};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressureFracBnd,*gmesh,index);
    nodesId = {10,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {8,6};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);

    nodesId = {9,15};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {14,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {8,2};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {2,3};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {3,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,ENoFluxFracBorder,*gmesh,index);

    nodesId = {8,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EIntersection,*gmesh,index);
    
    TPZManVector<int64_t,4> nodesIdVec = {6,7,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,EFracture,*gmesh,index);
    nodesIdVec = {8,9,11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,EFracture,*gmesh,index);
    nodesIdVec = {2,3,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,EFracture,*gmesh,index);
    nodesIdVec = {8,9,15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,EFracture,*gmesh,index);
    
    // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
//    gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure);
    gmesh = gen3d.BuildBoundaryElements(ENoFlux, ENoFlux, EFaceBCPressure, ENoFlux, EFaceBCPressure, ENoFlux);
    
    gmesh->BuildConnectivity();
    
    if (isRefMesh) {
        gRefDBase.InitializeUniformRefPattern(ECube);
        gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
        gRefDBase.InitializeUniformRefPattern(EOned);
        for (auto gel : gmesh->ElementVec()){
            TPZManVector<TPZGeoEl*,10> children;
            gel->Divide(children);
        }
        gmesh->BuildConnectivity();
    }
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *Create2DGeoMesh(const bool& isRefMesh, const bool& isFracIntersect) {
    
    // ----- Create Geo Mesh -----
    const TPZManVector<REAL,2> minX = {-1.,-1.};
    const TPZManVector<REAL,2> maxX = {1.,1.};
    const TPZManVector<int,2> nelDiv = {2,2};
    TPZGenGrid2D gen2d(nelDiv,minX,maxX);
    const MMeshType elType = MMeshType::EQuadrilateral;
    gen2d.SetElementType(elType);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen2d.Read(gmesh,EDomain);
//    for (int iside = 4; iside < 8; iside++) {
//        gen2d.SetBC(gmesh, iside, EFaceBCPressure);
//    }
    gen2d.SetBC(gmesh, 4, EFaceBCPressure);
    gen2d.SetBC(gmesh, 5, ENoFlux);
    gen2d.SetBC(gmesh, 6, EFaceBCPressure);
    gen2d.SetBC(gmesh, 7, ENoFlux);

    // ----- Fracture element and bcs -----
    int64_t index;
    // ===> BCs
    TPZManVector<int64_t,1> nodesId = {3};
    new TPZGeoElRefPattern<pzgeom::TPZGeoPoint>(nodesId,ENoFluxFracBorder,*gmesh,index);
    nodesId = {5};
    new TPZGeoElRefPattern<pzgeom::TPZGeoPoint>(nodesId,ENoFluxFracBorder,*gmesh,index);
    if (isFracIntersect) {
        nodesId = {1};
        new TPZGeoElRefPattern<pzgeom::TPZGeoPoint>(nodesId,EPressureFracBnd,*gmesh,index);
        nodesId = {7};
        new TPZGeoElRefPattern<pzgeom::TPZGeoPoint>(nodesId,EPressureFracBnd,*gmesh,index);
    }
    // ===> Fracture el
    TPZManVector<int64_t,2> nodesIdVec = {3,4};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec,EFracture,*gmesh,index);
    nodesIdVec = {4,5};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec,EFracture,*gmesh,index);
    if (isFracIntersect) {
        nodesIdVec = {1,4};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec,EFracture,*gmesh,index);
        nodesIdVec = {4,7};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec,EFracture,*gmesh,index);

        // ===> Intersection el
        nodesId = {4};
        new TPZGeoElRefPattern<pzgeom::TPZGeoPoint>(nodesId,EIntersection,*gmesh,index);
    }

    gmesh->BuildConnectivity();
    
    if (isRefMesh) {
        gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
        gRefDBase.InitializeUniformRefPattern(EOned);
        for (auto gel : gmesh->ElementVec()){
            TPZManVector<TPZGeoEl*,10> children;
            gel->Divide(children);
        }
        gmesh->BuildConnectivity();
    }
//    {
//        std::ofstream out("gmesh2d.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
//    }
    
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname) {
    std::set<int> matids;
    matids.insert(EDomain);
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences(); // compute integral in the multiphysics mesh
    TPZVec<STATE> vecint = cmesh->Integrate(varname, matids);
    if ((varname == "Pressure" && vecint.size() != 1) ||
        (varname == "Flux" && vecint.size() != 3)){
        DebugStop();
    }
    if (varname == "Pressure")
        return vecint[0];
    else if (varname == "Flux")
        return vecint[1];
    else
        DebugStop();
    
    return -100000; // default value so compiler does not complain
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void HybridizeIntersections(TPZVec<TPZCompMesh *>& meshvec_Hybrid, TPZHybridizeHDiv *hybridizer) {
    
    if (!hybridizer){
        DebugStop();
    }
    
    hybridizer->fHDivWrapMatid = EKd2;
//  hybridizer->fIdToHybridize = EFracture;
    hybridizer->fIdsToHybridize.insert(EFracture);
    
    // ===> Initializing variables
    TPZCompMesh* fluxmesh = meshvec_Hybrid[0];
    TPZGeoMesh* gmesh = fluxmesh->Reference();
    const int dim = gmesh->Dimension();
    int dimfrac = dim-1; // Fractures are always 2d in this case?
    
    // ===> Insert hybridization materials in flux and pressure meshes
    fluxmesh->SetDimModel(dimfrac); // to create hdivbounds around intersection
    hybridizer->InsertPeriferalMaterialObjects(meshvec_Hybrid);
    
    // ===> Find intersection geoels
    fluxmesh->LoadReferences();
    fluxmesh->SetAllCreateFunctionsHDiv();
    for (auto gel : gmesh->ElementVec()) {
        const int gelmatid = gel->MaterialId();
        if (gelmatid != EIntersection || gel->HasSubElement()) {
            continue;
        }
        if (gel->Dimension() != dim-2) {
            DebugStop();
        }
        
        // Search for first neighbor that that is domain
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        while(neigh != gelside){
            TPZGeoEl* gelneigh = neigh.Element();
            if (gelneigh->HasSubElement()) {
                neigh++;
                continue;
            }
            int neighmatid = gelneigh->MaterialId();
            int neighdim = gelneigh->Dimension();
            
            if (neighmatid == EFracture && neighdim == dimfrac) {
                cout << "\nElement with ID " << gel->Id() << " and index " << gel->Index() << " is an intersection element" << endl;
                cout << "===> Trying to split the connects of the flux mesh and create pressure element..." << endl;
                TPZCompEl* celneigh = gelneigh->Reference();
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (celneigh);
                if (!intel)
                    DebugStop();
                
                const int side = neigh.Side();
                TPZCompElSide celsideleft(intel, side);
                bool isNewInterface = hybridizer->HybridizeInterface(celsideleft,intel,side,meshvec_Hybrid,false/*isIntersectEnd*/);
                if (isNewInterface) {
                    cout << "=====> Connects splitted succesfuly!" << endl << endl;;
                    break;
                }
                else{
                    DebugStop();
                }
            }
            neigh = neigh.Neighbour();
        } // while
    }
    
    fluxmesh->SetDimModel(dim);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateIntersectionInterfaceElements(TPZMultiphysicsCompMesh* cmesh, TPZHybridizeHDiv *hybridizer) {
    TPZCompMesh* cmeshpressure = cmesh->MeshVector()[1];
    cmesh->Reference()->ResetReference();
    cmeshpressure->LoadReferences();
    const int lagrangematid = hybridizer->lagrangeInterfaceMatId();
    const int lagrangematidend = hybridizer->lagrangeInterfaceEndMatId();
    for (auto cel : cmeshpressure->ElementVec()) {
        if(!cel) continue;
        const int celmatid = cel->Material()->Id();
        if (celmatid != lagrangematid && celmatid != lagrangematidend) {
            continue;
        }
        TPZGeoEl* gel = cel->Reference();
        hybridizer->CreateInterfaceElementsForGeoEl(cmesh, cmesh->MeshVector(), gel);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

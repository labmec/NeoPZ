//
//  TestKernelHDiv.cpp
//  PZ
//
//  Created by Jeferson Fernandes on 16/01/22.
//  Copyright 2022 UNICAMP. All rights reserved.
//

#include <TPZGeoMeshTools.h>
#include <pzcmesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include "TPZMultiphysicsCompMesh.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"

using namespace pzshape;

#include <pzfstrmatrix.h>
#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h>
#include <TPZGmshReader.h>
#include "pzvec.h"
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include "pzmultiphysicscompel.h"
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"// for Hdiv problem
#include <TPZVTKGeoMesh.h>
#include "TPZAnalyticSolution.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElDisc.h"
#include "TPZCompElH1.h"
#include "TPZCompElKernelHDiv.h"
#include "TPZCompElKernelHDiv3D.h"

#include "TPZHCurlEquationFilter.h"

#include "pznoderep.h"
#include "pzlog.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#undef PZ_KERNELHDIV_DEBUG



enum MShapeType {EHDivKernel, EHDivConstant};
enum BCType {Dirichlet, Neumann};


auto exactSolKernel3D = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = x+y+z;
    gradU(0,0) = 1.;
    gradU(1,0) = 1.;
    gradU(2,0) = 1.;
};

auto exactSolKernel2D = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = x+y;
    gradU(0,0) = 1.;
    gradU(1,0) = 1.;
    gradU(2,0) = 0.;
};

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];

    // //Nabla u = 1
    u[0] = 1.;
};


auto exactSolConstant3D = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = (x*x+y*y+z*z)/6.;
    gradU(0,0) = x/3.;
    gradU(1,0) = y/3.;
    gradU(2,0) = z/3.;
};

auto exactSolConstant2D = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = 0.25*(x*x+y*y);
    gradU(0,0) = 0.5*x;
    gradU(1,0) = 0.5*y;
    gradU(2,0) = 0.;
};


/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] meshType element type to be created.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
TPZGeoMesh*
CreateGeoMesh(const MMeshType meshType, const TPZVec<int> &nDivs, const int dim,
              const int volId, const int bcId, MShapeType fShapeType);

/*
    Test the dimension of TPZHCurl approximation spaces with no high order
    gradient fields.
*/
template<class TSHAPE>
void TestHCurlNoGradsDim(const int &pOrder);

/*
    Test KernelHdiv problem
*/
template<class tshape>
void TestKernelHDiv(const int &xdiv, const int &pOrder, MShapeType fShapeType);

//Rotates a geometric by a given angle aroung a given axis
static void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);


TEST_CASE("HCurl no grads dimension", "[hdivkernel_mesh_tests]") {
    std::cout << "Testing dimension of Hcurl with no high order grads\n";
    const int pOrder = GENERATE(1,2,3);
  //SVD requires LAPACK
#ifndef PZ_USING_LAPACK
  return;
#endif
  TestHCurlNoGradsDim<pzshape::TPZShapeTriang>(pOrder);
  TestHCurlNoGradsDim<pzshape::TPZShapeQuad>(pOrder);
  TestHCurlNoGradsDim<pzshape::TPZShapeTetra>(pOrder);
  TestHCurlNoGradsDim<pzshape::TPZShapeCube>( pOrder);
  // TestHCurlNoGradsDim<pzshape::TPZShapePrism>(pOrder);
  std::cout << "Finish test dimension of HCurlNoGrads \n";
}

// Test 2D Kernel Hdiv: FOR DEBUGGING PURPOSES
TEST_CASE("HDiv Kernel", "[hdivkernel_mesh_tests]") {
  std::cout << "Testing HDiv Kernel \n";

  const int xdiv = GENERATE(2);
  const int pOrder = GENERATE(1,2);
  MShapeType shape = GENERATE(EHDivKernel);

  TestKernelHDiv<pzshape::TPZShapeTriang>(xdiv,pOrder,shape);
  TestKernelHDiv<pzshape::TPZShapeQuad>(xdiv,pOrder,shape);
  TestKernelHDiv<pzshape::TPZShapeTetra>(xdiv,pOrder,shape);
  TestKernelHDiv<pzshape::TPZShapeCube>(xdiv, pOrder, shape);
//   TestKernelHDiv<pzshape::TPZShapePrism>(xdiv,pOrder,shape);
  std::cout << "Finish test HDiv Kernel \n";
}

TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[3]["Domain"] = 1;
        stringtoint[2]["Surfaces"] = 2;

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);
    }

    return gmesh;
}

//Create 
TPZGeoMesh*
CreateGeoMesh(const MMeshType meshType, const TPZVec<int> &nDivs, const int dim,
              const int volId, const int bcId, MShapeType fShapeType)
{
    
    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    

    if (dim == 2 && fShapeType == EHDivKernel){
        auto gel = gmesh->ElementVec()[0];
        TPZGeoElSide geosidePoint(gel,0);
        TPZGeoElBC gelbcPoint(geosidePoint, -2);
        gel->ResetReference();        
    }

    return gmesh;
    
}

void CreateOrientedBoundaryElements(TPZGeoMesh *fGeoMesh)
{   
    for(auto gel : fGeoMesh->ElementVec())
    {
        if (gel->Dimension() < fGeoMesh->Dimension()) continue;
        
        int nSides = gel->NSides();
        //For tetrahedra only, loop over the surface sides
        for (int side = 0; side < nSides; side++){
            if (gel->SideDimension(side) != gel->Dimension()-1) continue;

            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            //Neighbour material id
            auto Nmatid = neighbour.Element()->MaterialId();

            /*  If the boundary has BC, delete the neighbour GeoElement and  
                create another one from TPZGeoElBC with the same material id
            */
            if (Nmatid != 2) continue;
            fGeoMesh->DeleteElement(neighbour.Element(),neighbour.Element()->Index()); 
            TPZGeoElBC gelbcWrap(gelside, Nmatid);
           
        }
    }
}

//Creates the flux mesh
TPZCompMesh * CreateFluxCMesh(TPZGeoMesh *fGeoMesh, int fDimension, int fDefaultPOrder, MShapeType fShapeType)
{   
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);

    //Inserts Null Materials
    std::set<int> allMat={1,2,-2};

    for (std::set<int>::iterator it=allMat.begin(); it!=allMat.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(fDimension);
        mat->SetBigNumber(1.e10);
    } 

    // //Creates computational elements
    if (fShapeType == EHDivKernel){
        cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivKernel);
        cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fDimension);
    } else if (fShapeType == EHDivConstant) {
        cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
        cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fDimension);
    }

    for (auto gel:fGeoMesh->ElementVec())
    {
        if (gel->Dimension() == 0){
            new TPZCompElH1<TPZShapePoint>(*cmesh,gel,H1Family::EH1Standard);
            // CreatePointEl(gel,*cmesh,H1Family::EH1Standard);
        }
    }
    
    cmesh->AutoBuild();

    return cmesh;
}

TPZCompMesh * CreatePressureCMesh(TPZGeoMesh *fGeoMesh, int fDimension, int fDefaultPOrder, MShapeType fShapeType)
{
    std::set<int> matIdVec = {0};
        
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(1.e10);
    }

    if (fShapeType == EHDivConstant){
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(1);//Domain elements
        mat->SetDimension(fDimension);
        cmesh->InsertMaterialObject(mat);
        mat -> SetBigNumber(1.e10);

        cmesh->SetAllCreateFunctionsDiscontinuous();

        cmesh->SetDefaultOrder(0);
        cmesh->SetDimModel(fDimension);
        cmesh->AutoBuild();

        int ncon = cmesh->NConnects();
        for(int i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cmesh->ConnectVec()[i]; 
            newnod.SetLagrangeMultiplier(1);
        }

        int nel = cmesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            if(!celdisc) continue;
            celdisc->SetConstC(1.);
            celdisc->SetTrueUseQsiEta();
            // espera-se elemento de pressao apenas para o contorno
            auto aaa = celdisc->Reference()->Dimension();
            // if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            // {
            //     DebugStop();
            // }
        }
    } else {

        cmesh->InitializeBlock();

        if(fDefaultPOrder == 0){
            cmesh->SetDefaultOrder(0);
            cmesh->SetDimModel(fDimension-1);
            cmesh->SetAllCreateFunctionsDiscontinuous();
        } else {
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
            cmesh->SetDefaultOrder(fDefaultPOrder-1);
            cmesh->SetDimModel(fDimension-1);
        }
        
        cmesh->AutoBuild(matIdVec);
            
        for(auto &newnod : cmesh->ConnectVec())
        {
            newnod.SetLagrangeMultiplier(1);
        }
    }
    

    return cmesh;
}

TPZMultiphysicsCompMesh * CreateMultiphysicsCMesh(TPZGeoMesh *fGeoMesh, int fDimension, int fDefaultPOrder, 
                                                  TPZVec<TPZCompMesh *> &meshvector,
                                                  BCType BC, MShapeType fShapeType)
{
    // gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);
   
    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyFlow(1,fDimension);//Domain
    mat->SetConstantPermeability(1.);
    if (fShapeType == EHDivConstant) mat->SetForcingFunction(forcefunction,1);

    // mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
    // mat->NStateVariables(3);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);

    //Boundary Conditions
    TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, 2, BC, val1, val2);
    if (fDimension == 2 && fShapeType == EHDivKernel)BCond->SetForcingFunctionBC(exactSolKernel2D,3);
    if (fDimension == 3 && fShapeType == EHDivKernel)BCond->SetForcingFunctionBC(exactSolKernel3D,3);
    if (fDimension == 2 && fShapeType == EHDivConstant)BCond->SetForcingFunctionBC(exactSolConstant2D,3);
    if (fDimension == 3 && fShapeType == EHDivConstant)BCond->SetForcingFunctionBC(exactSolConstant2D,3);
    cmesh->InsertMaterialObject(BCond);
   
    auto *matL2 = new TPZL2ProjectionCS<>(-2,0,1);
    cmesh->InsertMaterialObject(matL2);


    TPZManVector<MSpaceConfig> active(2,ETestTrial);
    active[0]=ETestTrial;
    active[1]=ETestTrial;
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    if (fShapeType == EHDivConstant) {
        TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);
    }
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 


    return cmesh;
}

TPZGeoMesh *CreateGeoMeshTetra(const MMeshType meshType, const TPZVec<int> &nDivs,
              const int volId, const int bcId)
{
    constexpr int dim{3};
    const TPZManVector<REAL,3> minX = {0,0,0};
    const TPZManVector<REAL,3> maxX = {1,1,1};
    constexpr bool createBoundEls{true};
    int nMats = 2*dim+1;
    
    //all bcs share the same id
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    
    TPZGeoMesh *gmesh = 
    TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,volId,createBoundEls);

    for (auto gel : gmesh->ElementVec())
    {
        if (gel->Dimension() == 3) continue;

        gel->SetMaterialId(bcId);

    }
    

  return gmesh;
}


#include "TPZShapeHCurl.h"
#include "TPZShapeHCurlNoGrads.h"
template<class TSHAPE>
void TestHCurlNoGradsDim(const int &pOrder){

    constexpr int nNodes = TSHAPE::NCornerNodes;
    constexpr int nSides = TSHAPE::NSides;
    constexpr int nConnects = nSides - nNodes;
    constexpr int dim = TSHAPE::Dimension;
    constexpr int curlDim = dim*2 - 3;
    
    TPZVec<int> connectorders(nConnects, pOrder);

    //let us count the shape functions
    int nunfiltfuncs{0}, nfiltfuncs{0};
    for(int ic = 0 ; ic < nConnects; ic++){
        const int conorder = connectorders[ic];
        nunfiltfuncs +=
            TPZShapeHCurl<TSHAPE>::ComputeNConnectShapeF(ic,conorder);
        nfiltfuncs +=
            TPZShapeHCurlNoGrads<TSHAPE>::ComputeNConnectShapeF(ic,conorder);
    }
    
    TPZVec<int64_t> ids(nNodes,0);
    for(int i = 0; i < nNodes; i++){ids[i] = i;}
    
    
    TPZShapeData data_unfilt, data_filt;
    TPZShapeHCurl<TSHAPE>::Initialize(ids, connectorders, data_unfilt);
    TPZShapeHCurlNoGrads<TSHAPE>::Initialize(ids, connectorders, data_filt);

    const int pord_int = TPZShapeHCurlNoGrads<TSHAPE>::MaxOrder(pOrder);
    typename TSHAPE::IntruleType intrule(2*pord_int);
    const int npts = intrule.NPoints();


    TPZFMatrix<REAL> curlcurl_filt(nfiltfuncs,nfiltfuncs,0.);
    TPZFMatrix<REAL> curlcurl_unfilt(nunfiltfuncs,nunfiltfuncs,0.);
    
    REAL weight;
    TPZVec<REAL> pt(dim,0.);

    
    for(int ipt = 0; ipt < npts; ipt++){
        intrule.Point(ipt,pt,weight);
        TPZFMatrix<REAL> phi_unfilt(dim,nunfiltfuncs);
        TPZFMatrix<REAL> curlphi_unfilt(curlDim, nunfiltfuncs);
        TPZShapeHCurl<TSHAPE>::Shape(pt,data_unfilt,phi_unfilt,curlphi_unfilt);
        TPZFMatrix<REAL> phi_filt (dim,nfiltfuncs);
        TPZFMatrix<REAL> curlphi_filt(curlDim, nfiltfuncs);
        TPZShapeHCurlNoGrads<TSHAPE>::Shape(pt,data_filt,phi_filt,curlphi_filt);
        
        //calc matrix curlphi_i curlphi_j
        //check if rank filt = rank unfilt
        //compare rank filt and dim filt
        for (int ish = 0; ish < nunfiltfuncs; ish++){
            for (int jsh = 0; jsh < nunfiltfuncs; jsh++){
                for (int d = 0; d < curlDim; d++){
                    curlcurl_unfilt(ish,jsh) += curlphi_unfilt(d,ish)*curlphi_unfilt(d,jsh) * weight;
                }
            }
        }
        for (int ish = 0; ish < nfiltfuncs; ish++){
            for (int jsh = 0; jsh < nfiltfuncs; jsh++){
                for (int d = 0; d < curlDim; d++){
                    curlcurl_filt(ish,jsh) += curlphi_filt(d,ish)*curlphi_filt(d,jsh) * weight;
                }
            }
        }
                
        
    }


    auto CalcRank = [](TPZFMatrix<STATE> & mat, const STATE tol){

        TPZFMatrix<STATE> S;
        TPZFMatrix<STATE> Udummy, VTdummy;
        mat.SVD(Udummy, S, VTdummy, 'N', 'N');

        int rank = 0;
        const int dimMat = S.Rows();
        for (int i = 0; i < dimMat; i++) {
            rank += S.GetVal(i, 0) > tol ? 1 : 0;
        }
        return std::make_tuple(S,rank);
    };
    

    std::ofstream out("CurlCheck.nb");
    curlcurl_unfilt.Print("Unf = ",out,EMathematicaInput);
    curlcurl_filt.Print("Fil = ",out,EMathematicaInput);
    static constexpr auto tol = std::numeric_limits<STATE>::epsilon()*10000;

    int rank_filt, rank_unfilt;
    TPZFMatrix<STATE> S_filt, S_unfilt;
    std::tie(S_filt, rank_filt) = CalcRank(curlcurl_filt,tol);
    std::tie(S_unfilt, rank_unfilt) = CalcRank(curlcurl_unfilt,tol);
    
    CAPTURE(pOrder);
    // CAPTURE(S_filt, S_unfilt);
    REQUIRE(rank_filt==rank_unfilt);
    REQUIRE(nfiltfuncs == rank_filt+TSHAPE::NCornerNodes-1);
}


template<class tshape>
void TestKernelHDiv(const int &xdiv, const int &pOrder, MShapeType shapeType){

    if (shapeType == EHDivConstant && tshape::Type() == EPrisma) return; // HDiv constant not implemented for prisms

    constexpr int volId{1};
    constexpr int bcId{2};
    MMeshType meshType;
    int DIM = tshape::Dimension;

    switch (tshape::Type())
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

    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    // ler uma malha do gmsh e ver se da certo. se der tem a ver com as orientações das faces.
    auto gmesh = CreateGeoMesh(meshType, nDivs, DIM, volId, bcId, shapeType);
    // TPZGeoMesh *gmesh = CreateGeoMeshTetra(meshType,nDivs,volId,bcId);
    // auto gmesh = ReadMeshFromGmsh("1tetra.msh");
    if (DIM == 3) CreateOrientedBoundaryElements(gmesh);

    //rotate the geoMesh
    int Axis;
    REAL theta, dump = 0.0;

    if (DIM == 3){
        theta = 48.0;
        Axis = 1;
        RotateGeomesh(gmesh, theta*dump, Axis);

        theta = -45.0;
        Axis = 2;
        RotateGeomesh(gmesh, theta*dump, Axis);
    }
    
    theta = 120.0;
    Axis = 3;
    RotateGeomesh(gmesh, theta*dump, Axis);

    TPZCompMesh * cmeshflux = CreateFluxCMesh(gmesh,DIM,pOrder,shapeType);
    TPZCompMesh * cmeshpres = CreatePressureCMesh(gmesh,DIM,pOrder,shapeType);

    TLaplaceExample1 LaplaceExact;
    LaplaceExact.fExact = TLaplaceExample1::ELaplace2D;


    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshflux;
    meshvector[1] = cmeshpres;
    // auto * cmesh = CreateMultiphysicsCMesh(gmesh,DIM,pOrder,meshvector,LaplaceExact.ExactSolution(),Dirichlet,shapeType);
    auto * cmesh = CreateMultiphysicsCMesh(gmesh,DIM,pOrder,meshvector,Dirichlet,shapeType);

    // Print mesh
    std::string txt =  "cmesh.txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);
#ifdef PZ_KERNELHDIV_DEBUG
    //Prints computational mesh properties
    std::string vtk_name = "cmesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);
#endif
    // Solve the problem
    TPZLinearAnalysis an(cmesh,RenumType::ENone);

    TPZSkylineStructMatrix<REAL> matskl(cmesh);

    //Filter equations
    if (DIM == 3 && shapeType == EHDivKernel){
        TPZHCurlEquationFilter<REAL> filter;

        TPZVec<int64_t> activeEqs;

        bool domainHybrid = false;
        if(filter.FilterEdgeEquations(cmesh, activeEqs, domainHybrid)){
            return;
        }
        const int neqs = activeEqs.size();
        
        matskl.EquationFilter().SetActiveEquations(activeEqs);
        std::cout << "Active equations = " << activeEqs.size() << std::endl;
    }

    an.SetStructuralMatrix(matskl);
    

    ///Setting a direct solver
    TPZStepSolver<REAL> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    // TPZStepSolver<STATE> jac;
    // REAL tol = 1.e-30;
    // jac.SetJacobi(100,tol,0);
    // jac.ShareMatrix(step);

    an.SetSolver(step);

    //Solve the problem
    an.Run();

    // ///solves the system
    // an.Solve();
    // auto mat = an.StructMatrix();

    ///Calculating approximation error  
    TPZManVector<REAL,5> error;

    // int64_t nelem = cmesh->NElements();
    // cmesh->LoadSolution(cmesh->Solution());
    // cmesh->ExpandSolution();
    // cmesh->ElementSolution().Redim(nelem, 5);
    
    // an.SetExact(LaplaceExact.ExactSolution(),pOrder);
    if (DIM == 2 && shapeType == EHDivKernel) an.SetExact(exactSolKernel2D,pOrder);
    if (DIM == 3 && shapeType == EHDivKernel) an.SetExact(exactSolKernel3D,pOrder);
    if (DIM == 2 && shapeType == EHDivConstant) an.SetExact(exactSolConstant2D,pOrder);
    if (DIM == 3 && shapeType == EHDivConstant) an.SetExact(exactSolConstant3D,pOrder);

    // an.PostProcessError(error);
    //Check element error
    int nel = cmesh->NElements();
    cmesh->LoadSolution(cmesh->Solution());
    for(int i = 0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension()<cmesh->Dimension()) continue;
        int ns = gel->NSides();
        int dim = gel->Dimension();
        TPZIntPoints *rule = gel->CreateSideIntegrationRule(ns-1, 4);//3
        int np = rule->NPoints();
        for(int ip=0; ip<np; ip++)
        {
            //TPZManVector<REAL,3> xi(2), xco(3), sol(1), exactsol(1);
            TPZManVector<REAL,3> xi(dim), xco(3);
            TPZManVector<STATE,3> sol(1), exactsol(1), solFlux(3);
            TPZFNMatrix<2,STATE> gradp(3,1);
            REAL weight;
            rule->Point(ip, xi, weight);
            gel->X(xi, xco);
            // cel->Solution(xi, 2, sol);
            cel->Solution(xi, 1, solFlux);
            
            if (DIM == 2 && shapeType == EHDivKernel) exactSolKernel2D(xco,exactsol,gradp);
            if (DIM == 3 && shapeType == EHDivKernel) exactSolKernel3D(xco,exactsol,gradp);
            if (DIM == 2 && shapeType == EHDivConstant) exactSolConstant2D(xco,exactsol,gradp);
            if (DIM == 3 && shapeType == EHDivConstant) exactSolConstant3D(xco,exactsol,gradp);

            STATE normError = std::sqrt((solFlux[0]+gradp(0,0))*(solFlux[0]+gradp(0,0)) + 
                                        (solFlux[1]+gradp(1,0))*(solFlux[1]+gradp(1,0)) + 
                                        (solFlux[2]+gradp(2,0))*(solFlux[2]+gradp(2,0))); // flux = -gradp

            if (normError > 1.e-6) {
                std::cout << "xi = " << xi <<std::endl;
                std::cout << "xco = " << xco <<std::endl;
                std::cout << "SOL = " << sol[0] << " " << exactsol[0] <<std::endl;
                std::cout << "flux = " << solFlux << "\n " << gradp <<std::endl;
            }

            // std::cout << "SOL = " << sol[0] << " " << exactsol[0] <<std::endl;
            // std::cout << "flux = " << solFlux << "\n " << gradp <<std::endl;

            REQUIRE(normError < 1.e-6); // Flux
            
        }
    }


#ifdef PZ_LOG
    // TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    // TPZManVector<std::string,10> scalnames(0), vecnames(2);
    // // scalnames[0] = "Pressure";
    // // scalnames[1] = "ExactPressure";
    // vecnames[0]= "Flux";
    // vecnames[1]= "ExactFlux";

    // int div = 0;
    // std::string plotfile = "solutionMDFB.vtk";
    // an.DefineGraphMesh(cmesh->Dimension(),scalnames,vecnames,plotfile);
    // an.PostProcess(div,cmesh->Dimension());
#endif
   
}


void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<REAL> RotationMatrix(3,3,0.0);

    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<REAL> iCoords(3,0.0);
    TPZVec<REAL> iCoordsRotated(3,0.0);
    
    //RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

    

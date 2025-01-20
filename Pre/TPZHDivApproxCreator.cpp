#include "TPZHDivApproxCreator.h"
#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "pzcmesh.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterialCS.h"

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "TPZBndCondT.h"
#include "TPZCompElDiscScaled.h"
#include "TPZCompMeshTools.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZCompElUnitaryLagrange.h"
#include "TPZCompElHDivDuplConnects.h"
#include "TPZCompElHDivDuplConnectsBound.h"

#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapetriang.h"
using namespace pzshape;

TPZHDivApproxCreator::TPZHDivApproxCreator(TPZGeoMesh *gmesh) : TPZApproxCreator(gmesh)
{ 
}

TPZHDivApproxCreator::~TPZHDivApproxCreator()
{
}

void TPZHDivApproxCreator::CheckSetupConsistency() {
    
    if (fProbType == ProblemType::EElastic) {
        std::cout << "WARNING! In elasticity problems, convergence has been checked only for Triangles and Tetrahedra. Please, check the convergence rates for the other topologies!\n";
    }

    if (fProbType == ProblemType::ENone) {
        std::cout << "You have to set a proper problem type!\n";
        DebugStop();
    }

    if (!fGeoMesh) {
        std::cout << "You have to set a GeoMesh!\n";
        DebugStop();
    }

    if (fMaterialVec.size() == 0){
        std::cout << "You have to set a Material!\n";
        DebugStop();
    }
 
    if (fIsRBSpaces && fHDivFam == HDivFamily::EHDivKernel){
        std::cout << "Are you sure about this?\n";
        DebugStop();
    }

    if (fHDivFam == HDivFamily::EHDivKernel){
        std::cout << "Pay attention to your material, the boundary conditions \n" <<
                     "need to be checked. Its implementation is not the same as \n" <<
                     "HDivStandard and HDivConstant spaces!!!! \n";
    }

    if (fExtraInternalPOrder < 0 || fExtraInternalPOrder > 2){
        std::cout << "Extra internal pOrder can only be 0, 1 or 2.\n";
        DebugStop();
    }
    
    if(fHybridType == HybridizationType::EStandardSquared){
        std::cout << "Standard Squared hybrization not implemented. Please implement if needed" << std::endl;
        DebugStop();
    }

    if(fHybridType == HybridizationType::ESemi && (fHDivFam != HDivFamily::EHDivConstant && fHDivFam != HDivFamily::EHDivOptimized)){
        std::cout << "The only HDiv spaces with available Semi hybridization is HDivConstant and HDivOptimized" << std::endl;
        DebugStop();
    }
  
    if (fProbType==ProblemType::EElastic && fHDivFam == HDivFamily::EHDivConstant && fShouldCondense && fHybridType != HybridizationType::ESemi) {
        std::cout << "CODE STOPPED! Condensing hdivconst spaces leads to singular K00. This happens because of the rotation space that has linear and constant functions."
        " One option is to separate the linear and constant functions is two different spaces. Please implement before using." << std::endl;
        DebugStop();
    }

}

void TPZHDivApproxCreator::CreateAtomicMeshes(TPZManVector<TPZCompMesh*,7>& meshvec, int& lagLevelCounter){
    if (fHybridType != HybridizationType::ENone){
        ComputePeriferalMaterialIds();
        AddHybridizationGeoElements();
        std::ofstream out("GeoMeshHybrid.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh, out);
    }

    bool isElastic = fProbType == ProblemType::EElastic;
    bool isDarcy = fProbType == ProblemType::EDarcy;
    if (isElastic){
        fNumMeshes = 3;
    } else if (isDarcy) {
        fNumMeshes = 2;
    } else {
        DebugStop();
    }
    if (fIsRBSpaces) fNumMeshes += 2;

    meshvec.Resize(fNumMeshes);
    int countMesh = 0;
    meshvec[countMesh++] = CreateHDivSpace();
    meshvec[countMesh-1]->SetName("HDiv");
    // lagLevelCounter = 1;
    meshvec[countMesh++] = CreateL2Space(fDefaultPOrder,lagLevelCounter++);
    meshvec[countMesh-1]->SetName("Pressure");

#ifdef PZDEBUG
    {
        std::ofstream out1("hdivmesh.txt");
        meshvec[0]->ComputeNodElCon();
        meshvec[0]->Print(out1);
        std::ofstream out("pressuremesh.txt");
        meshvec[1]->Print(out);
    }
#endif
    if (isElastic){
        meshvec[countMesh++] = CreateRotationSpace(fDefaultPOrder+fExtraInternalPOrder,lagLevelCounter++);
        meshvec[countMesh-1]->SetName("Rotation");
#ifdef PZDEBUG
        {
             std::ofstream out("rotation.txt");
             meshvec[countMesh-1]->Print(out);
        }
#endif
    }
    if (fIsRBSpaces){
        meshvec[countMesh++] = CreateConstantSpace(lagLevelCounter++);
        meshvec[countMesh-1]->SetName("DistFlux");
        meshvec[countMesh++] = CreateConstantSpace(lagLevelCounter++);
        meshvec[countMesh-1]->SetName("AvPressure");
#ifdef PZDEBUG
        {
             std::ofstream out("avpressure.txt");
             meshvec[countMesh-2]->Print(out);
            std::ofstream out2("distflux.txt");
            meshvec[countMesh-3]->Print(out2);
        }
#endif

    }
    
//    if(fShouldCondense && fHybridType != HybridizationType::ENone){
//        ChangeLagLevel(meshvec[1],lagLevelCounter);
//    }
    
//    std::ofstream outp("pmeshAfterLagLevelChange.txt");
//    meshvec[1]->Print(outp);

    if (countMesh != fNumMeshes) DebugStop();
}

void TPZHDivApproxCreator::CreateMultiPhysicsMesh(TPZManVector<TPZCompMesh*,7>& meshvec, int& lagLevelCounter, TPZMultiphysicsCompMesh*& cmeshmulti){
    bool isElastic = fProbType == ProblemType::EElastic;
    bool isDarcy = fProbType == ProblemType::EDarcy;
    
    cmeshmulti = CreateMultiphysicsSpace(meshvec);

    // PrintMeshElementsConnectInfo(cmeshmulti);
    if (fShouldCondense){
        if (isElastic && !fIsRBSpaces && fHybridType == HybridizationType::ENone){
            // In this case, the third (corresponding to the laglevCounter - 1) is the rotation mesh, whose can be condensed.
            // Displacements should be kept uncondensed
            TPZCompMeshTools::CondenseElements(cmeshmulti,lagLevelCounter-2,false);
        } else {
            if(fHybridType == HybridizationType::ENone){
                if (!fCondensePressure){
                    TPZCompMeshTools::CondenseElements(cmeshmulti,lagLevelCounter-1,false);
                }else{
                    TPZCompMeshTools::CondenseElements(cmeshmulti,lagLevelCounter,false);
                }
            }
            else{
                GroupAndCondenseElements(cmeshmulti);
            }
        }
    }
    
//    for(int i = 0 ; i < meshvec.size() ; i++) {
//        std::string str = "mesh_" + std::to_string(i) + ".txt";
//        std::ofstream out(str);
//        meshvec[i]->Print(out);
//        PrintMeshElementsConnectInfo(meshvec[i]);
//    }
//    PrintMeshElementsConnectInfo(cmeshmulti);
}

TPZMultiphysicsCompMesh * TPZHDivApproxCreator::CreateApproximationSpace(){
    std::cout << "\n---------------- Creating Space -----------------" << std::endl;
    CheckSetupConsistency();
    SetMeshElementType();

    int lagLevelCounter = 1;
    TPZManVector<TPZCompMesh*,7> meshvec(fNumMeshes);
    CreateAtomicMeshes(meshvec,lagLevelCounter);
    
    TPZMultiphysicsCompMesh *cmeshmulti = nullptr;
    CreateMultiPhysicsMesh(meshvec,lagLevelCounter,cmeshmulti);
    
    std::cout << "\n---------------- Finished Creating Space -----------------" << std::endl;
    
    return cmeshmulti;
}

void TPZHDivApproxCreator::PrintAllMeshes(TPZMultiphysicsCompMesh* mpcmesh){
    TPZVec<TPZCompMesh*>& meshvec = mpcmesh->MeshVector();
    for(int i = 0 ; i < meshvec.size() ; i++) {
        std::string str = "atomic_mesh_" + std::to_string(i) + ".txt";
        std::ofstream out(str);
        meshvec[i]->Print(out);
    }
    std::string str = "multiphysics_mesh.txt";
    std::ofstream out(str);
    mpcmesh->Print(out);
}

TPZCompMesh * TPZHDivApproxCreator::CreateHDivSpace(){
    
    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    if (fProbType == ProblemType::EDarcy || fProbType == ProblemType::EElastic){
        cmesh->SetDefaultOrder(fDefaultPOrder);
    } else {
        DebugStop();
    }
    
    cmesh->SetDimModel(dim);

    int nstate = 0;
    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
        if (!bnd && mat->Dimension() == dim){
            nstate = mat->NStateVariables(); // here we assume that all materials have same nstatevars
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),mat->Dimension(),mat->NStateVariables());
            cmesh->InsertMaterialObject(nullmat);
        } else if (bnd && bnd->Material()->Dimension() == dim){
            // if (bnd->Dimension() != dim-1) DebugStop();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim-1,mat->NStateVariables());
            cmesh->InsertMaterialObject(nullmat);
        }
    }
    
    if(nstate < 1) DebugStop();
    
    if(fHybridType != HybridizationType::ENone){
        // Then we insert the hdiv wrap materials
        const int wrapmatid = fHybridizationData.fWrapMatId;
        TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(wrapmatid,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    
    //Creates computational elements       
    cmesh->ApproxSpace().SetHDivFamily(fHDivFam);
    
    if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::EStandardSquared) {
        CreateHybridHDivMesh(cmesh);
    }
    else if (fHybridType == HybridizationType::ESemi){
        cmesh->ApproxSpace().SetAllCreateFunctionsHDivDuplConnects(dim);
        cmesh->AutoBuild();
        ActivateDuplicatedConnects(cmesh);
        SemiHybridizeDuplConnects(cmesh);
    } else {
        cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
        cmesh->AutoBuild();
    }
    
    cmesh->SetDimModel(dim);
    cmesh->InitializeBlock();
    
    if (fExtraInternalPOrder > 0){
        ChangeInternalOrder(cmesh,fDefaultPOrder+fExtraInternalPOrder);
    }
    
#ifdef PZDEBUG
//    std::ofstream out("fluxmesh_at_CreateHDivSpace.txt");
//    cmesh->Print(out);
//    PrintMeshElementsConnectInfo(cmesh);
#endif

    return cmesh;
}

void TPZHDivApproxCreator::CreateHybridHDivMesh(TPZCompMesh* cmesh) {
    const int dim = cmesh->Dimension();
    // First only build the vol compels disconnected
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    std::set<int> volmatids = GetVolumetricMatIds();
    cmesh->AutoBuild(volmatids);
    
    // Now, just the bcs but connected to the volumetric elements (therefore we do LoadReferences())
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    std::set<int> bcmatids = GetBCMatIds();
    cmesh->LoadReferences(); // So the bcs can see the volumetric elements during their creation
    cmesh->AutoBuild(bcmatids);

    // Last, we create the wraps
    fGeoMesh->ResetReference();
    for(auto cel : cmesh->ElementVec()) {
        if(!cel) DebugStop();
        TPZGeoEl* gel = cel->Reference();
        if(gel->Dimension() != dim) continue;
        const int firstSide = gel->FirstSide(dim-1);
        for(int is = firstSide ; is < gel->NSides()-1 ; is++) {
            TPZGeoElSide gelside(gel,is);
            if(gelside.HasNeighbour(bcmatids)) continue;
            cel->LoadElementReference();
            TPZGeoElSide neig = gelside.Neighbour();
            if(neig.Element()->MaterialId() != fHybridizationData.fWrapMatId) DebugStop();
            auto *celwrap = cmesh->CreateCompEl(neig.Element());
            neig.Element()->ResetReference();
            gel->ResetReference();
        }
    }
    
    FixSideOrientHydridMesh(cmesh);
}

void TPZHDivApproxCreator::FixSideOrientHydridMesh(TPZCompMesh* cmesh) {
    const int dim = cmesh->Dimension();
    int64_t nel = cmesh->NElements();
    std::set<int> bcmatids = GetBCMatIds();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        const int matid = gel->MaterialId();
        if(bcmatids.find(matid) != bcmatids.end()) continue;
//        if(gel->Dimension() != dim) continue;
        int firstside = gel->FirstSide(dim-1);
        for(int is = firstside; is< gel->NSides()-1; is++)
        {
            TPZGeoElSide gelside(gel,is);
            auto neigh = gelside.HasNeighbour(fHybridizationData.fWrapMatId);
            if(neigh){
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                intel->SetSideOrient(is, 1.);
            }
        }
    }
}

void TPZHDivApproxCreator::PrintMeshElementsConnectInfo(TPZCompMesh* cmesh) {
    std::cout << "\n===> Printing elements connect information" << std::endl;
    for (int i = 0; i < cmesh->NElements(); i++) {
        TPZCompEl *cel = cmesh->Element(i);
        if (!cel) continue;
        if (!cel->Reference()) continue;
        cel->LoadElementReference();
        int matid = cel->Reference()->MaterialId();
        auto nconnects = cel->NConnects();
        std::cout << "Element = " << i << ", dim = " << cel->Dimension() << ", mat = " << cel->Reference()->MaterialId() << ", nstate = " << cel->Material()->NStateVariables() << ", nconnects = " << nconnects << ": ";
        for (int j = 0; j < nconnects; j++)
        {
            std::cout << cel->ConnectIndex(j) << ", ";
        }
        std::cout << std::endl;
    }
}

TPZCompMesh * TPZHDivApproxCreator::CreateConstantSpace(const int lagLevel) {
    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDimModel(dim);

    int nstate = 0;
    if (fHDivFam != HDivFamily::EHDivKernel){
        for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
            TPZMaterial* mat = matpair.second;
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
            if (!bnd){
                if (mat->Dimension() != dim) continue;
                nstate = mat->NStateVariables();
                TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
                cmesh->InsertMaterialObject(nullmat);
                if(fProbType == ProblemType::EElastic && fIsRBSpaces && lagLevel > 2){
                    if(dim == 2) nullmat->SetNStateVariables(3);
                    else if(dim == 3) nullmat->SetNStateVariables(6);
                }
            }
        }
    }
    
    if(nstate < 1) DebugStop();

    cmesh->SetDefaultOrder(0);
    cmesh->SetAllCreateFunctionsDiscontinuous();

    cmesh->AutoBuild();

    {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cmesh->Element(el));
            if(!disc) DebugStop();
            disc->SetFalseUseQsiEta();
        }
    }
    int64_t ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagLevel);
    }

    return cmesh;
}

TPZCompMesh * TPZHDivApproxCreator::CreateL2Space(int pOrder, const int lagLevel){

    if (fProbType == ProblemType::EElastic){//Se triangular.
        pOrder = fDefaultPOrder+fExtraInternalPOrder;
    }

    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDimModel(dim);

    int nstate = 0;
    if (fHDivFam != HDivFamily::EHDivKernel){
        for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
            TPZMaterial* mat = matpair.second;
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
            if (!bnd){
                if (mat->Dimension() != dim) continue;
                nstate = mat->NStateVariables();
                TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
                cmesh->InsertMaterialObject(nullmat);
                if(fProbType == ProblemType::EElastic && fIsRBSpaces && lagLevel > 2){
                    if(dim == 2) nullmat->SetNStateVariables(3);
                    else if(dim == 3) nullmat->SetNStateVariables(6);
                }
            } 
        }
    }
    
    if(nstate < 1) DebugStop();

    //Creates computational elements
    switch (fHDivFam)
    {
    case HDivFamily::EHDivStandard:
    case HDivFamily::EHDivOptimized:
        if (pOrder > 0){
            cmesh->SetDefaultOrder(pOrder);
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        } else {
            cmesh->SetDefaultOrder(0);
            cmesh->SetAllCreateFunctionsDiscontinuous();
        }
        break;

    case HDivFamily::EHDivConstant:
        cmesh->SetDefaultOrder(0);
        cmesh->SetAllCreateFunctionsDiscontinuous();
        break;

    case HDivFamily::EHDivKernel:
        
        break;

    default:
        DebugStop();
        break;
    }
        
    cmesh->AutoBuild();
    
    if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::ESemi){
        // std::ofstream out("pmesh1.txt");
        // cmesh->Print(out);
        // Then we insert the lagrange materials into this l2 space cmesh
        const int lagmatid = fHybridizationData.fLagrangeMatId;
        TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(lagmatid,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
        
        int lagpord = pOrder;
        if (fProbType == ProblemType::EElastic) lagpord -= fExtraInternalPOrder;
        if (fHybridType == HybridizationType::ESemi) lagpord = 0;
        if (lagpord > 0){
            cmesh->SetDefaultOrder(lagpord);
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        } else {
            cmesh->SetDefaultOrder(0);
            cmesh->SetAllCreateFunctionsDiscontinuous();
        }
        
        // Now we only autobuild the langrange materials
        std::set<int> matidtobuild;
        matidtobuild.insert(lagmatid);
        cmesh->SetDimModel(dim-1);
        cmesh->AutoBuild(matidtobuild);
        // std::ofstream out2("pmesh2.txt");
        // cmesh->Print(out2);
    }

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(lagLevel);
    }
    
    // PrintMeshElementsConnectInfo(cmesh);

    return cmesh;
}


TPZMultiphysicsCompMesh * TPZHDivApproxCreator::CreateMultiphysicsSpace(const TPZVec<TPZCompMesh *> &meshvec){
    
    int dim = fGeoMesh->Dimension();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(dim);
    
    int nstate = 0;
    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *> (mat);
        if (!bnd){
            nstate = mat->NStateVariables();
            cmesh->InsertMaterialObject(mat);
        } else {
            cmesh->InsertMaterialObject(bnd);
        }
    }
    
    if(nstate < 1) DebugStop();
    
    if(fHybridType != HybridizationType::ENone){
        auto * nullmatWrap = new TPZNullMaterialCS<>(fHybridizationData.fWrapMatId,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmatWrap);

        auto * nullmatLag = new TPZNullMaterialCS<>(fHybridizationData.fLagrangeMatId,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmatLag);
    }

    TPZManVector<int> active(meshvec.size(),1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
    cmesh->BuildMultiphysicsSpace(active, meshvec);
    
    if(fHybridType != HybridizationType::ENone){
        InsertInterfaceMaterialObjects(cmesh);
        AddInterfaceComputationalElements(cmesh);
    }
    
    return cmesh;
}


TPZCompMesh * TPZHDivApproxCreator::CreateRotationSpace(const int pOrder, const int lagLevel){
    const int dim = fGeoMesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(fGeoMesh);
    
    // NS (Jan 2023): In discussion with Phil, it was agreed that the rotation space will always
    // be the same order of the approximation no matter what is the HDiv family
//    if (fHDivFam == HDivFamily::EHDivConstant){
//        cmesh->SetDefaultOrder(0);
//    } else if (fHDivFam == HDivFamily::EHDivStandard){
//        cmesh->SetDefaultOrder(pOrder);
//    }else {
//        DebugStop();
//    }
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    // This set will hold all the domain materials to later create the TPZCompElDiscScaled
    std::set<int> materialids;
    
    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
        if (!bnd){
            if (mat->Dimension() != dim) continue;
            const int matiddomain = mat->Id();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(matiddomain,dim,mat->NStateVariables());
            if(dim == 3) nullmat->SetNStateVariables(3);
            else if(dim == 2) nullmat->SetNStateVariables(1);
            else DebugStop();
            cmesh->InsertMaterialObject(nullmat);
            materialids.insert(matiddomain);
        }
    }
    
    
    fGeoMesh->ResetReference();
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if (!gel)continue;
        if(gel->HasSubElement()) continue;
        int matid = gel->MaterialId();
        if (materialids.find(matid) == materialids.end()) {
            continue;
        }
        new TPZCompElDiscScaled(*cmesh, gel);
        gel->ResetReference();
    }

    //cmesh->LoadReferences();
    //    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->AutoBuild();


    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagLevel);
    }

    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElDiscScaled *disc = dynamic_cast<TPZCompElDiscScaled *> (cel);
        if (!disc) {
            continue;
        }
        const REAL elementdim = cel->Reference()->ElementRadius();
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        disc->SetConstC(elementdim);
        disc->SetScale(1./elementdim);
    }
    cmesh->InitializeBlock();
    
//    PrintMeshElementsConnectInfo(cmesh);
    
    return cmesh;
}

void TPZHDivApproxCreator::GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) {
    if(mcmesh->GetSolType()==ESolType::EComplex){
        //we dont create condensed complex elements yet
        DebugStop();
    }
    const int dim = mcmesh->Dimension();
    
    TPZVec<int64_t> elementgroup; // it is resized inside AssociateElements()
    if (fHybridType == HybridizationType::ESemi){
        AssociateElementsDuplConnects(mcmesh,elementgroup);
    } else {
        AssociateElements(mcmesh, elementgroup);
    }

    int64_t nel = elementgroup.size();

    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = elementgroup[el];
        if(groupnum == -1) continue;
        auto iter = groupmap.find(groupnum);
        if (groupmap.find(groupnum) == groupmap.end()) {
            int64_t index;
            TPZElementGroup *elgr = new TPZElementGroup(*mcmesh);
            groupmap[groupnum] = elgr;
            elgr->AddElement(mcmesh->Element(el));
        }
        else
        {
            iter->second->AddElement(mcmesh->Element(el));
        }
    }
    mcmesh->ComputeNodElCon();
//    if(fHybridType == HybridizationType::EStandard)
//    {
//        int lagCTEspace2 = 4;
//        int64_t nconnects = mcmesh->NConnects();
//        for (int64_t ic = 0; ic<nconnects; ic++) {
//            TPZConnect &c = mcmesh->ConnectVec()[ic];
//            if(c.LagrangeMultiplier() == lagCTEspace2) c.IncrementElConnected();
//        }
//    }
    nel = mcmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mcmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompElT<STATE>(elgr);
            cond->SetKeepMatrix(false);
        }
    }
    mcmesh->CleanUpUnconnectedNodes();
    
    
//    int64_t nel = mcmesh->NElements();
//    mcmesh->LoadReferences();
//
//    TPZCompEl *cel;
//    TPZGeoEl *gel;
//    TPZGeoElSide gelside;
//    TPZGeoElSide neigh;
//    TPZCompEl *celTarget;
//
//    std::set<int64_t> targetMatIds;
//    targetMatIds.insert(fHybridizationData.fWrapMatId);
//    targetMatIds.insert(fHybridizationData.fInterfaceMatId);
//    std::set<int> bcmatids = GetBCMatIds();
//    for(auto it : bcmatids){
//        targetMatIds.insert(it);
//    }
//
//    //    std::cout << "Groups of connects " << groupindex << std::endl;
//    for (int64_t el = 0; el<nel; el++) {
//        cel = mcmesh->Element(el);
//        if(!cel){
//            continue;
//        }
//        if (mcmesh->Element(el)->Dimension() != dim) {
//            continue;
//        }
//        TPZElementGroup *elgr = new TPZElementGroup(*mcmesh);
//
//        gel = cel->Reference();
//
//        int firstSide;
//        if (gel->Dimension() >= 0) {
//            firstSide = gel->FirstSide(gel->Dimension() - 1);
//        } else {
//            DebugStop();
//        }
//
//        elgr->AddElement(cel);
//        for (int iside = firstSide; iside < gel->NSides() - 1; iside++) {
//            gelside = TPZGeoElSide(gel, iside);
//            neigh = gelside;
//            if(gelside.HasNeighbour(bcmatids)){
//                neigh = gelside.Neighbour();
//                if(bcmatids.find(neigh.Element()->MaterialId()) != bcmatids.end()){
//                    celTarget = neigh.Element()->Reference();
//                    elgr->AddElement(celTarget);
//                }
//                else{
//                    DebugStop(); // can this happen?
//                }
//            } else{
//
//                // TODO: Jeferson Nathan this methodology wont work for meshes with space restriction
//                // Another ideia is to create a vector of ints with the size of nconnects of the mesh
//                // and for each position we put the group to which that connect belongs
//                // Then, we use this information to group the connects!
//
//
//                neigh = gelside.Neighbour();
//                std::vector<int> targetMatIds= {fHybridizationData.fWrapMatId,fHybridizationData.fInterfaceMatId};
//                for (int ineigh = 0; ineigh < 2; ineigh++) {
//                    if (neigh.Element()->MaterialId() != targetMatIds[ineigh]) {
//                        DebugStop();
//                    }
//                    celTarget = neigh.Element()->Reference();
//                    elgr->AddElement(celTarget);
//                    neigh = neigh.Neighbour();
//                }
//
//            }
//        }
//    }
//    mcmesh->ComputeNodElCon();
//    int neoNel = mcmesh->NElements();
//    for (int64_t el = nel; el < neoNel; el++) {
//        TPZCompEl *cel = mcmesh->Element(el);
//        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
//        if (!elgr) {
//            DebugStop();
//        }
//        TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
//        cond->SetKeepMatrix(false);
//    }
//
//    mcmesh->CleanUpUnconnectedNodes();
    
}

void TPZHDivApproxCreator::AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mphys) {
    int numEl = mphys->NElements();

    auto fGeoMesh = mphys->Reference();
    fGeoMesh->ResetReference();
    mphys->LoadReferences();

    //In semi hybridization, only one interface compEl is required, but two interface CompEls are created in the Standard Hybridization
    //This map avoids the creation of an interface element at the same position twice
    std::map<int64_t,bool> semiHybridInterface;

    for(int iel = 0; iel < numEl; iel++){
        TPZCompEl *cel = mphys->ElementVec()[iel];
        if(!cel || !cel->Reference()){
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(matid != fHybridizationData.fWrapMatId){
            continue;
        }
        TPZGeoElSide gelside(gel);
        TPZGeoElSide LargeNeigh = gelside.HasLowerLevelNeighbour(fHybridizationData.fWrapMatId);
        TPZGeoElSide neighLag = gelside.HasNeighbour(fHybridizationData.fLagrangeMatId);
        if(!neighLag){
            continue;
        }
        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide cLagrange = neighLag.Reference();

        TPZGeoElSide ginterface = gelside.Neighbour();
        if(ginterface.Element()->MaterialId() != fHybridizationData.fInterfaceMatId){
            DebugStop();
        }
        
#ifdef PZDEBUG
//        std::cout << "===> Connecting elements with matids: " << celside.Element()->Reference()->MaterialId() << "\t" << cneigh.Element()->Reference()->MaterialId() << std::endl;
//        std::cout << "===> And indexes: " << celside.Element()->Reference()->Index() << "\t" << cneigh.Element()->Reference()->Index() << std::endl;
#endif
        if(LargeNeigh) {
            // create another interface
            TPZGeoElSide ginterface = neighLag.Neighbour();
            TPZCompElSide cLarge = LargeNeigh.Reference();
#ifdef PZDEBUG
            if(ginterface.Element()->MaterialId() != fHybridizationData.fInterfaceMatId) DebugStop();
#endif
            if (fHybridType == HybridizationType::ESemi){
                auto celint = new TPZCompElUnitaryLagrange(*mphys,ginterface.Element(),celside,cLagrange,semiHybridInterface[cLagrange.Element()->Index()]);
                semiHybridInterface[cLagrange.Element()->Index()] = true;
                
                TPZGeoElSide volside = cLarge.Reference().operator--();
                if (volside.Element()->Dimension() != cLarge.Element()->Dimension()+1) DebugStop();
                TPZCompEl *celvol = volside.Element()->Reference();

                //Sets the side orient
                TPZMultiphysicsElement *celmulti = dynamic_cast<TPZMultiphysicsElement *> (celvol); 
                TPZInterpolationSpace *celv = dynamic_cast<TPZInterpolationSpace *> (celmulti->Element(0)); 
                celint->SetSideOrient(celv->GetSideOrient(volside.Side()));
            } else {
                // cLarge.Element()->Print();
                TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*mphys,ginterface.Element(),cLarge,cLagrange);
            }
        }

        if (fHybridType == HybridizationType::ESemi){
            auto celint = new TPZCompElUnitaryLagrange(*mphys,ginterface.Element(),celside,cLagrange,semiHybridInterface[cLagrange.Element()->Index()]);
            semiHybridInterface[cLagrange.Element()->Index()] = true;

            TPZGeoElSide volside = celside.Reference().operator--();
            if (volside.Element()->Dimension() != celside.Element()->Dimension()+1) DebugStop();
            TPZCompEl *celvol = volside.Element()->Reference();

            //Sets the side orient
            TPZMultiphysicsElement *celmulti = dynamic_cast<TPZMultiphysicsElement *> (celvol); 
            TPZInterpolationSpace *celv = dynamic_cast<TPZInterpolationSpace *> (celmulti->Element(0)); 
            celint->SetSideOrient(celv->GetSideOrient(volside.Side()));
        } else {
            TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*mphys,ginterface.Element(),celside,cLagrange);
        }
        
        
    }
} 

void TPZHDivApproxCreator::ChangeLagLevel(TPZCompMesh* cmesh, const int newLagLevel) {
    if(cmesh->ApproxSpace().Style() != TPZCreateApproximationSpace::MApproximationStyle::EDiscontinuous &&
       cmesh->ApproxSpace().Style() != TPZCreateApproximationSpace::MApproximationStyle::EContinuous){
        DebugStop(); // Lagrange elements can only exist in L2 meshes
    }

    // Loop over all compels and pick the ones that have matid of lagrange multiplier
    for(auto cel : cmesh->ElementVec()){
        if(!cel) continue;
        const int matid = cel->Reference()->MaterialId();
        if(matid != fHybridizationData.fLagrangeMatId) continue;
        // loop over connects changing their laglevels
        const int ncon = cel->NConnects();
        for(int ic = 0 ; ic < ncon ; ic++){
            cel->Connect(ic).SetLagrangeMultiplier(newLagLevel);
        }
    }
}



void TPZHDivApproxCreator::ActivateDuplicatedConnects(TPZCompMesh *cmesh){

    for (int64_t i = 0; i < cmesh->NElements(); i++){

        auto celType = cmesh->Element(i)->Type();
        switch (celType)
        {
        case EPoint:
            {
                TPZCompElHDivDuplConnectsBound<TPZShapePoint> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapePoint> *> (cmesh->Element(i)); 
                if (celb) celb->ActiveDuplConnects(fConnDuplicated);
            }
        case EOned:
            {
                TPZCompElHDivDuplConnects<TPZShapeQuad> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeQuad> *> (cmesh->Element(i)); 
                if (celd) celd->ActiveDuplConnects(fConnDuplicated);
                TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *> (cmesh->Element(i)); 
                if (celb) celb->ActiveDuplConnects(fConnDuplicated);
            }
            break;

        case EQuadrilateral:
            {
                TPZCompElHDivDuplConnects<TPZShapeQuad> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeQuad> *> (cmesh->Element(i)); 
                if (celd) celd->ActiveDuplConnects(fConnDuplicated);
                TPZCompElHDivDuplConnectsBound<TPZShapeQuad> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeQuad> *> (cmesh->Element(i)); 
                if (celb) celb->ActiveDuplConnects(fConnDuplicated);
            }
            break;

        case ETriangle:
            {
                TPZCompElHDivDuplConnects<TPZShapeTriang> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTriang> *> (cmesh->Element(i)); 
                if (celd) celd->ActiveDuplConnects(fConnDuplicated);
                TPZCompElHDivDuplConnectsBound<TPZShapeTriang> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeTriang> *> (cmesh->Element(i)); 
                if (celb) celb->ActiveDuplConnects(fConnDuplicated);
            }
            break;

        case ECube:
            {
                TPZCompElHDivDuplConnects<TPZShapeCube> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeCube> *> (cmesh->Element(i)); 
                if (celd) celd->ActiveDuplConnects(fConnDuplicated);
            }
            break;

        case ETetraedro:
            {
                TPZCompElHDivDuplConnects<TPZShapeTetra> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTetra> *> (cmesh->Element(i)); 
                if (celd) celd->ActiveDuplConnects(fConnDuplicated);
            }
            break;

        default:
            DebugStop();
            break;
        }
    }

    PartitionDependMatrix(cmesh);
}

void TPZHDivApproxCreator::PartitionDependMatrix(TPZCompMesh *cmesh){

    // Now that all edge connects have been duplicated, we need to expand the dependency matrix of a connect if it exists
    // For this purpose we can use the map fConnDuplicated which relates the original connect and its corresponding duplicated one
    // PS: This operation is performed outside the element duplication of connects because all connects have to be already duplicated 
    // to the dependency matrix expansion to be possible
    std::map<int64_t,bool> fDepMatTreated;
    for (int64_t iCon = 0; iCon < cmesh->NConnects(); iCon++)
    {
        TPZConnect &c = cmesh->ConnectVec()[iCon];
        if (fDepMatTreated[iCon]) continue;
        if(c.HasDependency()){
            // std::cout << "Need to partition the dependency matrix." << std::endl;
            auto *ptr = c.FirstDepend();
            auto dept = dynamic_cast<TPZConnect::TPZDepend<STATE>*>(ptr);
            if(!dept){DebugStop();}
            const auto &depmat = dept->GetDepMatrix();
            while(ptr) {
                int64_t cIndex_old1 = iCon;
                int64_t cIndex_old2 = ptr->fDepConnectIndex;
                int64_t cIndex_new1 = fConnDuplicated[cIndex_old1];
                int64_t cIndex_new2 = fConnDuplicated[cIndex_old2];

                int rows = depmat.Rows();
                int cols = depmat.Cols();

                TPZFMatrix<REAL> DepMat00(1,1,0.), DepMat01(1,cols-1,0.), DepMat10(rows-1,1,0.), DepMat11(rows-1,cols-1,0.);
                DepMat00(0,0) = depmat(0,0);
                for (int i = 1; i < rows; i++){
                    DepMat10(i-1,0) = depmat(i,0);
                    for (int j = 1; j < cols; j++){
                        if (i == 1) DepMat01(0,j-1) = depmat(0,j);
                        DepMat11(i-1,j-1) = depmat(i,j);
                    }                  
                }
                // std::cout << "K00 " << DepMat00 << std::endl;
                // std::cout << "K01 " << DepMat01 << std::endl;
                // std::cout << "K10 " << DepMat10 << std::endl;
                // std::cout << "K11 " << DepMat11 << std::endl;
            
                ptr = ptr->fNext;
                
                c.RemoveDepend(cIndex_old1,cIndex_old2);
                
                fDepMatTreated[cIndex_old1] = true;
                fDepMatTreated[cIndex_new1] = true;
                
                //Get the duplicated connect and set the dependency matrix for all
                TPZConnect &c2 = cmesh->ConnectVec()[cIndex_new1];
                c2.RemoveDepend();

                //Dependency 00 - old1 + old2
                TPZConnect::TPZDependBase *depend00 = c.AddDependency(cIndex_old1, cIndex_old2, DepMat00, 0, 0, 1, 1);

                //Dependency 01 - old1 + new2
                TPZConnect::TPZDependBase *depend01 = c.AddDependency(cIndex_old1, cIndex_new2, DepMat01, 0, 0, 1, cols-1);

                //Dependency 10 - new1 + old2
                TPZConnect::TPZDependBase *depend10 = c2.AddDependency(cIndex_new1, cIndex_old2, DepMat10, 0, 0, rows-1, 1);

                //Dependency 11 - new1 + new2
                TPZConnect::TPZDependBase *depend11 = c2.AddDependency(cIndex_new1, cIndex_new2, DepMat11, 0, 0, rows-1, cols-1);

            }
        }
    }

    cmesh->InitializeBlock(); 
}

void TPZHDivApproxCreator::DisableDuplicatedConnects(TPZCompMesh *cmesh){


    for (int64_t i = 0; i < cmesh->NElements(); i++){

        auto celType = cmesh->Element(i)->Type();
        switch (celType)
        {
        case EPoint:
            {
                TPZCompElHDivDuplConnectsBound<TPZShapePoint> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapePoint> *> (cmesh->Element(i)); 
                if (celb) celb->InactiveDuplConnects();
            }
        case EOned:
            {
                TPZCompElHDivDuplConnects<TPZShapeQuad> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeQuad> *> (cmesh->Element(i)); 
                if (celd) celd->InactiveDuplConnects();
                TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *> (cmesh->Element(i)); 
                if (celb) celb->InactiveDuplConnects();
            }
            break;

        case EQuadrilateral:
            {
                TPZCompElHDivDuplConnects<TPZShapeQuad> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeQuad> *> (cmesh->Element(i)); 
                if (celd) celd->InactiveDuplConnects();
                TPZCompElHDivDuplConnectsBound<TPZShapeQuad> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeQuad> *> (cmesh->Element(i)); 
                if (celb) celb->InactiveDuplConnects();
            }
            break;

        case ETriangle:
            {
                TPZCompElHDivDuplConnects<TPZShapeTriang> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTriang> *> (cmesh->Element(i)); 
                if (celd) celd->InactiveDuplConnects();
                TPZCompElHDivDuplConnectsBound<TPZShapeTriang> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeTriang> *> (cmesh->Element(i)); 
                if (celb) celb->InactiveDuplConnects();
            }
            break;

        case ECube:
            {
                TPZCompElHDivDuplConnects<TPZShapeCube> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeCube> *> (cmesh->Element(i)); 
                if (celd) celd->InactiveDuplConnects();
            }
            break;

        case ETetraedro:
            {
                TPZCompElHDivDuplConnects<TPZShapeTetra> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTetra> *> (cmesh->Element(i)); 
                if (celd) celd->InactiveDuplConnects();
            }
            break;

        default:
            DebugStop();
            break;
        }
    }

    GroupDependMatrix(cmesh);
}

void TPZHDivApproxCreator::GroupDependMatrix(TPZCompMesh *cmesh){

    // Now that we have disabled the duplicated connect, we also need to restore the corresponding dependency matrix
    // We again use the map fConnDuplicated which relates the original connect and its corresponding duplicated one
    std::map<int64_t,bool> fDepMatTreated;

    for (int64_t iCon = 0; iCon < cmesh->NConnects(); iCon++){
        TPZConnect &c = cmesh->ConnectVec()[iCon];

        if (fDepMatTreated[iCon]) continue;

        if(c.HasDependency()){
            // std::cout << "Need to partition the dependency matrix." << std::endl;
            auto *ptr = c.FirstDepend();
            
            while(ptr) {

                int64_t cIndex_old1 = iCon;
                int64_t cIndex_old2;
                if (fConnDuplicated.find(ptr->fDepConnectIndex) != fConnDuplicated.end()){
                    cIndex_old2 = fConnDuplicated[ptr->fDepConnectIndex];
                } else {
                    int64_t findVal = ptr->fDepConnectIndex;    
                    auto it = find_if(fConnDuplicated.begin(), fConnDuplicated.end(), [findVal](const std::map<int64_t,int64_t>::value_type & p) {
                        return p.second == findVal;
                    });
                    cIndex_old2 = it->first;
                }

                int64_t cIndex_new1 = fConnDuplicated[cIndex_old1];
                int64_t cIndex_new2 = fConnDuplicated[cIndex_old2];

                TPZConnect &c2 = cmesh->ConnectVec()[cIndex_new1];
                auto *ptr2 = c2.FirstDepend();

                auto dep00 = ptr->HasDepend(cIndex_old2);
                auto dep01 = ptr->HasDepend(cIndex_new2);

                auto dep10 = ptr2->HasDepend(cIndex_old2);
                auto dep11 = ptr2->HasDepend(cIndex_new2);

                if (!dep00 || !dep01 || !dep10 || !dep11) DebugStop();

                auto GetDepMat= [](TPZConnect::TPZDependBase *dep)
                -> TPZFMatrix<STATE> const &{
                    auto dept = dynamic_cast<TPZConnect::TPZDepend<STATE>*>(dep);
                    if(!dept){DebugStop();}
                    return dept->GetDepMatrix();
                };


                const auto &matdep00 = GetDepMat(dep00);
                const auto &matdep01 = GetDepMat(dep01);
                const auto &matdep10 = GetDepMat(dep10);
                const auto &matdep11 = GetDepMat(dep11);
                
                int rows = matdep00.Rows() + matdep11.Rows();
                int cols = matdep00.Cols() + matdep11.Cols();
                
                TPZFMatrix<REAL> DepMat(rows,cols,0.);

                DepMat(0,0) = matdep00(0,0);
                for (int i = 1; i < rows; i++){
                    DepMat(i,0) = matdep10(i-1,0);
                    for (int j = 1; j < cols; j++){
                        if (i == 1) DepMat(0,j) = matdep01(0,j-1);
                        DepMat(i,j) = matdep11(i-1,j-1);
                    }                  
                }
                // std::cout << "K00 " << matdep00 << std::endl;
                // std::cout << "K01 " << dep01->fDepMatrix << std::endl;
                // std::cout << "K10 " << dep10->fDepMatrix << std::endl;
                // std::cout << "K11 " << dep11->fDepMatrix << std::endl;
                // std::cout << "K " << DepMat << std::endl;
                
                c.RemoveDepend(cIndex_old1,cIndex_old2);
                c.RemoveDepend(cIndex_old1,cIndex_new2);
                c2.RemoveDepend(cIndex_new1,cIndex_old2);
                c2.RemoveDepend(cIndex_new1,cIndex_new2);

                //Only one dependency - old1 + old2
                TPZConnect::TPZDependBase *depend = c.AddDependency(cIndex_old1, cIndex_old2, DepMat, 0, 0, rows, cols);
                               
                fDepMatTreated[cIndex_old1] = true;
                fDepMatTreated[cIndex_new1] = true;
                
                ptr = ptr->fNext; 
            }
        }
    }

}

void TPZHDivApproxCreator::SemiHybridizeDuplConnects(TPZCompMesh *cmesh) {
    int dim = cmesh->Dimension();
    int nel = cmesh->NElements();
    std::set<int> matBCId = GetBCMatIds();

    //For each element interface, duplicate a connect to hybridize the constant flux based on the side orientation
    //Then sets the same connect indexes to the neighbour wrap element.S
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int64_t index;

        if (cel->Dimension() != dim) continue;
        TPZInterpolatedElement *eldc = dynamic_cast<TPZInterpolatedElement *> (cel);

        auto nsides = gel->NSides();
        int nfacets = cel->Reference()->NSides(cel->Dimension()-1);
        int nEdges = 0;
        if (cmesh->Dimension() == 3){
            nEdges = gel->NSides(1);
        }
        int ncorner = cel->Reference()->NCornerNodes();

        for (int iface = 0; iface < nfacets; iface++)
        {
            TPZGeoElSide gelside(gel,iface+ncorner+nEdges);            
            TPZGeoElSide neighbour = gelside.Neighbour();
            if (matBCId.find(neighbour.Element()->MaterialId())!=matBCId.end()) continue;//Ignore the elements with BC neighbours

            int sideOrient = eldc->GetSideOrient(iface+ncorner+nEdges);

            // std::cout << "Element = " << el << ",iface = " << iface << " Neigh matid - " << neighbour.Element()->MaterialId() << ", side orient - " << sideOrient << std::endl;
            if (sideOrient == -1){
                TPZConnect &c = cel->Connect(2*iface);

                // We do not need to explicitly hybridize if the element is a sub-element
                if (c.HasDependency()) continue;

                // If the side has hanging nodes, we do not need to explicitly hybridize the connect
                TPZStack<TPZCompElSide > compElConnected;
                gelside.ConnectedCompElementList(compElConnected,1,1);
                if (compElConnected.size() > 2) continue;

                auto pOrder = cmesh->GetDefaultOrder();
                int nshape = c.NShape();//It is updated in the next loop
                int nstate = c.NState();//It can possibly change
                int64_t newConnect = cmesh->AllocateNewConnect(nshape,nstate,pOrder);
                cel->SetConnectIndex(2*iface,newConnect);
               
                if (neighbour.Element()->MaterialId() == fHybridizationData.fWrapMatId){
                    neighbour.Element()->Reference()->SetConnectIndex(0,newConnect);
                } else {
                    DebugStop();
                }
            }
        }
    }

    cmesh->ExpandSolution();

}

void TPZHDivApproxCreator::AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup) {
    
    // First, associate all the connects of a certain element to its index
    int64_t nel = cmesh->NElements();
    elementgroup.Resize(nel, -1);
    elementgroup.Fill(-1);
    int64_t nconnects = cmesh->NConnects();
    TPZVec<int64_t> groupindex(nconnects, -1);
    int dim = cmesh->Dimension();
    //The group index of connects belonging to volumetric elements equals its index.
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        elementgroup[cel->Index()] = cel->Index();
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto cindex : connectlist) {
#ifdef PZDEBUG
            if (groupindex[cindex] != -1) {
                DebugStop();
            }
#endif
            groupindex[cindex] = cel->Index();
        }
    }

    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int64_t celindex = cel->Index();
        
        TPZVec<int> connectgroup(connectlist.size());
        for(int i=0; i<connectlist.size(); i++) connectgroup[i] = groupindex[connectlist[i]];
        int64_t groupfound = -1;
        for (auto cindex : connectlist) {
            if (groupindex[cindex] != -1) {
                elementgroup[celindex] = groupindex[cindex];
                //two connects in the same element can't belong to different computational element groups,
                //but in interface elements, some connects might belong to a certain element group, while others might not be initialized.
                if(groupfound != -1 && groupfound != groupindex[cindex])
                {
                    DebugStop();
                }
                groupfound = groupindex[cindex];
            }
        }
    }
}

void TPZHDivApproxCreator::AssociateElementsDuplConnects(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup)
{
    
    int64_t nel = cmesh->NElements();
    elementgroup.Resize(nel, -1);
    elementgroup.Fill(-1);
    std::set<int> matBCId = GetBCMatIds();

    int64_t nconnects = cmesh->NConnects();
    TPZVec<int64_t> groupindex(nconnects, -1);
    int dim = cmesh->Dimension();
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        cel->LoadElementReference();
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        elementgroup[cel->Index()] = cel->Index();
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto cindex : connectlist) {
#ifdef PZDEBUG
            if (groupindex[cindex] != -1) {
                continue;
                // DebugStop();
            }
#endif
            groupindex[cindex] = cel->Index();
        }
    }

    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int64_t celindex = cel->Index();
        
        TPZVec<int> connectgroup(connectlist.size());
        for(int i=0; i<connectlist.size(); i++) connectgroup[i] = groupindex[connectlist[i]];
        int64_t groupfound = -1;
        //Already sets the element group to the volumetric elemet
        if (cel->Dimension() == dim) elementgroup[celindex] = celindex; 
        for (auto cindex : connectlist) {
            //Gets only the first connect to avoid trouble with the duplicated connect
            if (groupindex[cindex] != -1 && elementgroup[celindex] == -1) {
                elementgroup[celindex] = groupindex[cindex];
                //two connects in the same element can't belong to different computational element groups,
                //but in interface elements, some connects might belong to a certain element group, while others might not be initialized.
                if(groupfound != -1 && groupfound != groupindex[cindex])
                {
                    // DebugStop();
                }
                groupfound = groupindex[cindex];
            }
        }
    }
}

//
// Created by victor on 21/09/2022.
//

#include "TPZH1ApproxCreator.h"
#include "pzcmesh.h"
#include "TPZMaterial.h"
#include "TPZNullMaterial.h"
#include "pzintel.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZCompMeshTools.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("H1ApproxCreator"));
#endif

TPZH1ApproxCreator::TPZH1ApproxCreator(TPZGeoMesh * gmesh): TPZApproxCreator(gmesh){

}

void TPZH1ApproxCreator::CheckSetupConsistency(){

    if (fProbType == ProblemType::ENone) {
        std::cout << "You have to set a proper problem type!\n";
        DebugStop();
    }
    else if (fProbType != ProblemType::EDarcy || fHybridType == HybridizationType::ESemi){
        std::cout <<"Case not yet implemented.\n";
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

    if (fHybridType != HybridizationType::ENone){
        if(fGeoMesh->Dimension() == 3 && fExtraInternalPOrder < 2){
            std::cout << "Problem not well-formed";
            DebugStop();
        }
        if(fGeoMesh->Dimension() == 3 && fExtraInternalPOrder < 3){
            std::cout << "Problem not well-formed";
            DebugStop();
        }
    }
    else if(fHybridType == HybridizationType::ENone){
        if(fHybridizationData.fHybridizeBCLevel > 0) {
            std::cout << "You can't hybridize the boundary of this space";
            DebugStop();
        } else if(fIsRBSpaces){
            std::cout << "Currently, there is no support for enhanced spaces on classic H1 approximations.\n\tSetting fIsEnhancedSpaces = 0.\n";
            fIsRBSpaces = 0;
        }
        else if(fExtraInternalPOrder){
            std::cout << "Defining extra polynomial orders for classic H1 spaces makes little sense.\n\tSetting fExtraInternalPOrder = 0.\n";
            fExtraInternalPOrder = 0;
        }
    }

    if(fHybridType == HybridizationType::EStandard && fHybridizationData.fHybridizeBCLevel > 1){
        std::cout << "You can't hybridize twice the boundary of this space";
        DebugStop();
    }
}

TPZMultiphysicsCompMesh * TPZH1ApproxCreator::CreateApproximationSpace(){

    CheckSetupConsistency();

    if(fHybridType == HybridizationType::ENone){
        std::cout << __PRETTY_FUNCTION__ << "Perhaps TPZH1ApproxCreator::CreateClassicH1ApproximationSpace is more suited for your needs\n.";
        DebugStop();
    }

    if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::EStandardSquared) {
        ComputePeriferalMaterialIds();
        AddHybridizationGeoElements();
    }


    bool isDarcy = fProbType == ProblemType::EDarcy;
    if (isDarcy) {
            fNumMeshes = 2;
    }
    else {
        DebugStop();
    }
    if (fIsRBSpaces) fNumMeshes += 2;

    TPZManVector<TPZCompMesh*,7> meshvec(fNumMeshes);
    int countMesh = 0;
    if(HybridType() != HybridizationType::ENone)
        meshvec[countMesh++] = CreateBoundaryHDivSpace();
    meshvec[countMesh++] = CreateL2Space();

    if (fIsRBSpaces){
        int lagMult1 = 2, lagMult2 = 4;
        meshvec[countMesh++] = CreateConstantSpace(lagMult1);
        meshvec[countMesh++] = CreateConstantSpace(lagMult2);
    }

    if (countMesh != fNumMeshes) DebugStop();

    TPZMultiphysicsCompMesh *cmeshmulti = CreateMultiphysicsSpace(meshvec);

    return cmeshmulti;
}

TPZCompMesh * TPZH1ApproxCreator::CreateClassicH1ApproximationSpace() {

    CheckSetupConsistency();

    if(fHybridType != HybridizationType::ENone){
        std::cout << __PRETTY_FUNCTION__ << "Perhaps TPZH1ApproxCreator::CreateApproximationSpace is more suited for your needs\n.";
        DebugStop();
    }

    bool isDarcy = fProbType == ProblemType::EDarcy;
    if (!isDarcy)
        DebugStop();


    TPZCompMesh *H1mesh = CreateL2Space();

    if(fShouldCondense) {
        //GroupAndCondenseElements(H1mesh);
        //TPZCompMeshTools::CondenseElements(H1mesh,0,false);
        TPZCompMeshTools::CreatedCondensedElements(H1mesh,0,false);

    }
    return H1mesh;
}

TPZMultiphysicsCompMesh * TPZH1ApproxCreator::CreateMultiphysicsSpace(TPZManVector<TPZCompMesh *> meshvec){

    int dim = fGeoMesh->Dimension();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(dim);

    for (std::pair<int, TPZMaterial*> mat: fMaterialVec)
    {
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *> (mat.second);
        if (!bnd){
            cmesh->InsertMaterialObject(mat.second);
        } else {
            cmesh->InsertMaterialObject(bnd);
        }
    }

    bool isHybrid = fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::EStandardSquared;
    if(isHybrid)
        InsertWrapAndLagrangeMaterialObjects(cmesh);

    TPZManVector<MSpaceConfig> active(fNumMeshes,ETestTrial);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
    cmesh->BuildMultiphysicsSpace(active, meshvec);

    if(isHybrid) {
        InsertInterfaceMaterialObjects(cmesh);
        AddInterfaceComputationalElements(cmesh);
    }
    if(fShouldCondense)
        GroupAndCondenseElements(cmesh);

    cmesh->InitializeBlock();
    cmesh->ComputeNodElCon();

    return cmesh;
}

TPZCompMesh* TPZH1ApproxCreator::CreateConstantSpace(const int &lagMult){
    TPZCompMesh *constMesh = new TPZCompMesh(fGeoMesh);
    constMesh->SetAllCreateFunctionsDiscontinuous();
    constMesh->SetDefaultOrder(0);

    int meshDim = fGeoMesh->Dimension();
    std::set<int> volumetricMatIds = GetVolumetricMatIds();
    for(auto it : volumetricMatIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim);
        constMesh->InsertMaterialObject(mat);
    }

    constMesh->AutoBuild(volumetricMatIds);

    int64_t nc = constMesh->NConnects();
    for(int icon = 0; icon < nc; icon++){
        constMesh->ConnectVec()[icon].SetLagrangeMultiplier(lagMult);
    }
    return constMesh;
}

TPZCompMesh *TPZH1ApproxCreator::CreateL2Space()
{
    TPZCompMesh *L2mesh = new TPZCompMesh(fGeoMesh);
    InsertL2MaterialObjects(L2mesh);
    L2mesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    if(fHybridType != HybridizationType::ENone){
        L2mesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    int L2Order = fDefaultPOrder +fExtraInternalPOrder;

    L2mesh->SetDefaultOrder(L2Order);

    std::set<int> volumeMatIds = GetVolumetricMatIds();
    std::set<int> bcMatIds = GetBCMatIds();
    L2mesh->AutoBuild(volumeMatIds);

    int lagMult1 = 1, lagMult2 =1;
    if(fShouldCondense && !fIsRBSpaces){
        lagMult2 = 3;
        if(fHybridType == HybridizationType::ENone) {
            lagMult1 = lagMult2 =  0;
        }
    }

    int64_t nelem = L2mesh->NElements();
    for (int64_t el = 0; el<nelem; el++) {
        TPZCompEl *cel = L2mesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int nconnects = cel->NConnects();
        cel->Connect(0).SetLagrangeMultiplier(lagMult2);
        for (int ic=1; ic<nconnects; ic++) {
            cel->Connect(ic).SetLagrangeMultiplier(lagMult1);
        }
    }
    #ifdef LOG4CXX
    if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Number of volumetric pressure elements created " << nelem << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
    #endif
    // se nao condensar tem que mudar o nivel de lagrange multiplier de um connect
    if(fHybridType == HybridizationType::EStandardSquared || fHybridizationData.fHybridizeBCLevel == 0)
    {
        std::set<int> matids;
        std::set<int> bcMatIds = GetBCMatIds();
        if(fHybridizationData.fHybridizeBCLevel == 2 || fHybridizationData.fHybridizeBCLevel == 0)
            matids = bcMatIds;
        if(fHybridizationData.fHybridizeBCLevel == 2 )
            matids.insert(fHybridizationData.fSecondLagrangeMatId);
        L2mesh->SetDefaultOrder(fDefaultPOrder);
        L2mesh->AutoBuild(matids);
    }

    int lagMult3 = 5;
    if (fHybridType == HybridizationType::ENone){
        return L2mesh;
    }
    int64_t nelem_big = L2mesh->NElements();
    for (int64_t el = nelem; el<nelem_big; el++) {
        TPZCompEl *cel = L2mesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int nconnects = cel->NConnects();
        for (int ic=0; ic<nconnects; ic++) {
            cel->Connect(ic).SetLagrangeMultiplier(lagMult3);
        }
    }
    #ifdef LOG4CXX
    if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Number of lower dimensional pressure elements created " << nelem_big-nelem << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
    #endif

    //Wrap creation. Its connects are contained in the volumetric element, and since the H1 space is broken,
    //it requires special attention.
    if (fHybridType == HybridizationType::EStandard ||  fHybridType == HybridizationType::EStandardSquared) {

#ifdef LOG4CXX
        std::map<int,int> numcreated;
#endif

        auto HasBCNeighbour = [](const TPZGeoElSide &gelside,const std::set<int> &matbc) {
            TPZGeoElSide neighbour(gelside.Neighbour());
            while (neighbour != gelside) {
                int matid = neighbour.Element()->MaterialId();
                if (matbc.find(matid) != matbc.end()) return neighbour;
                neighbour = neighbour.Neighbour();
            }
            return TPZGeoElSide();
        };
        // create elements connected to the pressure elements, the neighbour of each pressure
        // element should be either matwrap or boundary condition
        L2mesh->Reference()->ResetReference();
        L2mesh->ApproxSpace().CreateDisconnectedElements(false);
        int meshdim =  fGeoMesh->Dimension();
        int64_t nel = L2mesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = L2mesh->Element(el);
            if (!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            int porder = intel->GetPreferredOrder();
            TPZGeoEl *gel = cel->Reference();
            if (!gel) DebugStop();
            int matid = gel->MaterialId();
            // loop only over volumetric elements
            if (volumeMatIds.find(matid) == volumeMatIds.end()) continue;
            if (gel->Dimension() != meshdim) {
                DebugStop();
            }
            int nsides = gel->NSides();
            int ncorner = gel->NCornerNodes();
            for (int side = ncorner; side < nsides; side++) {
                if (gel->SideDimension(side) != meshdim - 1) continue;
                TPZGeoElSide gelside(gel, side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                TPZGeoElSide bcneighbour = HasBCNeighbour(gelside, bcMatIds);
                // the boundary condition element should be the first neighbour
                if (fHybridizationData.fHybridizeBCLevel == 0 && bcneighbour && bcneighbour != neighbour) DebugStop();
                if (!bcneighbour && neighbour.Element()->MaterialId() != fHybridizationData.fWrapMatId) DebugStop();
                L2mesh->SetDefaultOrder(porder);
                int dir = gel->NormalOrientation(side);
                int wrapmat = fHybridizationData.fWrapMatId;
                // if neighbour exists then create the wrap conditionally
                {
                    // load the element reference so that the created element will share the connects
                    cel->LoadElementReference();
                    int64_t index;
                    TPZCompEl *bc_cel = L2mesh->ApproxSpace().CreateCompEl(neighbour.Element(), *L2mesh);
                    // reset the references so that future elements will not share connects
#ifdef LOG4CXX
                    numcreated[neighbour.Element()->MaterialId()]++;
#endif
                    gel->ResetReference();
                    bc_cel->Reference()->ResetReference();
                }
            }
        }

#ifdef LOG4CXX
        if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << std::endl;
        for(auto it : numcreated) sout << "Material ID " << it.first << " number of elements created " << it.second << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    }
    return L2mesh;
}

void TPZH1ApproxCreator::AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mphys)
{
#ifdef LOG4CXX
    std::map<int,int> numcreated;
#endif
    TPZGeoMesh *gmesh = mphys->Reference();
    gmesh->ResetReference();
    mphys->LoadReferences();
    int64_t nel = mphys->NElements();
    std::set<int> bcMatIds = GetBCMatIds();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mphys->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
//        if(matid != fH1Hybrid.fMatWrapId.first && matid != fH1Hybrid.fMatWrapId.second) continue;
        if(matid == fHybridizationData.fWrapMatId)
        {
            TPZCompEl *fluxel = FindHDivNeighbouringElement(cel);
            TPZGeoEl *fluxgel = fluxel->Reference();
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighmat = neighbour.Element()->MaterialId();
            if(neighmat != fHybridizationData.fInterfaceMatId)
            {
                DebugStop();
            }
            // determine if the interface should be positive or negative...
            int interfacematid = neighmat;
            int64_t index;
            TPZCompElSide celwrap(cel,gel->NSides()-1);
            TPZGeoElSide fluxgelside(fluxgel);
            TPZCompElSide fluxside = fluxgelside.Reference();
//            std::cout << "Creating interface from wrap element " << gel->Index() << " using neighbour " << neighbour.Element()->Index() <<
//             " and flux element " << fluxgel->Index() << std::endl;
            if(neighbour.Element()->Reference()) DebugStop();
            //TODO: Should we support more general interface classes?
            new TPZMultiphysicsInterfaceElement(*mphys,neighbour.Element(),celwrap,fluxside);
#ifdef LOG4CXX
            numcreated[neighmat]++;
#endif
        }
        if(fHybridType == HybridizationType::EStandardSquared && matid == fHybridizationData.fLagrangeMatId)
        {
            TPZGeoElSide gelsideflux(gel);
            TPZGeoElSide neighbour = gelsideflux.Neighbour();
            // we only handle the flux element neighbour to second lagrange multiplier
            if(neighbour.Element()->MaterialId() != fHybridizationData.fSecondInterfaceMatId) continue;
            TPZGeoElSide anteriorSecondInterface = neighbour; // Second interface before second Lagrange Multiplier
            TPZGeoElSide secondLagrange = anteriorSecondInterface.Neighbour();
            if(secondLagrange.Element()->MaterialId() != fHybridizationData.fSecondLagrangeMatId)
            {
                // If there is no other interface neighour, this must be a boundary element
                int pressmatid = secondLagrange.Element()->MaterialId();
                if(bcMatIds.find(pressmatid) == bcMatIds.end()) DebugStop();
            }
            else {
                TPZGeoElSide posteriorSecondInterface = secondLagrange.Neighbour();
                if(posteriorSecondInterface.Element()->MaterialId() !=  fHybridizationData.fSecondInterfaceMatId) DebugStop();
                // now we have to find the second flux element
                TPZGeoElSide lagrangeCandidate = posteriorSecondInterface.HasNeighbour(fHybridizationData.fLagrangeMatId);
                // if the fluxelement found is the first flux element
                if(lagrangeCandidate == gelsideflux) {
                    // we have to find a larger (lower level) flux element
                    lagrangeCandidate = gelsideflux.HasLowerLevelNeighbour(fHybridizationData.fLagrangeMatId);
                    if(!lagrangeCandidate)
                    {
                        TPZManVector<REAL,3> x(3,0.);
                        gelsideflux.CenterX(x);
                        std::cout << "gelsideflux center " << x << std::endl;
                        lagrangeCandidate = gelsideflux.HasLowerLevelNeighbour(fHybridizationData.fLagrangeMatId);
                        DebugStop();
                    }
#ifdef PZDEBUG
                    if(lagrangeCandidate == gelsideflux)
                    {
                        DebugStop();
                    }
#endif
                }
                {
                    TPZCompElSide celflux = lagrangeCandidate.Reference();
                    TPZCompElSide pressure = secondLagrange.Reference();
                    if(!celflux || !pressure) DebugStop();
                    int64_t index;
                    if(posteriorSecondInterface.Element()->Reference()) DebugStop();
                    //TODO: Should we support more general interface classes?
                    new TPZMultiphysicsInterfaceElement(*mphys,posteriorSecondInterface.Element(),celflux,pressure);
#ifdef LOG4CXX
                    numcreated[posteriorSecondInterface.Element()->MaterialId()]++;
#endif
                }
            }
            {
                TPZCompElSide celflux = gelsideflux.Reference();
                TPZCompElSide pressure = secondLagrange.Reference();
                if(!celflux || !pressure) DebugStop();
                int64_t index;
#ifdef LOG4CXX
                numcreated[firstlagrange.Element()->MaterialId()]++;
#endif
                if(anteriorSecondInterface.Element()->Reference()) DebugStop();
                //TODO: Should we support more general interface classes?
                new TPZMultiphysicsInterfaceElement(*mphys,anteriorSecondInterface.Element(),celflux,pressure);
            }
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << "Number of computational interface elements created by material id\n";
        for(auto it : numcreated) sout << "Material id " << it.first << " number of elements created " << it.second << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

TPZCompEl *TPZH1ApproxCreator::FindHDivNeighbouringElement(TPZCompEl *wrapcel)
{
    TPZGeoEl *gel = wrapcel->Reference();
    std::set<int> bcMatIds = GetBCMatIds();
    int nsides = gel->NSides();
    TPZGeoElSide gelside(gel,nsides-1);
    TPZStack<TPZCompElSide> celstack;
    gelside.EqualLevelCompElementList(celstack, 0, 0);
    int nelstack = celstack.size();
    for (int st = 0; st<nelstack; st++) {
        TPZCompElSide celside = celstack[st];
        TPZCompEl *cel = celside.Element();
        TPZGeoEl *gelneigh = cel->Reference();
        int matid = gelneigh->MaterialId();
        if (matid == fHybridizationData.fLagrangeMatId) {
            return cel;
        }
        if (fHybridizationData.fHybridizeBCLevel == 1 && bcMatIds.find(matid) != bcMatIds.end()) {
            return cel;
        }
    }
    TPZCompElSide cellarge = gelside.LowerLevelCompElementList2(0);
    if(!cellarge){
        gelside.Print(std::cout);
        cellarge = gelside.LowerLevelCompElementList2(0);
        DebugStop();
    }
    {
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide gellarge = cellarge.Reference();
        if(gellarge.Element()->MaterialId() == fHybridizationData.fLagrangeMatId)
        {
            return cellarge.Element();
        }
        gellarge.EqualLevelCompElementList(celstack, 0, 0);
        int nelst = celstack.size();
        for (int ist = 0; ist < nelst; ist++) {
            TPZCompElSide celside = celstack[ist];
            TPZGeoElSide gelside = celside.Reference();
            TPZGeoEl *candidate = gelside.Element();
            if(candidate->MaterialId() == fHybridizationData.fLagrangeMatId)
            {
                return celside.Element();
            }
        }
    }
    return NULL;
}

TPZCompMesh *TPZH1ApproxCreator::CreateBoundaryHDivSpace()
{
    TPZCompMesh *hdivmesh = new TPZCompMesh(fGeoMesh);

    //Inserting HDiv material
    if (fHybridType!= HybridizationType::EStandard || fHybridType!= HybridizationType::EStandardSquared) {
        int matid = fHybridizationData.fLagrangeMatId;
        auto nullmat = new TPZNullMaterial(matid);
        nullmat->SetDimension(fGeoMesh->Dimension()-1);
        nullmat->SetNStateVariables(1);
        hdivmesh->InsertMaterialObject(nullmat);
    }
    else {
        DebugStop();
    }

    std::set<int> bcMatIds = GetBCMatIds();
    if(fHybridizationData.fHybridizeBCLevel == 1)
    {
        for (auto matid:bcMatIds) {
            auto nullmat = new TPZNullMaterial(matid);
            nullmat->SetDimension(fGeoMesh->Dimension()-1);
            nullmat->SetNStateVariables(1);
            hdivmesh->InsertMaterialObject(nullmat);
        }
    }

    int lagMult = 0;
    hdivmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fGeoMesh->Dimension());
    hdivmesh->ApproxSpace().CreateDisconnectedElements(true);
    hdivmesh->SetDefaultOrder(fDefaultPOrder);
    hdivmesh->AutoBuild();
    int64_t nconnects = hdivmesh->NConnects();
    for (int ic=0; ic<nconnects; ic++) {
        hdivmesh->ConnectVec()[ic].SetLagrangeMultiplier(lagMult);
    }
    return hdivmesh;
}

void TPZH1ApproxCreator::InsertL2MaterialObjects(TPZCompMesh * L2Mesh){

    int dim = fGeoMesh->Dimension();
    auto insertMat = [](int id,int dim, TPZCompMesh *cmesh){
        auto nullmat = new TPZNullMaterial(id);
        nullmat->SetDimension(dim);
        nullmat->SetNStateVariables(1);
        cmesh->InsertMaterialObject(nullmat);
    };
    if(fHybridType == HybridizationType::ENone){
        for (auto matpair:fMaterialVec) {
            auto *mat = matpair.second;
            L2Mesh->InsertMaterialObject(mat);
        }
    }
    else {
        std::set<int> volMatIds = GetVolumetricMatIds();

        for (auto matid:volMatIds) {
            insertMat(matid, dim, L2Mesh);
        }
        if ((fHybridizationData.fHybridizeBCLevel == 2) || (fHybridizationData.fHybridizeBCLevel == 0)) {
            std::set<int> bcMatIds = GetBCMatIds();
            for (auto matid:bcMatIds) {
                insertMat(matid, dim - 1, L2Mesh);
            }
        }
        if (fHybridType != HybridizationType::ENone) {
            insertMat(fHybridizationData.fWrapMatId, dim - 1, L2Mesh);
        }
        if (fHybridType == HybridizationType::EStandardSquared) {
            insertMat(fHybridizationData.fSecondLagrangeMatId, dim - 1, L2Mesh);
        }
    }
}

void TPZH1ApproxCreator::GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh){
    if(mcmesh->GetSolType()==ESolType::EComplex){
        //we dont create condensed complex elements yet
        DebugStop();
    }
    /// same procedure as hybridize hdiv
    int64_t nel = mcmesh->NElements();
    TPZVec<int64_t> groupnumber(nel,-1);
    /// compute a groupnumber associated with each element
    if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::EStandardSquared)
        AssociateElements(mcmesh, groupnumber);

    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = groupnumber[el];
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
    if(fHybridType == HybridizationType::EStandard)
    {
        int lagCTEspace2 = 4;
        int64_t nconnects = mcmesh->NConnects();
        for (int64_t ic = 0; ic<nconnects; ic++) {
            TPZConnect &c = mcmesh->ConnectVec()[ic];
            if(c.LagrangeMultiplier() == lagCTEspace2) c.IncrementElConnected();
        }
    }
    nel = mcmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mcmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompElT<STATE>(elgr);
            cond->SetKeepMatrix(false);
        }
    }
}

void TPZH1ApproxCreator::AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup)
{
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

    int numloops = 1;
    if( fHybridType == HybridizationType::EStandardSquared) numloops = 2;
    // this loop will associate a first layer of interface elements to the group
    // this loop will associate the wrap elements with the group
    // if HybridSquared the connects of interface elements with matid fLagrangeMatId will be added to the group
    //    in the second pass :
    //    incorporate the HDiv lagrange elements in the group
    //    incorporate boundary HDiv elements
    //    incorporate the interface elements to the pressure lagrange DOFs in the group
    for (int iloop = 0; iloop < numloops; iloop++) for (TPZCompEl *cel : cmesh->ElementVec())
        {
            if (!cel || !cel->Reference()) {
                continue;
            }
            TPZStack<int64_t> connectlist;
            cel->BuildConnectList(connectlist);
            int matid = cel->Reference()->MaterialId();
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
            //Before this step, the connects that interface elements share with wrap elements, belong to an element group,
            //While connects shared with lagrange elements, belongs to no group.
            if(fHybridType == HybridizationType::EStandardSquared && matid == fHybridizationData.fInterfaceMatId)
            {
                for(auto cindex : connectlist)
                {
                    groupindex[cindex] = groupfound;
                }
            }
        }
}

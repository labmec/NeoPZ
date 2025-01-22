//
//  TPZHybridizeHDiv.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 16/05/18.
//

#include "TPZHybridizeHDiv.h"
#include "pzconnect.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZMaterial.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzintel.h"
#include "pzgeoelbc.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzgeoelside.h"
#include "pztrnsform.h"
#include "tpzintpoints.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"

#include "TPZMultiphysicsCompMesh.h"

#include <algorithm>

using namespace std;

TPZHybridizeHDiv::TPZHybridizeHDiv(TPZVec<TPZCompMesh*>& meshvec_Hybrid) {
    ComputePeriferalMaterialIds(meshvec_Hybrid);
    ComputeNState(meshvec_Hybrid);
}

void TPZHybridizeHDiv::ComputeNState(TPZVec<TPZCompMesh*>& meshvec_Hybrid) {
    if (meshvec_Hybrid.size() > 1 && meshvec_Hybrid[1]->NMaterials() > 0) {
        fNState = meshvec_Hybrid[1]->MaterialVec().begin()->second->NStateVariables();
    }
}

void TPZHybridizeHDiv::ComputePeriferalMaterialIds(TPZVec<TPZCompMesh*>& meshvec_Hybrid) {
    int maxMatId = std::numeric_limits<int>::min();
    for (auto &mesh : meshvec_Hybrid) {
        
        if(!mesh) continue;
        
        for (auto &mat : mesh->MaterialVec()) {
            maxMatId = std::max(maxMatId, mat.first);
        }
    }
    if(meshvec_Hybrid.size())
    {
        
        
        TPZGeoMesh *gmesh = meshvec_Hybrid[1]->Reference();
        int64_t nel = gmesh->NElements();
        for(int64_t el=0; el<nel; el++)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel) maxMatId = std::max(maxMatId,gel->MaterialId());
        }
    }
    if (maxMatId == std::numeric_limits<int>::min()) {
        maxMatId = 0;
    }
    fHDivWrapMatid = maxMatId + 1;
    fLagrangeInterface = maxMatId + 2;
    fInterfaceMatid.first = maxMatId + 3;
    fInterfaceMatid.second = maxMatId + 4;
    fLagrangeInterfaceEnd = maxMatId + 5;
}

/// split the connect between two neighbouring elements

std::tuple<int64_t, int> TPZHybridizeHDiv::SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    if (fHDivWrapMatid == 0 || fLagrangeInterface == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }
    TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
    //TPZCompMesh *pressuremesh = meshvec[1];
    //TPZGeoMesh *gmesh = fluxmesh->Reference();
    TPZGeoElSide gleft(left.Reference());
    TPZGeoElSide gright(right.Reference());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());
    TPZInterpolatedElement *intelright = dynamic_cast<TPZInterpolatedElement *> (right.Element());
    intelleft->SetSideOrient(left.Side(), 1);
    intelright->SetSideOrient(right.Side(), 1);
    TPZStack<TPZCompElSide> equalright;
    TPZConnect &cleft = intelleft->SideConnect(0, left.Side());

    if (cleft.HasDependency()) {
        // check whether the wrap of the large element was already created
        gright.EqualLevelCompElementList(equalright,1,0);
        // only one wrap element should exist
#ifdef PZDEBUG
        if(equalright.size() > 1)
        {
            DebugStop();
        }
        if(equalright.size()==1)
        {
            TPZGeoEl *equalgel = equalright[0].Element()->Reference();
            if(equalgel->Dimension() != fluxmesh->Dimension()-1)
            {
                DebugStop();
            }
        }
#endif
        // reset the reference of the wrap element
        if(equalright.size()) equalright[0].Element()->Reference()->ResetReference();
        cleft.RemoveDepend();
    }
    else
    {
        left.SplitConnect(right);
    }
    int sideorder = cleft.Order();
    fluxmesh->SetDefaultOrder(sideorder);
    // create HDivBound on the sides of the elements
    TPZCompEl *wrap1, *wrap2;
    {
        intelright->Reference()->ResetReference();
        intelleft->LoadElementReference();
        intelleft->SetPreferredOrder(sideorder);
        TPZGeoElBC gbc(gleft, fHDivWrapMatid);
        wrap1 = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh);
        if(cleft.Order() != sideorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap1);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelleft->Reference()->ResetReference();
        wrap1->Reference()->ResetReference();
    }
    // if the wrap of the large element was not created...
    if(equalright.size() == 0)
    {
        intelleft->Reference()->ResetReference();
        intelright->LoadElementReference();
        TPZConnect &cright = intelright->SideConnect(0,right.Side());
        int rightprevorder = cright.Order();
        intelright->SetPreferredOrder(cright.Order());
        TPZGeoElBC gbc(gright, fHDivWrapMatid);
        wrap2 = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh);
        if(cright.Order() != rightprevorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap2);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelright->Reference()->ResetReference();
        wrap2->Reference()->ResetReference();
    }
    else
    {
        wrap2 = equalright[0].Element();
    }
    wrap1->LoadElementReference();
    wrap2->LoadElementReference();
    int64_t pressureindex;
    int lagmatid = fLagrangeInterface;
    int pressureorder;
    {
        TPZGeoElBC gbc(gleft, lagmatid);
        pressureindex = gbc.CreatedElement()->Index();
        pressureorder = sideorder;
    }
    intelleft->LoadElementReference();
    intelright->LoadElementReference();
    return std::make_tuple(pressureindex, pressureorder);
}

/// find an element which shares the connect and has the same dimension
// else return 0

TPZCompElSide TPZHybridizeHDiv::RightElement(TPZInterpolatedElement *intel, int side) {
    bool isrestrained = false;
    {
        TPZConnect &c = intel->SideConnect(0, side);
        if (c.HasDependency()) {
            isrestrained = true;
        }
    }
    TPZGeoEl *gel = intel->Reference();
    TPZGeoElSide gelside(gel, side);

    if (isrestrained == true) {
        /// if the side is restrained we will hybridize between the element and the larger element
        TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
        if (!celside) DebugStop();
        TPZGeoEl *neigh = celside.Element()->Reference();
        /// we assume that a larger element should not be a boundary element
        if (neigh->Dimension() != gel->Dimension()) {
            DebugStop();
        }
        return celside;
    } else {
        // if the connect is not restrained
        // - there should be only one neighbour
        //   if there is more than one neighbour the element side has already been hybridized
        // - the neighbour should be of the same dimension
        //   if the neighbour is of lower dimension it is a boundary element
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        if (celstack.size() == 1) {
            TPZGeoEl *neigh = celstack[0].Element()->Reference();
            if (neigh->Dimension() == gel->Dimension()) {
                return celstack[0];
            }
        }
    }
    return TPZCompElSide();
}

std::tuple<int64_t, int> TPZHybridizeHDiv::SplitConnects(const TPZCompElSide &left, const TPZStack<TPZCompElSide> &cellsidestack, TPZVec<TPZCompMesh *> &meshvec_Hybrid, const bool isIntersectEnd) {
    if (fHDivWrapMatid == 0 || fLagrangeInterface == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }
    
    TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
    
    // All left connect needed variables
    TPZGeoElSide gleft(left.Reference());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());
    intelleft->SetSideOrient(left.Side(), 1);
    
    TPZConnect &cleft = intelleft->SideConnect(0, left.Side());
    
    TPZManVector<TPZInterpolatedElement*> intelvec(cellsidestack.size());
    int count = 0;
    if (cleft.HasDependency()){
        DebugStop(); // Please implement me!
    }
    else {
        for (auto& celside : cellsidestack){
            left.SplitConnect(celside);
            
            TPZInterpolatedElement *inteltochange = dynamic_cast<TPZInterpolatedElement *> (celside.Element());
            if (!inteltochange)
                DebugStop();
            inteltochange->SetSideOrient(celside.Side(), 1);
            intelvec[count++] = inteltochange;
        }
    }
    int sideorder = cleft.Order();
    fluxmesh->SetDefaultOrder(sideorder);
            
    // Create HdivBound first for left element
    TPZCompEl *wrapleft;
    TPZManVector<TPZCompEl*> wrapvec(cellsidestack.size());
    {
        for (TPZCompEl* cel : intelvec)
            cel->Reference()->ResetReference();
        intelleft->LoadElementReference();
        intelleft->SetPreferredOrder(sideorder);
        TPZGeoElBC gbc(gleft, fHDivWrapMatid); // creates geoelbc for side
        wrapleft = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh);
        if(cleft.Order() != sideorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrapleft);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelleft->Reference()->ResetReference();
        wrapleft->Reference()->ResetReference();
    }
    
    // Create HdivBound then for all others
    intelleft->Reference()->ResetReference();
    int i = 0;
    for (TPZCompElSide celside : cellsidestack) {
        TPZInterpolatedElement* intel = intelvec[i];
        intel->LoadElementReference();
        TPZConnect &con = intel->SideConnect(0, celside.Side());
        int prevorder = con.Order();
        intel->SetPreferredOrder(con.Order());
        
        TPZGeoElSide gelside(celside.Reference());
        TPZGeoElBC gbc(gelside,fHDivWrapMatid);
        TPZCompEl* wrap = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh);
        wrapvec[i++] = wrap;
        if (con.Order() != prevorder)
            DebugStop();
        
        TPZInterpolatedElement* intelwrap = dynamic_cast<TPZInterpolatedElement*>(wrap);
        int wrapside = gbc.CreatedElement()->NSides()-1;
        intelwrap->SetSideOrient(wrapside, 1);
        intel->Reference()->ResetReference();
        wrap->Reference()->ResetReference();
    }
    wrapleft->LoadElementReference();
    for (TPZCompEl* wrap : wrapvec)
        wrap->LoadElementReference();
    
    int64_t pressureindex;
    int pressureorder;
    int lagmatid = fLagrangeInterface;
    if (isIntersectEnd) {
        lagmatid = fLagrangeInterfaceEnd;
    }
    {
      TPZGeoElBC gbc(gleft, lagmatid);
      pressureindex = gbc.CreatedElement()->Index();
      pressureorder = sideorder;
    }
    
    intelleft->LoadElementReference();
    for (TPZInterpolatedElement* intel : intelvec)
        intel->LoadElementReference();
    
    return std::make_tuple(pressureindex,pressureorder);
}

bool TPZHybridizeHDiv::HybridizeInterface(TPZCompElSide& celsideleft, TPZInterpolatedElement *intelleft, int side, TPZVec<TPZCompMesh*>& meshvec_Hybrid,
                                          const bool isIntersectEnd) {
    
    // ==> Getting meshes
    TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    TPZGeoMesh *gmesh = fluxmesh->Reference();
    
    
    // ==> Splitting flux mesh connect
    // this is extremely expensive
    gmesh->ResetReference();
	if (fIdsToHybridize.size()){
		for (auto cel : fluxmesh->ElementVec()) {
			if(!cel) continue;
			const int celmatid = cel->Reference()->MaterialId();
			if (fIdsToHybridize.find(celmatid) != fIdsToHybridize.end()) {
				cel->LoadElementReference();
			}
		}
	}
    else {
        fluxmesh->LoadReferences();
    }    
    
    const bool isFractureIntersectionMesh = true;
    std::tuple<int64_t, int> pindexporder;
    if (isFractureIntersectionMesh) {
        TPZStack<TPZCompElSide> celsidestack;
        GetAllConnectedCompElSides(intelleft, side, celsidestack);
        pindexporder = SplitConnects(celsideleft, celsidestack, meshvec_Hybrid, isIntersectEnd);
    }
    else{
        TPZCompElSide celsideright = RightElement(intelleft, side);
        if (celsideright) {
            pindexporder = SplitConnects(celsideleft, celsideright, meshvec_Hybrid);
        }
        else {
    #ifdef PZDEBUG
            cout << "Cannot find right side connect. "
            "Interface could be already hybridized, skipping..." << endl;
    #endif
            return false;
        }
    }
    
    
    fluxmesh->InitializeBlock();
    fluxmesh->ComputeNodElCon();
    
    // ==> Creating lagrange multiplier element
    gmesh->ResetReference();
    pressuremesh->SetDimModel(gmesh->Dimension()-1);
    int64_t elindex;
    int order;
    std::tie(elindex, order) = pindexporder;
    TPZGeoEl *gel = gmesh->Element(elindex);
    TPZCompEl *cel = pressuremesh->ApproxSpace().CreateCompEl(gel, *pressuremesh);
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
    TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
    if (intel){
        intel->PRefine(order);
    } else if (intelDisc) {
        intelDisc->SetDegree(order);
        intelDisc->SetTrueUseQsiEta();
    } else {
        DebugStop();
    }
    int n_connects = cel->NConnects();
    for (int i = 0; i < n_connects; ++i) {
        cel->Connect(i).SetLagrangeMultiplier(2);
    }
    gel->ResetReference();
    
    pressuremesh->InitializeBlock();
    pressuremesh->SetDimModel(gmesh->Dimension());
    
    // Sets the fluxmesh as reference again if user wants to split another interface
    gmesh->ResetReference();
    fluxmesh->LoadReferences();
    
    return true;
}



/// split the connects between flux elements and create a dim-1 pressure element

void TPZHybridizeHDiv::HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    InsertPeriferalMaterialObjects(meshvec_Hybrid);
    
        TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
        TPZGeoMesh *gmesh = fluxmesh->Reference();
        int dim = gmesh->Dimension();
        gmesh->ResetReference();
        fluxmesh->LoadReferences();
        int64_t nel = fluxmesh->NElements();
        std::list<std::tuple<int64_t, int> > pressures;
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = fluxmesh->Element(el);
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
            if (!intel || intel->Reference()->Dimension() != dim) {
                continue;
            }
            // loop over the side of dimension dim-1
            TPZGeoEl *gel = intel->Reference();
            for (int side = gel->NCornerNodes(); side < gel->NSides() - 1; side++) {
                if (gel->SideDimension(side) != dim - 1) {
                    continue;
                }
                TPZCompElSide celside(intel, side);
                TPZCompElSide neighcomp = RightElement(intel, side);
                if (neighcomp) {
                    // SplitConnects returns the geometric element index and interpolation order
                    pressures.push_back(SplitConnects(celside, neighcomp, meshvec_Hybrid));
                }
            }
        }
        fluxmesh->InitializeBlock();
        fluxmesh->ComputeNodElCon();
        
    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    gmesh->ResetReference();
    pressuremesh->SetDimModel(gmesh->Dimension()-1);
    for (auto pindex : pressures) {
        int64_t elindex;
        int order;
        std::tie(elindex, order) = pindex;
        TPZGeoEl *gel = gmesh->Element(elindex);
        TPZCompEl *cel = pressuremesh->ApproxSpace().CreateCompEl(gel, *pressuremesh);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
        if (intel){
            intel->PRefine(order);
//            intel->SetSideOrder(gel->NSides() - 1, order);
        } else if (intelDisc) {
            intelDisc->SetDegree(order);
            intelDisc->SetTrueUseQsiEta();
        } else {
            DebugStop();
        }
        int n_connects = cel->NConnects();
        for (int i = 0; i < n_connects; ++i) {
            cel->Connect(i).SetLagrangeMultiplier(2);
        }
        gel->ResetReference();
    }
    pressuremesh->InitializeBlock();
    pressuremesh->SetDimModel(gmesh->Dimension());
}

void TPZHybridizeHDiv::CreateInterfaceElementsForGeoEl(TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid, TPZGeoEl *gel) {
    
    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    int dim = gel->Dimension()+1;
    
    cmesh_Hybrid->Reference()->ResetReference();
    cmesh_Hybrid->LoadReferences();
    
    TPZStack<TPZCompElSide> celstack;
//    TPZGeoEl *gel = cel->Reference();
    TPZCompEl *mphysics = gel->Reference();
    TPZGeoElSide gelside(gel, gel->NSides() - 1);
    TPZCompElSide celside = gelside.Reference();
    gelside.EqualLevelCompElementList(celstack, 0, 0);
    int count = 0;
    for (auto &celstackside : celstack) {
        if (celstackside.Reference().Element()->Dimension() == dim - 1) {
            int matid = fInterfaceMatid.first;
            if(count == 1) matid = fInterfaceMatid.second;
            TPZGeoElBC gbc(gelside, matid);
            // check if the right side has a dependency
            TPZCompEl *celneigh = celstackside.Element();
            if (celneigh->NConnects() != 1) {
                DebugStop();
            }
            
            TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh_Hybrid, gbc.CreatedElement(), celside, celstackside);
            count++;
        }
    }
    if (count == 1)
    {
        TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
        if(!clarge) DebugStop();
        TPZGeoElSide glarge = clarge.Reference();
        if (glarge.Element()->Dimension() == dim) {
            TPZGeoElSide neighbour = glarge.Neighbour();
            while(neighbour != glarge)
            {
                if (neighbour.Element()->Dimension() < dim) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if(neighbour == glarge) DebugStop();
            glarge = neighbour;
        }
        clarge = glarge.Reference();
        if(!clarge) DebugStop();
        TPZGeoElBC gbc(gelside, fInterfaceMatid.second);
        
        TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh_Hybrid, gbc.CreatedElement(), celside, clarge);
        count++;
    }
    
    pressuremesh->InitializeBlock();
}

void TPZHybridizeHDiv::CreateInterfaceElements(TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    if (fInterfaceMatid.first == 0 || fInterfaceMatid.second == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }
    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    int dim = pressuremesh->Dimension();
    TPZVec<TPZCompEl *> celpressure(pressuremesh->NElements(), 0);
    for (int64_t el = 0; el < pressuremesh->NElements(); el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if (cel && cel->Reference() && cel->Reference()->Dimension() == dim - 1) {
            celpressure[el] = cel;
        }
    }
    cmesh_Hybrid->Reference()->ResetReference();
    cmesh_Hybrid->LoadReferences();
    for (auto cel : celpressure) {
        if (!cel) continue;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoEl *gel = cel->Reference();
        TPZCompEl *mphysics = gel->Reference();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        int count = 0;
        for (auto &celstackside : celstack) {
            if (celstackside.Reference().Element()->Dimension() == dim - 1) {
                int matid = fInterfaceMatid.first;
                if(count == 1) matid = fInterfaceMatid.second;
                TPZGeoElBC gbc(gelside, matid);
                // check if the right side has a dependency
                TPZCompEl *celneigh = celstackside.Element();
                if (celneigh->NConnects() != 1) {
                    DebugStop();
                }
                
                TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh_Hybrid, gbc.CreatedElement(), celside, celstackside);
                count++;
            }
        }
        if (count == 1)
        {
            TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
            if(!clarge) DebugStop();
            TPZGeoElSide glarge = clarge.Reference();
            if (glarge.Element()->Dimension() == dim) {
                TPZGeoElSide neighbour = glarge.Neighbour();
                while(neighbour != glarge)
                {
                    if (neighbour.Element()->Dimension() < dim) {
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == glarge) DebugStop();
                glarge = neighbour;
            }
            clarge = glarge.Reference();
            if(!clarge) DebugStop();
            TPZGeoElBC gbc(gelside, fInterfaceMatid.second);
            
            TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh_Hybrid, gbc.CreatedElement(), celside, clarge);
            count++;
        }
        if (count != 2 && count != 0) {
            DebugStop();
        }
    }
    pressuremesh->InitializeBlock();
}

void TPZHybridizeHDiv::CreateInterfaceElements(TPZMultiphysicsCompMesh *cmesh_Hybrid)
{
    CreateInterfaceElements(cmesh_Hybrid, cmesh_Hybrid->MeshVector());
}

TPZCompMesh * TPZHybridizeHDiv::CreateMultiphysicsMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvector_Hybrid, double Lagrange_term_multiplier /* = 1. */) {
    TPZGeoMesh *gmesh = cmesh_HDiv->Reference();
    TPZCompMesh *cmesh_Hybrid = new TPZCompMesh(gmesh);
    cmesh_HDiv->CopyMaterials(*cmesh_Hybrid);
    InsertPeriferalMaterialObjects(cmesh_Hybrid, Lagrange_term_multiplier);
    cmesh_Hybrid->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh_Hybrid->AutoBuild();

    TPZBuildMultiphysicsMesh::AddElements(meshvector_Hybrid, cmesh_Hybrid);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector_Hybrid, cmesh_Hybrid);
    cmesh_Hybrid->LoadReferences();
    return cmesh_Hybrid;
}

TPZCompMesh * TPZHybridizeHDiv::CreateMultiphysicsMesh(TPZMultiphysicsCompMesh *cmesh_HDiv, double Lagrange_term_multiplier /* = 1. */) {
    TPZManVector<TPZCompMesh *, 3> meshvec_Hybrid(cmesh_HDiv->MeshVector().size(), 0);
    for (int i = 0; i < cmesh_HDiv->MeshVector().size(); i++) {
        meshvec_Hybrid[i] = cmesh_HDiv->MeshVector()[i]->Clone();
    }
    InsertPeriferalMaterialObjects(meshvec_Hybrid);
    HybridizeInternalSides(meshvec_Hybrid);

    TPZGeoMesh *gmesh = cmesh_HDiv->Reference();
    TPZMultiphysicsCompMesh *cmesh_Hybrid = new TPZMultiphysicsCompMesh(gmesh);
    cmesh_HDiv->CopyMaterials(*cmesh_Hybrid);
    InsertPeriferalMaterialObjects(cmesh_Hybrid, Lagrange_term_multiplier);

    TPZManVector<MSpaceConfig,5> active(meshvec_Hybrid.size(),ETestTrial);
    cmesh_Hybrid->BuildMultiphysicsSpace(active, meshvec_Hybrid);
    return cmesh_Hybrid;
}

/// create a multiphysics hybridized mesh based on and input mesh
void TPZHybridizeHDiv::ReCreateMultiphysicsMesh(TPZMultiphysicsCompMesh *cmesh_HDiv, double Lagrange_term_multiplier)
{
    TPZManVector<TPZCompMesh *, 5> meshvec_Hybrid = cmesh_HDiv->MeshVector();
    InsertPeriferalMaterialObjects(cmesh_HDiv, Lagrange_term_multiplier);
    InsertPeriferalMaterialObjects(meshvec_Hybrid);
    TPZManVector<MSpaceConfig> active = cmesh_HDiv->GetActiveApproximationSpaces();
    HybridizeInternalSides(meshvec_Hybrid);
    cmesh_HDiv->BuildMultiphysicsSpace(active, meshvec_Hybrid);

}

/// Associate elements with a volumetric element
// elementgroup[el] = index of the element with which the element should be grouped
// this method only gives effective result for hybridized hdiv meshes
void TPZHybridizeHDiv::AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup)
{
    int64_t nel = cmesh->NElements();
    elementgroup.Resize(nel, -1);
    elementgroup.Fill(-1);
    int64_t nconnects = cmesh->NConnects();
    TPZVec<int64_t> groupindex(nconnects, -1);
    int dim = cmesh->Dimension();
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
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
//        std::cout << "Analysing element " << cel->Index();
        int64_t groupfound = -1;
        for (auto cindex : connectlist) {
            if (groupindex[cindex] != -1) {
                // assign the element to the group
                elementgroup[cel->Index()] = groupindex[cindex];
                if(groupfound != -1 && groupfound != groupindex[cindex])
                {
                    DebugStop();
                }
//                if(groupfound == -1)
//                {
//                    std::cout << " added to " << groupindex[cindex];
//                }

                groupfound = groupindex[cindex];
            }
        }
//        std::cout << std::endl;
    }

}



/// group and condense the elements

void TPZHybridizeHDiv::GroupandCondenseElements(TPZCompMesh *cmesh) {

    const bool real_sol = cmesh->GetSolType()!=ESolType::EComplex;
    if(!real_sol){
        //we are not creating cplx condensed els yet
        DebugStop();
    }
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> groupnumber(nel,-1);
    /// compute a groupnumber associated with each element
    AssociateElements(cmesh, groupnumber);
    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = groupnumber[el];
        if(groupnum == -1) continue;
        auto iter = groupmap.find(groupnum);
        if (groupmap.find(groupnum) == groupmap.end()) {
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
            groupmap[groupnum] = elgr;
            elgr->AddElement(cmesh->Element(el));
        }
        else
        {
            iter->second->AddElement(cmesh->Element(el));
        }
//        std::cout << std::endl;
    }
    cmesh->ComputeNodElCon();
    nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompElT<STATE>(elgr);
            cond->SetKeepMatrix(false);
        }
    }
}

void TPZHybridizeHDiv::GroupandCondenseElements(TPZCompMesh *cmesh, int lagrange_keep) {

    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> groupnumber(nel,-1);
    /// compute a groupnumber associated with each element
    AssociateElements(cmesh, groupnumber);
    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = groupnumber[el];
        if(groupnum == -1) continue;
        auto iter = groupmap.find(groupnum);
        if (groupmap.find(groupnum) == groupmap.end()) {
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
            groupmap[groupnum] = elgr;
            elgr->AddElement(cmesh->Element(el));
        }
        else
        {
            iter->second->AddElement(cmesh->Element(el));
        }
//        std::cout << std::endl;
    }
    cmesh->ComputeNodElCon();
    nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            int nc = elgr->NConnects();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = elgr->Connect(ic);
                if(c.NElConnected() == 1 && c.LagrangeMultiplier() == lagrange_keep)
                {
                    c.IncrementElConnected();
                    break;
                }
            }
            TPZCondensedCompEl *cond = new TPZCondensedCompElT<STATE>(elgr);
            cond->SetKeepMatrix(false);
        }
    }
}


/// insert the material objects for HDivWrap and LagrangeInterface in the atomic meshes

void TPZHybridizeHDiv::InsertPeriferalMaterialObjects(TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    if (fLagrangeInterface == 0 || fHDivWrapMatid == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }


    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    int dim = pressuremesh->Dimension();
    
    if (!pressuremesh->FindMaterial(fLagrangeInterface)) {
        auto matPerif = new TPZNullMaterial(fLagrangeInterface);
        matPerif->SetDimension(dim-1);
        matPerif->SetNStateVariables(fNState);
        pressuremesh->InsertMaterialObject(matPerif);
    }
    if (!pressuremesh->FindMaterial(fLagrangeInterfaceEnd)) {
        auto matPerif = new TPZNullMaterial(fLagrangeInterfaceEnd);
        matPerif->SetDimension(dim-1);
        matPerif->SetNStateVariables(fNState);
        pressuremesh->InsertMaterialObject(matPerif);
    }
    
    if(meshvec_Hybrid[0]!=NULL){

        TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
        if (!fluxmesh->FindMaterial(fHDivWrapMatid)) {
            auto matPerif = new TPZNullMaterial(fHDivWrapMatid);
            matPerif->SetDimension(dim-1);
            matPerif->SetNStateVariables(fNState);
            fluxmesh->InsertMaterialObject(matPerif);
        }
    }
}

/// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the multiphysics mesh

void TPZHybridizeHDiv::InsertPeriferalMaterialObjects(TPZCompMesh *cmesh_Hybrid, double Lagrange_term_multiplier /* = 1. */) {
    if (fLagrangeInterface == 0 || fHDivWrapMatid == 0 || fInterfaceMatid.first == 0 ||
        fInterfaceMatid.second == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }
    TPZFNMatrix<1, STATE> xk(fNState, fNState, 0.), xb(fNState, fNState, 0.), xc(fNState, fNState, 0.), xf(fNState, 1, 0.);
    int dim = cmesh_Hybrid->Dimension();
    
    if (!cmesh_Hybrid->FindMaterial(fLagrangeInterface)) {
//        std::cout<<"LagrangeInterface MatId "<<fLagrangeInterface<<std::endl;
        auto matPerif = new TPZNullMaterialCS<STATE>(fLagrangeInterface);
        matPerif->SetNStateVariables(fNState);
        matPerif->SetDimension(dim-1);
        cmesh_Hybrid->InsertMaterialObject(matPerif);
    }
    if (!cmesh_Hybrid->FindMaterial(fLagrangeInterfaceEnd)) {
//        std::cout<<"LagrangeInterface MatId "<<fLagrangeInterfaceEnd<<std::endl;
        auto matPerif = new TPZNullMaterialCS<STATE>(fLagrangeInterfaceEnd);
        matPerif->SetNStateVariables(fNState);
        matPerif->SetDimension(dim-1);
        cmesh_Hybrid->InsertMaterialObject(matPerif);
    }
    if (!cmesh_Hybrid->FindMaterial(fHDivWrapMatid)) {
//        std::cout<<"HDivWrapMatid MatId "<<fHDivWrapMatid<<std::endl;
        auto matPerif = new TPZNullMaterialCS<STATE>(fHDivWrapMatid);
        matPerif->SetNStateVariables(fNState);
        matPerif->SetDimension(dim-1);
        cmesh_Hybrid->InsertMaterialObject(matPerif);
    }
    if (!cmesh_Hybrid->FindMaterial(fInterfaceMatid.first)) {
        
//        std::cout<<"InterfaceMatid MatId "<<fInterfaceMatid<<std::endl;
        auto *matleft = new TPZLagrangeMultiplierCS<STATE>(fInterfaceMatid.first, dim - 1, fNState);
        matleft->SetMultiplier(Lagrange_term_multiplier);
        cmesh_Hybrid->InsertMaterialObject(matleft);
    }
    if (!cmesh_Hybrid->FindMaterial(fInterfaceMatid.second)) {
        
//        std::cout<<"InterfaceMatid MatId "<<fInterfaceMatid<<std::endl;
        auto *matleft = new TPZLagrangeMultiplierCS<STATE>(fInterfaceMatid.second, dim - 1, fNState);
        matleft->SetMultiplier(Lagrange_term_multiplier);
        cmesh_Hybrid->InsertMaterialObject(matleft);
    }
}

std::tuple<TPZCompMesh*, TPZVec<TPZCompMesh*> > TPZHybridizeHDiv::Hybridize(TPZCompMesh* cmesh_HDiv, TPZVec<TPZCompMesh*>& meshvec_HDiv, bool group_elements, double Lagrange_term_multiplier /* = 1. */) {
    TPZManVector<TPZCompMesh *, 3> meshvec_Hybrid(meshvec_HDiv.size(), 0);
    for (int i = 0; i < meshvec_HDiv.size(); i++) {
        meshvec_Hybrid[i] = meshvec_HDiv[i]->Clone();
    }
    ComputePeriferalMaterialIds(meshvec_Hybrid);
    ComputeNState(meshvec_Hybrid);
    /// insert the material objects for HDivWrap and LagrangeInterface
    InsertPeriferalMaterialObjects(meshvec_Hybrid);
    HybridizeInternalSides(meshvec_Hybrid);
    TPZCompMesh *cmesh_Hybrid = CreateMultiphysicsMesh(cmesh_HDiv, meshvec_Hybrid, Lagrange_term_multiplier);
    CreateInterfaceElements(cmesh_Hybrid, meshvec_Hybrid);
    if (group_elements){
        GroupandCondenseElements(cmesh_Hybrid);
    }
    return std::make_tuple(cmesh_Hybrid, meshvec_Hybrid);
}

/// make a hybrid mesh from a H(div) multiphysics mesh
TPZMultiphysicsCompMesh *TPZHybridizeHDiv::Hybridize(TPZMultiphysicsCompMesh *multiphysics, bool group_elements, double Lagrange_term_multiplier)
{
    ComputePeriferalMaterialIds(multiphysics->MeshVector());
    ComputeNState(multiphysics->MeshVector());
//    InsertPeriferalMaterialObjects(multiphysics->MeshVector());
//    HybridizeInternalSides(multiphysics->MeshVector());
    TPZCompMesh *cmesh = CreateMultiphysicsMesh(multiphysics);
    TPZMultiphysicsCompMesh *result = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    CreateInterfaceElements(result);
    if(group_elements)
    {
      GroupandCondenseElements(result);
    }
    if(!result) DebugStop();
    return result;

}

/// make a hybrid mesh from a H(div) multiphysics mesh
void TPZHybridizeHDiv::HybridizeGivenMesh(TPZMultiphysicsCompMesh &multiphysics, bool group_elements, double Lagrange_term_multiplier)
{
    ComputePeriferalMaterialIds(multiphysics.MeshVector());
    ComputeNState(multiphysics.MeshVector());
    ReCreateMultiphysicsMesh(&multiphysics);
    CreateInterfaceElements(&multiphysics);
    if (group_elements) {
        GroupandCondenseElements(&multiphysics);
    }
}



static void CompareFluxes(TPZCompElSide &left, TPZCompElSide &right, std::ostream &out)
{
    TPZGeoElSide gleft = left.Reference();
    TPZGeoElSide gright = right.Reference();
    TPZTransform<REAL> trleft = gleft.Element()->SideToSideTransform(gleft.Side(),gleft.Element()->NSides()-1);
    TPZTransform<REAL> trright = gright.Element()->SideToSideTransform(gright.Side(),gright.Element()->NSides()-1);
    int sidedim = gleft.Dimension();
    int meshdim = gleft.Element()->Mesh()->Dimension();
    TPZTransform<REAL> tr(sidedim); 
    gleft.SideTransform3(gright,tr);
    TPZIntPoints *intpts = gleft.Element()->CreateSideIntegrationRule(gleft.Side(),5);
    TPZManVector<REAL,3> ptleft(sidedim,0.),ptright(sidedim,0.),ptleftvol(meshdim,0),ptrightvol(meshdim,0);
    REAL weight;
    int npoints = intpts->NPoints();
    for(int ip = 0; ip<npoints; ip++)
    {
        intpts->Point(ip,ptleft,weight);
        tr.Apply(ptleft,ptright);
        trleft.Apply(ptleft,ptleftvol);
        trright.Apply(ptright,ptrightvol);
        TPZManVector<STATE,3> fluxleft(3,0.),fluxright(3,0.);
        left.Element()->Solution(ptleftvol,0,fluxleft);
        right.Element()->Solution(ptrightvol,0,fluxright);
        TPZManVector<REAL,3> normal(3,0.);
        {
            TPZFNMatrix<9,REAL> axes(sidedim,3),jacobian(sidedim,sidedim),jacinv(sidedim,sidedim);
            REAL detjac;
            gleft.Jacobian(ptleft,jacobian,axes,detjac,jacinv);
            normal[0] = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
            normal[1] = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
            normal[2] = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
        }
        STATE fluxnormalleft = 0.,fluxnormalright = 0.;
        for(int i=0; i<3; i++)
        {
            fluxnormalleft += fluxleft[i]*normal[i];
            fluxnormalright += fluxright[i]*normal[i];
        }
        REAL diff = fabs(fluxnormalleft-fluxnormalright);
        if(diff > 1.e-6)
        {
            TPZManVector<REAL,3> leftx(3),rightx(3);
            gleft.X(ptleft,leftx);
            gright.X(ptright,rightx);
            out << "Left geo element index " << gleft.Element()->Index();
            out << " Right geo element index " << gright.Element()->Index();
            out << " xleft " << leftx << " xright " << rightx << std::endl;
            out << " fluxleft " << fluxleft << std::endl;
            out << " fluxright" << fluxright << std::endl;
            out << "normal " << normal << std::endl;
            out << "fluxnormalleft " << fluxnormalleft << std::endl;
            out << "fluxnormalright " << fluxnormalright << std::endl;
        }
    }
    delete intpts;
}

/// verify the consistency of the solution of the flux mesh
void TPZHybridizeHDiv::VerifySolutionConsistency(TPZCompMesh *fluxmesh, std::ostream &out)
{
    int64_t nel = fluxmesh->NElements();
    fluxmesh->Reference()->ResetReference();
    fluxmesh->LoadReferences();
    int meshdim = fluxmesh->Dimension();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != meshdim) continue;
        int nsides = gel->NSides();
        for(int side=0; side<nsides; side++)
        {
            if(gel->SideDimension(side) != meshdim-1) continue;
            TPZStack<TPZCompElSide> celsides;
            TPZGeoElSide gelside(gel,side);
            gelside.EqualLevelCompElementList(celsides,1,0);

            {
                TPZCompElSide celneigh = gelside.LowerLevelCompElementList2(1);
                if(celneigh) 
                {
                    celsides.push_back(celneigh);
                    celneigh.Reference().EqualLevelCompElementList(celsides,1,0);
                }
            }
            TPZCompElSide celneighselect;
            int nummatch = 0;
            for(int i=0; i<celsides.size(); i++)
            {
                if(celsides[i].Reference().Element()->Dimension() == meshdim) 
                {
                    celneighselect = celsides[i];
                    nummatch++;
                }
            }
            if(nummatch > 1) DebugStop();
            TPZCompElSide celside(gelside.Reference());
            if(celneighselect)
            {
                CompareFluxes(celside,celneighselect,out);
            }
        }
    }
}

void TPZHybridizeHDiv::GetAllConnectedCompElSides(TPZInterpolatedElement *intel, int side, TPZStack<TPZCompElSide> &celsidestack) {
    bool isrestrained = false;
    {
        TPZConnect &c = intel->SideConnect(0, side);
        if (c.HasDependency()) {
            isrestrained = true;
        }
    }
    TPZGeoEl *gel = intel->Reference();
    TPZGeoElSide gelside(gel, side);

    if (isrestrained == true) {
        /// if the side is restrained we will hybridize between the element and the larger element
        DebugStop();
    } else {
        // if the connect is not restrained
        // - the neighbour should be of the same dimension
        //   if the neighbour is of lower dimension it is a boundary element
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        for (auto& cel : celstack) {
            TPZGeoEl *neigh = cel.Element()->Reference();
            if (neigh->Dimension() == gel->Dimension()) {
                const int celmatid = cel.Element()->Reference()->MaterialId();
                if (fIdsToHybridize.size()) {
                    if (fIdsToHybridize.find(celmatid) != fIdsToHybridize.end()) {
                        celsidestack.push_back(cel);
                    }
                }
                else {
                    celsidestack.push_back(cel);
                }
                
            }
        } // cel
    } // else
}

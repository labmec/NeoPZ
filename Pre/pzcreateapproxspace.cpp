/**
 * @file
 * @brief Contains the implementation of the functions to creates computational elements.
 */

#include "pzcreateapproxspace.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzcondensedcompel.h"
#include "pzinterpolationspace.h" 
#include "TPZCompElH1.h"

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"

#include "pzgeopoint.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "TPZMaterial.h"
#include "TPZInterfaceEl.h"

#include "pzcompelwithmem.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzcreateapproximationspace");
#endif

/** @brief Creates computational point element */
TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);
/** @brief Creates computational linear element */
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);
/** @brief Creates computational quadrilateral element */
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);
/** @brief Creates computational triangular element */
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);
/** @brief Creates computational cube element */
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);
/** @brief Creates computational prismal element */
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);
/** @brief Creates computational pyramidal element */
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);
/** @brief Creates computational tetrahedral element */
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam);


using namespace pzshape;

TPZCompEl *CreateNoElement(TPZGeoEl *gel,TPZCompMesh &mesh) {
#ifdef PZ_LOG
    if (logger.isWarnEnabled()) {
        std::stringstream sout;
        sout << "Mesh dimension " << mesh.Dimension() << " gel dimension " << gel->Dimension() << " will not create a computational element\n";
        LOGPZ_WARN(logger, sout.str())
    }
#endif
    
	return NULL;
}


TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
    {
		return new TPZCompElH1<TPZShapePoint>(mesh,gel,h1fam);
    }
    
	return NULL;
}
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
    {
		TPZCompEl *result = new TPZCompElH1<TPZShapeLinear>(mesh,gel,h1fam);
        return result;//new TPZCondensedCompel(result);
    }
    
	return NULL;
}
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
    {
		return new TPZCompElH1<TPZShapeQuad>(mesh,gel,h1fam);
    }
    
	return NULL;
}

TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElH1<TPZShapeTriang>(mesh,gel,h1fam);
    
	return NULL;
}
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElH1<TPZShapeCube>(mesh,gel,h1fam);
    
	return NULL;
}
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElH1<TPZShapePrism>(mesh,gel,h1fam);
    
	return NULL;
}
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElH1<TPZShapePiram>(mesh,gel,h1fam);
    
	return NULL;
}
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,const H1Family h1fam) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElH1<TPZShapeTetra>(mesh,gel,h1fam);
    
	return NULL;
}


// with mem
TPZCompEl *CreatePointElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZCompElH1<TPZShapePoint> >(mesh,gel) ;
    
	return NULL;
}
TPZCompEl *CreateLinearElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZCompElH1<TPZShapeLinear> >(mesh,gel);
    
	return NULL;
}
TPZCompEl *CreateQuadElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZCompElH1<TPZShapeQuad> >(mesh,gel);
    
	return NULL;
}
TPZCompEl *CreateTriangleElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZCompElH1<TPZShapeTriang> >(mesh,gel);
    
	return NULL;
}
TPZCompEl *CreateCubeElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZCompElH1<TPZShapeCube> >(mesh,gel);
    
	return NULL;
}
TPZCompEl *CreatePrismElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZCompElH1<TPZShapePrism> >(mesh,gel);
    
	return NULL;
}
TPZCompEl *CreatePyramElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZCompElH1<TPZShapePiram> >(mesh,gel);
    
	return NULL;
}
TPZCompEl *CreateTetraElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZCompElH1<TPZShapeTetra> >(mesh,gel);
    
	return NULL;
}

/**
 * @brief Creates the computational elements, and the degree of freedom nodes
 */ 
/** Only element of material id in the set<int> will be created */
void TPZCreateApproximationSpace::BuildMesh(TPZCompMesh &cmesh, const TPZVec<int64_t> &gelindexes) const {
    TPZVec<TPZGeoEl *> elvec(gelindexes.size());
    TPZGeoMesh *gmesh = cmesh.Reference();
    for(int64_t el = 0; el<gelindexes.size(); el++)
    {
        elvec[el] = gmesh->Element(gelindexes[el]);
    }
	int64_t i, nelem = elvec.NElements();
	int64_t neltocreate = 0;
	int64_t index;
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
        if (gel->Reference()) {
            continue;
        }
	}
	std::set<int> matnotfound;
	int64_t nbl = cmesh.Block().NBlocks();
    if(neltocreate > nbl) 
    {
        cmesh.Block().SetNBlocks(neltocreate);
    }
    cmesh.Block().SetNBlocks(nbl);
	
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel || gel->Reference()) continue;
        int matid = gel->MaterialId();
        TPZMaterial * mat = cmesh.FindMaterial(matid);
        if(!mat)
        {
            matnotfound.insert(matid);
            continue;
        }
        
        if(!gel->Reference() && gel->NumInterfaces() == 0)
        {
            index = CreateCompEl(gel,cmesh)->Index();
            if (fCreateHybridMesh) {
                cmesh.ElementVec()[index]->Reference()->ResetReference();
            }
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                
                std::stringstream sout;
                if (index < 0) {
                    if(gel->Dimension() == 0){
                        sout << "Zero dimensional element, is your approximation space Hdiv?. " << std::endl;
                        gel->Print(sout);
                        sout << "No computational element was created. " << std::endl;
                    }
                }else{
                    cmesh.ElementVec()[index]->Print(sout);
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
        }
	}
	cmesh.InitializeBlock();
}

/**
 * @brief Creates the computational elements, and the degree of freedom nodes
 */
/** Only element of material id in the set<int> will be created */
void TPZCreateApproximationSpace::BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs) const {
    TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh.Reference()->ElementVec();
    int64_t i, nelem = elvec.NElements();
    int64_t neltocreate = 0;
    int64_t index;
    for(i=0; i<nelem; i++) {
        TPZGeoEl *gel = elvec[i];
        if(!gel) continue;
        if (gel->Reference()) {
            continue;
        }
        if(!gel->HasSubElement()) {
            neltocreate++;
        }
    }
    std::set<int> matnotfound;
    int64_t nbl = cmesh.Block().NBlocks();
    if(neltocreate > nbl)
    {
        cmesh.Block().SetNBlocks(neltocreate);
    }
    cmesh.Block().SetNBlocks(nbl);
    
    for(i=0; i<nelem; i++) {
        TPZGeoEl *gel = elvec[i];
        if(!gel || gel->Reference()) continue;
        if(!gel->HasSubElement()) {
            int matid = gel->MaterialId();
            TPZMaterial * mat = cmesh.FindMaterial(matid);
            if(!mat)
            {
                matnotfound.insert(matid);
                continue;
            }
            int printing = 0;
            if (printing) {
                gel->Print(std::cout);
            }
            
            //checking material in MaterialIDs
            std::set<int>::const_iterator found = MaterialIDs.find(matid);
            if (found == MaterialIDs.end())
            {
                continue;
            }
            
            if(!gel->Reference() && gel->NumInterfaces() == 0)
            {
                index = CreateCompEl(gel,cmesh)->Index();
                if (fCreateHybridMesh) {
                    cmesh.ElementVec()[index]->Reference()->ResetReference();
                }
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    
                    std::stringstream sout;
                    if (index < 0) {
                        if(gel->Dimension() == 0){
                            sout << "Zero dimensional element, is your approximation space Hdiv?. " << std::endl;
                            gel->Print(sout);
                            sout << "No computational element was created. " << std::endl;
                        }
                    }else{
                        cmesh.ElementVec()[index]->Print(sout);
                    }
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
        }
    }
    cmesh.InitializeBlock();
}

/** @brief Creates the computational elements, and the degree of freedom nodes */
void TPZCreateApproximationSpace::BuildMesh(TPZCompMesh &cmesh) const {
    std::set<int> materialids;
    std::map<int,TPZMaterial *>::iterator it;
    for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
        materialids.insert(it->first);
    }
    BuildMesh(cmesh,materialids);
}

void TPZCreateApproximationSpace::CreateInterfaces(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs){
    for(int64_t el = 0; el < cmesh.ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh.ElementVec()[el];
		if(!compEl) continue;
        int matid = compEl->Reference()->MaterialId();
        
        //checking if matid is in MaterialIDs
        std::set<int>::iterator it;
        it = MaterialIDs.find(matid);
        if(it!=MaterialIDs.end())
        {
            TPZInterpolationSpace * interpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
            if(!interpEl) continue;
            interpEl->CreateInterfaces();
        }
    }
    
}

void TPZCreateApproximationSpace::CreateInterfaces(TPZCompMesh &cmesh){
    for(int64_t el = 0; el < cmesh.ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh.ElementVec()[el];
		if(!compEl) continue;
        
        TPZInterpolationSpace * interpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
        if(!interpEl) continue;
        interpEl->CreateInterfaces();
    }
    
}

void TPZCreateApproximationSpace::SetAllCreateFunctions(TPZCompEl &cel, TPZCompMesh *mesh){
	cel.SetCreateFunctions(mesh);
}



void TPZCreateApproximationSpace::SetAllCreateFunctionsDiscontinuous(){
	
    fStyle = EDiscontinuous;
    fp[EPoint] = TPZCompElDisc::CreateDisc;
    fp[EOned] = TPZCompElDisc::CreateDisc;
    fp[ETriangle] = TPZCompElDisc::CreateDisc;
    fp[EQuadrilateral] = TPZCompElDisc::CreateDisc;
    fp[ETetraedro] = TPZCompElDisc::CreateDisc;
    fp[EPiramide] = TPZCompElDisc::CreateDisc;
    fp[EPrisma] = TPZCompElDisc::CreateDisc;
    fp[ECube] = TPZCompElDisc::CreateDisc;
    
    /*
     pzgeom::TPZGeoPoint::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoLinear::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoQuad::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoTriangle::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoPrism::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoTetrahedra::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoPyramid::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoCube::fp = TPZCompElDisc::CreateDisc;
     */
}


void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuous(){
    fStyle = EContinuous;
    const H1Family &h1fam = this->fh1fam;
    fp[EPoint] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreatePointEl(gel,mesh,h1fam);};
    fp[EOned] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateLinearEl(gel,mesh,h1fam);};
    fp[EQuadrilateral] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateQuadEl(gel,mesh,h1fam);};
    fp[ETriangle] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateTriangleEl(gel,mesh,h1fam);};
    fp[ETetraedro] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateTetraEl(gel,mesh,h1fam);};
    fp[EPiramide] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreatePyramEl(gel,mesh,h1fam);};
    fp[EPrisma] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreatePrismEl(gel,mesh,h1fam);};
    fp[ECube] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateCubeEl(gel,mesh,h1fam);};
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreatePointEl;
     pzgeom::TPZGeoLinear::fp =  CreateLinearEl;
     pzgeom::TPZGeoQuad::fp = CreateQuadEl;
     pzgeom::TPZGeoTriangle::fp =  CreateTriangleEl;
     pzgeom::TPZGeoPrism::fp = CreatePrismEl;
     pzgeom::TPZGeoTetrahedra::fp = CreateTetraEl;
     pzgeom::TPZGeoPyramid::fp = CreatePyramEl;
     pzgeom::TPZGeoCube::fp = CreateCubeEl;
     */
	
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuousWithMem()
{
    fStyle = EContinuous;
	// These ones will be always continuous for viscoelasticity
	fp[EPoint] = CreatePointElWithMem;
	fp[EOned] = CreateLinearElWithMem;
	fp[ETriangle] = CreateTriangleElWithMem;
	fp[EQuadrilateral] = CreateQuadElWithMem;	
    fp[ETetraedro] = CreateTetraElWithMem;
    fp[EPiramide] = CreatePyramElWithMem;
    fp[EPrisma] = CreatePrismElWithMem;
    fp[ECube] = CreateCubeElWithMem;
    
}


#include "pzelchdiv.h"
#include "pzelchdivbound2.h"
#include "TPZCompElKernelHDiv3D.h"
#include "TPZCompElKernelHDiv.h"
#include "TPZCompElHDivDuplConnects.h"
#include "TPZCompElHDivDuplConnectsBound.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDiv(int dimension){

    fStyle = EHDiv;
    const HDivFamily &hdivfam = this->fhdivfam;

    if (fhdivfam == HDivFamily::EHDivKernel){
        switch (dimension) {
            case 1:
                std::cout << "HDivKernel family not implemented for 1D problems" << std::endl;
                DebugStop();//
                break;
            case 2:
                fp[EPoint] = CreateNoElement;
                fp[EOned] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelBoundLinearEl(gel,mesh,hdivfam);};
                fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelTriangleEl(gel,mesh,hdivfam);};
                fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelQuadEl(gel,mesh,hdivfam);};
                fp[ETetraedro] = CreateNoElement;
                fp[EPiramide] = CreateNoElement;
                fp[EPrisma] = CreateNoElement;
                fp[ECube] = CreateNoElement;
                break;
            case 3:
                fp[EPoint] = CreateNoElement;
                fp[EOned] = CreateNoElement;
                fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelBoundTriangleEl(gel,mesh,hdivfam);};
                fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelBoundQuadEl(gel,mesh,hdivfam);};
                fp[ETetraedro] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelTetraEl(gel,mesh,hdivfam);};
                fp[EPiramide] = CreateNoElement;
                fp[EPrisma] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelPrismEl(gel,mesh,hdivfam);};
                fp[ECube] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivKernelCubeEl(gel,mesh,hdivfam);};
                break;
            default:
                DebugStop();
                break;
        }
    } else { // EHDivStandard or EHDivConstant
        switch (dimension) {
            case 1:
                fp[EPoint] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundPointEl(gel,mesh,hdivfam);};
                fp[EOned] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivLinearEl(gel,mesh,hdivfam);};
                fp[ETriangle] = CreateNoElement;
                fp[EQuadrilateral] = CreateNoElement;
                fp[ETetraedro] = CreateNoElement;
                fp[EPiramide] = CreateNoElement;
                fp[EPrisma] = CreateNoElement;
                fp[ECube] = CreateNoElement;
                break;
            case 2:
                fp[EPoint] = CreateNoElement;
                fp[EOned] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundLinearEl(gel,mesh,hdivfam);};
                fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivTriangleEl(gel,mesh,hdivfam);};
                fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivQuadEl(gel,mesh,hdivfam);};
                fp[ETetraedro] = CreateNoElement;
                fp[EPiramide] = CreateNoElement;
                fp[EPrisma] = CreateNoElement;
                fp[ECube] = CreateNoElement;
                break;
            case 3:
                fp[EPoint] = CreateNoElement;
                fp[EOned] = CreateNoElement;
                fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundTriangleEl(gel,mesh,hdivfam);};
                fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundQuadEl(gel,mesh,hdivfam);};
                fp[ETetraedro] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivTetraEl(gel,mesh,hdivfam);};
                fp[EPiramide] = CreateNoElement;
                fp[EPrisma] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPrismEl(gel,mesh,hdivfam);};
                fp[ECube] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivCubeEl(gel,mesh,hdivfam);};
                break;
            default:
                DebugStop();
                break;
        }
    }
    
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreateHDivPointEl;
     pzgeom::TPZGeoLinear::fp =  CreateHDivLinearEl;
     pzgeom::TPZGeoQuad::fp = CreateHDivQuadEl;
     pzgeom::TPZGeoTriangle::fp =  CreateHDivTriangleEl;
     pzgeom::TPZGeoPrism::fp = CreateHDivPrismEl;
     pzgeom::TPZGeoTetrahedra::fp = CreateHDivTetraEl;
     pzgeom::TPZGeoPyramid::fp = CreateHDivPyramEl;
     pzgeom::TPZGeoCube::fp = CreateHDivCubeEl;
     */
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDivDuplConnects(int dimension){

    fStyle = EHDiv;
    const HDivFamily &hdivfam = this->fhdivfam;

    if (hdivfam != HDivFamily::EHDivConstant && hdivfam != HDivFamily::EHDivOptimized){
        std::cout << "HDiv Dupl Connects not implemented yet for this HDiv family!" << std::endl;
        DebugStop();
    }

    switch (dimension) {
        case 1:
            fp[EPoint] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsBoundPointEl(gel,mesh,hdivfam);};
            fp[EOned] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsLinearEl(gel,mesh,hdivfam);};
            fp[ETriangle] = CreateNoElement;
            fp[EQuadrilateral] = CreateNoElement;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 2:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsBoundLinearEl(gel,mesh,hdivfam);};
            fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsTriangleEl(gel,mesh,hdivfam);};
            fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsQuadEl(gel,mesh,hdivfam);};
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateNoElement;
            fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsBoundTriangEl(gel,mesh,hdivfam);};
            fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsBoundQuadEl(gel,mesh,hdivfam);};
            fp[ETetraedro] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsTetraEl(gel,mesh,hdivfam);};
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivDuplConnectsCubeEl(gel,mesh,hdivfam);};
            break;
        default:
            DebugStop();
            break;
    }

}

#include <TPZCompElHCurl.h>

#define HCURL_EL_NOT_AVAILABLE \
    PZError<<__PRETTY_FUNCTION__;\
    PZError<<"Element not available.\n";\
    PZError<<"Aborting...\n";\
    DebugStop();\
    return nullptr;

void TPZCreateApproximationSpace::SetAllCreateFunctionsHCurl(int dimension){

    fStyle = EHCurl;
    const HCurlFamily &hcurlfam = this->fhcurlfam;

    
    switch (dimension) {
    case 1:
        fp[EPoint] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            HCURL_EL_NOT_AVAILABLE
        };
        fp[EOned] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapeLinear>(mesh,gel,hcurlfam);
        };
        fp[ETriangle] = CreateNoElement;
        fp[EQuadrilateral] = CreateNoElement;
        fp[ETetraedro] = CreateNoElement;
        fp[EPiramide] = CreateNoElement;
        fp[EPrisma] = CreateNoElement;
        fp[ECube] = CreateNoElement;
        break;
    case 2:
        fp[EPoint] = CreateNoElement;
        fp[EOned] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapeLinear>(mesh,gel,hcurlfam);
        };
        fp[ETriangle] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return  new TPZCompElHCurl<pzshape::TPZShapeTriang>(mesh, gel, hcurlfam);
        };
        fp[EQuadrilateral] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapeQuad>(mesh,gel,hcurlfam);
        };
        fp[ETetraedro] = CreateNoElement;
        fp[EPiramide] = CreateNoElement;
        fp[EPrisma] = CreateNoElement;
        fp[ECube] = CreateNoElement;
        break;
    case 3:
        fp[EPoint] = CreateNoElement;
        fp[EOned] = CreateNoElement;
        fp[ETriangle] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapeTriang>(mesh,gel,hcurlfam);
        };
        fp[EQuadrilateral] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapeQuad>(mesh,gel,hcurlfam);
        };
        fp[ETetraedro] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapeTetra>(mesh,gel,hcurlfam);
        };
        fp[EPiramide] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            HCURL_EL_NOT_AVAILABLE
                };
        fp[EPrisma] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapePrism>(mesh,gel,hcurlfam);
        };
        fp[ECube] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            return new TPZCompElHCurl<pzshape::TPZShapeCube>(mesh,gel,hcurlfam);
        };
        break;
    default:
        DebugStop();
        break;
    }
    
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsHCurlWithMem(int dimension){

    fStyle = EHCurl;
    const HCurlFamily &hcurlfam = this->fhcurlfam;

    
    switch (dimension) {
    case 1:
        fp[EPoint] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            HCURL_EL_NOT_AVAILABLE
        };
        fp[EOned] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeLinear>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[ETriangle] = CreateNoElement;
        fp[EQuadrilateral] = CreateNoElement;
        fp[ETetraedro] = CreateNoElement;
        fp[EPiramide] = CreateNoElement;
        fp[EPrisma] = CreateNoElement;
        fp[ECube] = CreateNoElement;
        break;
    case 2:
        fp[EPoint] = CreateNoElement;
        fp[EOned] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeLinear>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[ETriangle] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                  new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeTriang>>(mesh, gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[EQuadrilateral] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                 new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeQuad>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[ETetraedro] = CreateNoElement;
        fp[EPiramide] = CreateNoElement;
        fp[EPrisma] = CreateNoElement;
        fp[ECube] = CreateNoElement;
        break;
    case 3:
        fp[EPoint] = CreateNoElement;
        fp[EOned] = CreateNoElement;
        fp[ETriangle] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                 new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeTriang>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[EQuadrilateral] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                 new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeQuad>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[ETetraedro] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                 new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeTetra>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[EPiramide] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            HCURL_EL_NOT_AVAILABLE
                };
        fp[EPrisma] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                 new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapePrism>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        fp[ECube] = [hcurlfam](TPZGeoEl *gel,TPZCompMesh &mesh) {
            auto *cel =
                 new TPZCompElWithMem<TPZCompElHCurl<pzshape::TPZShapeCube>>(mesh,gel);
            cel->SetFamily(hcurlfam);
            return cel;
        };
        break;
    default:
        DebugStop();
        break;
    }
    
}

#undef HCURL_EL_NOT_AVAILABLE

#if defined(USING_LAPACK) && !defined(STATE_COMPLEX)

#include "TPZSBFemVolume.h"
#include "TPZSBFemVolumeMultiphysics.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsSBFem(int dimension){
    
    fStyle = ESBFem;
    
    // NOTE: I dont think we support multiple h1 family spaces for SBFEM. Please pass as argument to function in case it is needed.
    const H1Family &h1fam = this->fh1fam;
    switch (dimension) {
        case 1:
            DebugStop();
            break;
        case 2:
            fp[EPoint] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreatePointEl(gel,mesh,h1fam);};
            fp[EOned] = CreateSBFemCompEl;
            fp[ETriangle] = CreateNoElement;
            fp[EQuadrilateral] = CreateSBFemCompEl;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreatePointEl(gel,mesh,h1fam);};
            fp[EOned] = CreateSBFemCompEl;
            fp[ETriangle] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateTriangleEl(gel,mesh,h1fam);};
            fp[EQuadrilateral] = [h1fam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateQuadEl(gel,mesh,h1fam);};
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateSBFemCompEl;
            fp[ECube] = CreateSBFemCompEl;
            break;
        default:
            DebugStop();
            break;
    }
}
void TPZCreateApproximationSpace::SetAllCreateFunctionsSBFemMultiphysics(int dimension){
    fStyle = EMultiphysicsSBFem;
    switch (dimension) {
        case 1:
            DebugStop();
            break;
        case 2:
            fp[EPoint] = CreateMultiphysicsPointEl;
            fp[EOned] = CreateMultiphysicsLinearEl;
            fp[ETriangle] = CreateMultiphysicsTriangleEl;
            fp[EQuadrilateral] = CreateSBFemMultiphysicsQuadEl; // SBFem
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateNoElement;
            fp[ETriangle] = CreateNoElement;
            fp[EQuadrilateral] = CreateNoElement;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateSBFemMultiphysicsPrismaEl;
            fp[ECube] = CreateSBFemMultiphysicsCubeEl;
            break;
        default:
            DebugStop();
            break;
    }
}
#else
void TPZCreateApproximationSpace::SetAllCreateFunctionsSBFem(int dimension){
    PZError<<__PRETTY_FUNCTION__<<" depends on LAPACK!\n";
    PZError<<"Please configure NeoPZ with USING_LAPACK or USING_MKL =ON\n";
    PZError<<"Aborting..."<<std::endl;
    DebugStop();
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsSBFemMultiphysics(int dimension){
    PZError<<__PRETTY_FUNCTION__<<" depends on LAPACK!\n";
    PZError<<"Please configure NeoPZ with USING_LAPACK or USING_MKL =ON\n";
    PZError<<"Aborting..."<<std::endl;
    DebugStop();
}
#endif

////#include "pzhdivfull.h"
//
//void TPZCreateApproximationSpace::SetAllCreateFunctionsHDivFull(int dimension){
//    
//    switch (dimension) {
//        case 1:
//            fp[EPoint] = CreateHDivBoundPointEl;
//            fp[EOned] = CreateHDivFullLinearEl;
//            fp[ETriangle] = CreateNoElement;
//            fp[EQuadrilateral] = CreateNoElement;
//            fp[ETetraedro] = CreateNoElement;
//            fp[EPiramide] = CreateNoElement;
//            fp[EPrisma] = CreateNoElement;
//            fp[ECube] = CreateNoElement;
//            break;
//        case 2:
//            fp[EPoint] = CreateNoElement;
//            fp[EOned] = CreateHDivBoundLinearEl;
//            fp[ETriangle] = CreateHDivFullTriangleEl;
//            fp[EQuadrilateral] = CreateHDivFullQuadEl;
//            fp[ETetraedro] = CreateNoElement;
//            fp[EPiramide] = CreateNoElement;
//            fp[EPrisma] = CreateNoElement;
//            fp[ECube] = CreateNoElement;
//            break;
//        case 3:
//            fp[EPoint] = CreateNoElement;
//            fp[EOned] = CreateNoElement;
//            fp[ETriangle] = CreateHDivBoundTriangleEl;
//            fp[EQuadrilateral] = CreateHDivBoundQuadEl;
//            fp[ETetraedro] = CreateHDivFullTetraEl;
//            fp[EPiramide] = CreateHDivFullPyramEl;
//            fp[EPrisma] = CreateHDivFullPrismEl;
//            fp[ECube] = CreateHDivFullCubeEl;
//            break;
//        default:
//            DebugStop();
//            break;
//    }
//    
//}



#ifndef STATE_COMPLEX

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDivPressure(int dimension){
    // This function has been deprecated since TPZCompElHDivPressure has been deprecated (Jan 31/2022)
    // TODO: Delete me in the future?
    /*
    fStyle = EHDiv;
    switch (dimension) {
        case 1:
            fp[EPoint] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundPointEl(gel,mesh,hdivfam);};
            fp[EOned] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPressureLinearEl(gel,mesh,hdivfam);};
            fp[ETriangle] = CreateNoElement;
            fp[EQuadrilateral] = CreateNoElement;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 2:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundLinearEl(gel,mesh,hdivfam);};
            fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPressureTriangleEl(gel,mesh,hdivfam);};
            fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPressureQuadEl(gel,mesh,hdivfam);};
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateNoElement;
            fp[ETriangle] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundTriangleEl(gel,mesh,hdivfam);};
            fp[EQuadrilateral] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivBoundQuadEl(gel,mesh,hdivfam);};
            fp[ETetraedro] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPressureTetraEl(gel,mesh,hdivfam);};
            fp[EPiramide] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPressurePyramEl(gel,mesh,hdivfam);};
            fp[EPrisma] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPressurePrismEl(gel,mesh,hdivfam);};
            fp[ECube] = [hdivfam](TPZGeoEl *gel,TPZCompMesh &mesh) {return CreateHDivPressureCubeEl(gel,mesh,hdivfam);};
            break;
        default:
            DebugStop();
            break;
    }
     */
    
}
#endif

#include "pzmultiphysicscompel.h"
void TPZCreateApproximationSpace::SetAllCreateFunctionsMultiphysicElem(){
	
    fStyle = EMultiphysics;
    fp[EPoint] = CreateMultiphysicsPointEl;
    fp[EOned] = CreateMultiphysicsLinearEl;
    fp[ETriangle] = CreateMultiphysicsTriangleEl;
    fp[EQuadrilateral] = CreateMultiphysicsQuadEl;
    fp[ETetraedro] = CreateMultiphysicsTetraEl;
    fp[EPiramide] = CreateMultiphysicsPyramEl;
    fp[EPrisma] = CreateMultiphysicsPrismEl;
    fp[ECube] = CreateMultiphysicsCubeEl;
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreateMultiphysicsPointEl;
     pzgeom::TPZGeoLinear::fp =  CreateMultiphysicsLinearEl;
     pzgeom::TPZGeoTriangle::fp =  CreateMultiphysicsTriangleEl;
     pzgeom::TPZGeoQuad::fp = CreateMultiphysicsQuadEl;
     pzgeom::TPZGeoCube::fp = CreateMultiphysicsCubeEl;
     pzgeom::TPZGeoPrism::fp = CreateMultiphysicsPrismEl;
     pzgeom::TPZGeoTetrahedra::fp = CreateMultiphysicsTetraEl;
     pzgeom::TPZGeoPyramid::fp = CreateMultiphysicsPyramEl;
     */
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsMultiphysicElemWithMem()
{

    fStyle = EMultiphysics;
	fp[EPoint] = CreateMultiphysicsPointElWithMem;
	fp[EOned] = CreateMultiphysicsLinearElWithMem;
	fp[ETriangle] = CreateMultiphysicsTriangleElWithMem;
	fp[EQuadrilateral] = CreateMultiphysicsQuadElWithMem;
	fp[ETetraedro] = CreateMultiphysicsTetraElWithMem;
	fp[EPiramide] = CreateMultiphysicsPyramElWithMem;
	fp[EPrisma] = CreateMultiphysicsPrismElWithMem;
	fp[ECube] = CreateMultiphysicsCubeElWithMem;
	
}

/*
 * @brief Create a computational element using the function pointer for the topology
 */
TPZCompEl *TPZCreateApproximationSpace::CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh) const
{
    switch (gel->Type()) {
        case EPoint:
            return fp[EPoint](gel,mesh);
            break;
        case EOned:
            return fp[EOned](gel,mesh);
            break;
        case EQuadrilateral:
            return fp[EQuadrilateral](gel,mesh);
            break;
        case ETriangle:
            return fp[ETriangle](gel,mesh);
            break;
        case EPiramide:
            return fp[EPiramide](gel,mesh);
            break;
        case EPrisma:
            return fp[EPrisma](gel,mesh);
            break;
        case ETetraedro:
            return fp[ETetraedro](gel,mesh);
            break;
        case ECube:
            return fp[ECube](gel,mesh);
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

/**
 * @brief Set custom function pointers
 */

void TPZCreateApproximationSpace::SetCreateFunctions(TPZVec<TCreateFunction> &createfuncs)
{
    fStyle = ECustom;
    fp[EPoint] = createfuncs[EPoint];
    fp[EOned] = createfuncs[EOned];
    fp[ETriangle] = createfuncs[ETriangle];
    fp[EQuadrilateral] = createfuncs[EQuadrilateral];
    fp[ETetraedro] = createfuncs[ETetraedro];
    fp[EPiramide] = createfuncs[EPiramide];
    fp[EPrisma] = createfuncs[EPrisma];
    fp[ECube] = createfuncs[ECube];
    
}

#include "pzcondensedcompel.h"

/**
 * @brief Encapsulate the elements in condensed computational elements
 */
void TPZCreateApproximationSpace::CondenseLocalEquations(TPZCompMesh &cmesh)
{
    const bool real_sol = cmesh.GetSolType()!=ESolType::EComplex;

    int64_t nel = cmesh.NElements();
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh.ElementVec()[iel];
        if(!cel) {
            continue;
        }
        if(real_sol){
            new TPZCondensedCompElT<STATE>(cel);
        }else{
            new TPZCondensedCompElT<CSTATE>(cel);
        }
    }
    
}

/**
 * @brief Undo the encapsulate elements
 */
void TPZCreateApproximationSpace::UndoCondenseLocalEquations(TPZCompMesh &cmesh)
{
    int64_t nel = cmesh.NElements();
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh.ElementVec()[iel];
        TPZCondensedCompEl *condel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if(!condel) {
            continue;
        }
        condel->Unwrap();
    }    
}

/**
 * @brief transform in low order Raviar Tomas
 */
void TPZCreateApproximationSpace::MakeRaviartThomas(TPZCompMesh &cmesh)
{
    int64_t numcell = cmesh.NElements();
    int64_t el;
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            continue;
        }
        intel->SetPreferredOrder(1);
    }
    cmesh.ExpandSolution();
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        int geldim = gel->Dimension();
        int is;
        int nsides = gel->NSides();
        for (is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != geldim-1) {
                continue;
            }
            int nsconnects = intel->NSideConnects(is);
            // only interested in HDiv elements
            if (nsconnects != 1) {
                continue;
            }
            int64_t cindex = intel->SideConnectIndex(0, is);
            TPZConnect &c = intel->Connect(intel->SideConnectLocId(0,is));
            if (c.HasDependency()) {
                continue;
            }
            int nshape = 1;
            int nstate = 1;
            int order = 0;
            int64_t cindex2 = cmesh.AllocateNewConnect(nshape, nstate, order);
            //            TPZConnect &c2 = cmesh.ConnectVec()[cindex];
            TPZFNMatrix<2,REAL> depmat(2,1,1.);
            c.AddDependency(cindex, cindex2, depmat, 0, 0, 2, 1);
        }
    }
    cmesh.ExpandSolution();
}

/**
 * @brief transform in low order Raviar Tomas
 */
void TPZCreateApproximationSpace::UndoMakeRaviartThomas(TPZCompMesh &cmesh)
{
    int64_t numcell = cmesh.NElements();
    int64_t el;
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        int geldim = gel->Dimension();
        int is;
        int nsides = gel->NSides();
        for (is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != geldim-1) {
                continue;
            }
            int nsconnects = intel->NSideConnects(is);
            // only interested in HDiv elements
            if (nsconnects != 1) {
                continue;
            }
            //            int cindex = intel->SideConnectIndex(0, is);
            TPZConnect &c = intel->Connect(intel->SideConnectLocId(0,is));
            if (c.HasDependency()) {
                c.RemoveDepend();
            }
        }
    }
    cmesh.ExpandSolution();
    cmesh.CleanUpUnconnectedNodes();
}

/** @brief Create interface elements between the computational elements */
void TPZCreateApproximationSpace::CreateInterfaceElements(TPZCompMesh *mesh, bool betweencontinuous, bool multiphysics)
{
    TPZChunkVector<TPZCompEl *> compelvec = mesh->ElementVec();
    int64_t nel = compelvec.NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = compelvec[el];
        if (!cel) {
            continue;
        }
        if(!multiphysics)
        {
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (intel) {
                intel->CreateInterfaces(betweencontinuous);
            }
        }
        else {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mfcel) {
                continue;
            }
            mfcel->CreateInterfaces();
        }
        
    }
}

/// this method will substitute all interface elements with materialid within the set by three elements : one H1 element and two interface elements
void TPZCreateApproximationSpace::Hybridize(TPZCompMesh &cmesh,const std::set<int> &matids, bool isconnectedElem)
{
    cmesh.ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh.Reference()->ResetReference();
    int nel = cmesh.NElements();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!face) {
            continue;
        }
        int matid = face->Material()->Id();
        if (matids.find(matid) == matids.end()) {
            continue;
        }
        TPZCompElSide left, right;
        left = face->LeftElementSide();
        right = face->RightElementSide();
        int leftmatid = left.Element()->Material()->Id();
        int rightmatid = right.Element()->Material()->Id();
        TPZGeoEl *gelface = face->Reference();
        gelface = gelface->CreateBCGeoEl(gelface->NSides()-1, matid);
        delete face;
        
        //hp mesh
        TPZInterpolationSpace *leftint = dynamic_cast<TPZInterpolationSpace *>(left.Element());
        TPZInterpolationSpace *rightint = dynamic_cast<TPZInterpolationSpace *>(right.Element());
        int leftorder = leftint->MaxOrder();
        int rightorder = rightint->MaxOrder();
        int neworder =  MIN(leftorder,rightorder);
        if(neworder<0) neworder = MAX(leftorder,rightorder);
        cmesh.SetDefaultOrder(neworder);
//        int neworder = MAX(leftorder,rightorder);
//        cmesh.SetDefaultOrder(neworder);
       
        
        TPZCompEl *newcel = cmesh.ApproxSpace().CreateCompEl(gelface, cmesh);
//        TPZCompEl *newcel = cmesh.ElementVec()[index];
        if(!isconnectedElem){
            gelface->ResetReference();
        }
        TPZCompElSide center(newcel,gelface->NSides()-1);
        //        
        TPZGeoEl *leftgelface = gelface->CreateBCGeoEl(gelface->NSides()-1, leftmatid);
        TPZGeoEl *rightgelface = gelface->CreateBCGeoEl(gelface->NSides()-1, rightmatid);

        new TPZInterfaceElement(cmesh,leftgelface,left,center);
        new TPZInterfaceElement(cmesh,rightgelface,right,center);
        
//        TPZInterfaceElement *faceleft = new TPZInterfaceElement(cmesh,leftgelface,index,left,center);
//        TPZInterfaceElement *faceright = new TPZInterfaceElement(cmesh,rightgelface,index,right,center);
        
    }
    cmesh.CleanUpUnconnectedNodes();
    cmesh.ExpandSolution();
}

int TPZCreateApproximationSpace::ClassId() const {
    return Hash("TPZCreateApproximationSpace");
}

void TPZCreateApproximationSpace::Read(TPZStream& buf, void* context) { //ok
    buf.Read(fCreateHybridMesh);
    buf.Read(fCreateLagrangeMultiplier);
    buf.Read(fCreateWithMemory);
}

void TPZCreateApproximationSpace::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(fCreateHybridMesh);
    buf.Write(fCreateLagrangeMultiplier);
    buf.Write(fCreateWithMemory);
}

template class TPZRestoreClass<TPZCreateApproximationSpace>;

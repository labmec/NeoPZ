#include "TPZCompElHDivDuplConnectsBound.h"
#include "TPZMaterial.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstantBound.h"
#include "TPZShapeHDivOptimized.h"
#include "pzlog.h"
#include "pzcmesh.h"
#include <sstream>

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

template<class TSHAPE>
TPZCompElHDivDuplConnectsBound<TSHAPE>::TPZCompElHDivDuplConnectsBound(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivDuplConnectsBound::ClassId), TPZCompElHDivBound2<TSHAPE>(mesh,gel,hdivfam){
    
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::NConnects() const {
    if (fDuplicationActive){
        return this->fConnectIndexes.size();
    }else{
        return TPZCompElHDivBound2<TSHAPE>::NConnects();
    }
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::NSideConnects(int side) const{
	if(side == TSHAPE::NSides-1)
	{
		if (fDuplicationActive){
            return 2;
        } else {
            return TPZCompElHDivBound2<TSHAPE>::NSideConnects(side);
        };
	}
	return 0;
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::NConnectShapeF(int connect, int connectorder)const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif

    if(connectorder == 0){
        if (connect == 0){
            return 1;
        }else {
            return 0;
        }
    }
     

	switch (this->fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        {
            
            TPZManVector<int,22> order(TSHAPE::NSides-TSHAPE::NCornerNodes,connectorder);
            int nshape = TSHAPE::NShapeF(order);
            if (!fDuplicationActive) return nshape;
            int conCorrect = connect/2;
            int res = connect % 2;
            if (res == 1){ 
                nshape -= 1;
            } else {
                nshape = 1;
            }

            return nshape;
        }    
        break;
    case HDivFamily::EHDivConstant:
    case HDivFamily::EHDivOptimized:
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape = TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(connect,connectorder);
            if (!fDuplicationActive) return nshape;
            if (res == 1){ 
                nshape -= 1;
            } else {
                nshape = 1;
            }

            return nshape;
        }
        break;
    
    default:
        DebugStop();//You shoud choose an approximation space
        break;
    }
    
    return -1;
 }

template<class TSHAPE>
int64_t TPZCompElHDivDuplConnectsBound<TSHAPE>::ConnectIndex(int con) const{
// #ifndef PZNODEBUG
// 	if(con<0 || con > 2) {
// 		std::cout << "TPZCompElHDivDuplConnectsBound::ConnectIndex wrong parameter connect " << con <<
// 		" NConnects " << TSHAPE::NFacets << std::endl;
// 		DebugStop();
// 		return -1;
// 	}

// #endif
    if (fDuplicationActive){
        return this->fConnectIndexes[con];
    } else {
        return TPZCompElHDivBound2<TSHAPE>::ConnectIndex(con);
    }  
	
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::SideConnectLocId(int node, int side) const
{
    if (fDuplicationActive){
	    if(side == TSHAPE::NSides-1 && node <2){
		    return node;
	    }else{
	        return -1;
	    }
    }else{
        return TPZCompElHDivBound2<TSHAPE>::SideConnectLocId(node,side);
    }	
}

template<class TSHAPE>
void TPZCompElHDivDuplConnectsBound<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
    if (fDuplicationActive){
        this->fConnectIndexes[i] = connectindex;
    } else {
        TPZCompElHDivBound2<TSHAPE>::SetConnectIndex(i,connectindex);
    }  

}

template<class TSHAPE>
void TPZCompElHDivDuplConnectsBound<TSHAPE>::ActiveDuplConnects(std::map<int64_t,int64_t> &fConnDuplicated){

    if (fDuplicationActive == true) return;
    
    int conOrder = this->ConnectOrder(0);

    fDuplicationActive = true;

    this->fConnectIndexes.Resize(2);

    auto nFacets = this->Reference()->NSides(this->Dimension());
    //Loop over the element facets - which are the connects the be duplicated (edges in 2D and faces in 3D)
    for (int iFacet = 0; iFacet < nFacets; iFacet++)
    {
        // Algorithm description: for each element facet, checks if the corresponding original connect is in the map fConnDuplicated.
        // If Yes, just sets the returning value from fConnDuplicated to the duplicated connect in the current element;
        // If No, allocate a new connect and inserts its index to fConnDuplicated using the original connect as key

        int64_t conn = ConnectIndex(2*iFacet);           

        if (fConnDuplicated.find(conn) == fConnDuplicated.end()){
            //not found, so allocate a new connect
            auto pOrder = conOrder;
            int nshape = 0;//It is updated in the next loop
            int nstate = 1;//It can possibly change
            int64_t newConnect = this->Mesh()->AllocateNewConnect(nshape,nstate,pOrder);
            fConnDuplicated[conn] = newConnect;
            SetConnectIndex(2*iFacet+1,newConnect);
        } else {
            //found, so just set the proper index of the duplicated connect
            SetConnectIndex(2*iFacet+1,fConnDuplicated[conn]);
        }
    }

    TPZInterpolatedElement *celHybBound = dynamic_cast<TPZInterpolatedElement *> (this);
    if (celHybBound){
        for (int icon = 0; icon < celHybBound->NConnects(); icon++)
        {
            TPZConnect &c = celHybBound->Connect(icon);
            int nShapeF = celHybBound->NConnectShapeF(icon,c.Order());
            c.SetNShape(nShapeF);
            int64_t seqnum = c.SequenceNumber();
            int nvar = 1;
            TPZMaterial * mat = celHybBound->Material();
            if (mat) nvar = mat->NStateVariables();
            c.SetNState(nvar);
            celHybBound->Mesh()->Block().Set(seqnum, nvar * nShapeF);
        }
    }

}

template<class TSHAPE>
void TPZCompElHDivDuplConnectsBound<TSHAPE>::InactiveDuplConnects(){

    if (fDuplicationActive == false) return;

    TPZConnect &c2 = this->Connect(1);
    c2.SetNShape(0);

    fDuplicationActive = false;
    this->fConnectIndexes.Resize(1);
    TPZConnect &c = this->Connect(0);
    int nshape = NConnectShapeF(0,c.Order());
    // nshape++;
    c.SetNShape(nshape);
    
    // int64_t seqnum = c.SequenceNumber();
    // int nvar = 1;
    // TPZMaterial * mat = this->Material();
    // if (mat) nvar = mat->NStateVariables();
    // c.SetNState(nvar);
    // this->Mesh()->Block().Set(seqnum, nvar * nshape);
    
}


#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
using namespace pzshape;

template class TPZCompElHDivDuplConnectsBound<TPZShapePoint>;
template class TPZCompElHDivDuplConnectsBound<TPZShapeLinear>;
template class TPZCompElHDivDuplConnectsBound<TPZShapeQuad>;
template class TPZCompElHDivDuplConnectsBound<TPZShapeTriang>;




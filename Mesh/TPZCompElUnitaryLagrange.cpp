#include "TPZCompElUnitaryLagrange.h"
#include "pzcompel.h"
#include "pzcmesh.h"
#include "TPZElementMatrixT.h"

#include "pzmultiphysicselement.h"

TPZCompElUnitaryLagrange::TPZCompElUnitaryLagrange(TPZCompMesh &mesh, TPZGeoEl *reference, TPZCompElSide &wrapSide, TPZCompElSide &lagrangeSide, bool allocateNewConnect) :
    TPZCompElDisc(mesh, reference) {
    this->fConnectIndexes.resize(2);

    //Initializes the connect indexes structure - set the first connect as the same as Lagrange element and the second as the wrap element
    fConnectIndexes[0] = wrapSide.Element()->ConnectIndex(0);
    fConnectIndexes[1] = lagrangeSide.Element()->ConnectIndex(0);

    TPZConnect &c = wrapSide.Element()->Connect(0);
    if (c.HasDependency() && allocateNewConnect){
        int64_t newConnect = this->Mesh()->AllocateNewConnect(1,c.NState(),c.Order());
        fConnectIndexes[0] = newConnect;     

        TPZGeoElSide volside = wrapSide.Reference().operator--();
        int nSides = volside.Element()->NSides();
        int nFacets = volside.Element()->NSides(volside.Element()->Dimension()-1);
        int effectiveSide = volside.Side() - (nSides - nFacets - 1);
        volside.Element()->Reference()->SetConnectIndex(2*effectiveSide,newConnect);
        this->Mesh()->InitializeBlock();
    }
   
}

int TPZCompElUnitaryLagrange::NConnects() const {
	return 2;
}

int TPZCompElUnitaryLagrange::NSideConnects(int side) const{

	return -1;

}

int TPZCompElUnitaryLagrange::NConnectShapeF(int connect, int order)const
{
    return 1;
}


int64_t TPZCompElUnitaryLagrange::ConnectIndex(int con) const{
    if (con < 0 || con >=2) DebugStop();
	return this->fConnectIndexes[con];
}

int TPZCompElUnitaryLagrange::SideConnectLocId(int node,int side) const {
   
    return -1;
}


void TPZCompElUnitaryLagrange::SetConnectIndex(int i, int64_t connectindex)
{
	this->fConnectIndexes[i] = connectindex;
}

template<class TVar>
void TPZCompElUnitaryLagrange::CalcStiffInternal(TPZElementMatrixT<TVar> &ek,TPZElementMatrixT<TVar> &ef){
    this->InitializeElementMatrix(ek, ef);

    //Compute the stiffness matrix (constant Lagrange multiplier), scaled by the dependency matrix.
    //In this case it corresponds to a scalar relating the areas of the wrap and large elements
    int fSize = ek.fMat.Rows();
    int nvar = fSize/2;
    REAL val = 1.;
    for (int i = 0; i < nvar; i++){
        // if (i == nvar) val = -1.;
        ek.fMat(i,nvar + i) = (TVar)fSideOrient;
        ek.fMat(nvar + i,i) = (TVar)fSideOrient;
    }

    // std::cout << "Unitary Lagrange Matrix = " << ek.fMat << std::endl;

    ef.fMat.Zero();
}


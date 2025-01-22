//
// Created by Gustavo A. Batistela on 01/06/2021.
//

#include "TPZNullMaterialCS.h"
#include "TPZMaterialDataT.h"

template<class TVar>
int TPZNullMaterialCS<TVar>::ClassId() const {
    return Hash("TPZNullMaterial") ^ TBase::ClassId() << 1;
}

template<class TVar>
void TPZNullMaterialCS<TVar>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    if (fDim < 1 || fDim > 3) {
        DebugStop();
    }
    buf.Write(&fDim);
    buf.Write(&fNState);
}

template<class TVar>
void TPZNullMaterialCS<TVar>::Read(TPZStream &buf, void *context) {
    TBase::Read(buf, context);
    buf.Read(&fDim);
    buf.Read(&fNState);
}


template<class TVar>
void TPZNullMaterialCS<TVar>::FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const {
    for (auto i = 0; i < datavec.size(); i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fActiveApproxSpace = ENotActive;
    }
}

template<class TVar>
void TPZNullMaterialCS<TVar>::SetDimension(int dim) {
#ifdef PZDEBUG
    if (fDim < 1 || fDim > 3) { DebugStop(); }
#endif
    fDim = dim;
}

template<class TVar>
void TPZNullMaterialCS<TVar>::Print(std::ostream &out) const {
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMaterial::Print(out);
    out << "Dimension " << fDim << std::endl;
    out << "NState " << fNState << std::endl;
}

template class TPZNullMaterialCS<STATE>;
template class TPZNullMaterialCS<CSTATE>;

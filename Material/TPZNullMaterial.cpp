//
// Created by Gustavo Batistela on 5/4/21.
//

#include "TPZNullMaterial.h"
#include "TPZStream.h"
#include "TPZMaterialDataT.h"

template<class TVar>
int TPZNullMaterial<TVar>::ClassId() const {
    return Hash("TPZNullMaterial") ^ TBase::ClassId() << 1;
}

template<class TVar>
void TPZNullMaterial<TVar>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    if (fDim < 1 || fDim > 3) {
        DebugStop();
    }
    buf.Write(&fDim);
    buf.Write(&fNState);
}

template<class TVar>
void TPZNullMaterial<TVar>::Read(TPZStream &buf, void *context) {
    TBase::Read(buf, context);
    buf.Read(&fDim);
    buf.Read(&fNState);
}


template<class TVar>
void TPZNullMaterial<TVar>::FillDataRequirements(TPZMaterialData& data) const{
    data.SetAllRequirements(false);
    data.fActiveApproxSpace = ENotActive;
}

template<class TVar>
void TPZNullMaterial<TVar>::GetSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col) const{
    u_len=du_row=du_col=0;
}

template<class TVar>
void TPZNullMaterial<TVar>::SetDimension(int dim) {
#ifdef PZDEBUG
    if (fDim < 1 || fDim > 3) { DebugStop(); }
#endif
    fDim = dim;
}

template<class TVar>
void TPZNullMaterial<TVar>::Print(std::ostream &out) const {
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMaterial::Print(out);
    out << "Dimension " << fDim << std::endl;
    out << "NState " << fNState << std::endl;
}

template<class TVar>
TPZMaterial *TPZNullMaterial<TVar>::NewMaterial() const {
    return new TPZNullMaterial<TVar>(*this);
}

template<class TVar>
void TPZNullMaterial<TVar>::Solution(const TPZMaterialDataT<TVar> &data, int var, TPZVec<TVar> &sol) {
    if (var == 0) {
        sol = data.sol[0];
    } else if (var == 1) {
        for (int i = 0; i < 3; i++) {
            sol[i] += data.dsol[0].GetVal(i, 0);
        }
    } else {
        DebugStop();
    }
}

template<class TVar>
int TPZNullMaterial<TVar>::NSolutionVariables(int var) const {
    if (var == 0) {
        return fNState;
    } else if (var == 1) {
        return 3;
    } else {
        DebugStop();
        return 0;
    }
}

template class TPZNullMaterial<STATE>;
template class TPZNullMaterial<CSTATE>;

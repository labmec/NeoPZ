#include "TPZMatLoadCases.h"
#include "TPZStream.h"
#include "TPZMaterial.h"
#include "TPZBndCondT.h"
#include "pzfmatrix.h"
#include "Hash/TPZHash.h"

void TPZMatLoadCasesBase::SetNumLoadCases(int numloadcases){
    if(numloadcases <= 0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nnumload cases: "<<numloadcases;
        PZError<<" must be a positive number\n";
        PZError<<"Aborting...\n";
        DebugStop();
    }
    fNumLoadCases = numloadcases;
}

void TPZMatLoadCasesBase::SetPostProcessIndex(int index)
{
#ifdef PZDEBUG
    if (index < 0 || index >= fNumLoadCases)
        {
            DebugStop();
        }
#endif
    fPostProcIndex = index;
}

int TPZMatLoadCasesBase::ClassId() const
{
    return Hash("TPZMatLoadCases");
}

void TPZMatLoadCasesBase::Read(TPZStream& buf, void* context)
{
    buf.Read(&fNumLoadCases);
    buf.Read(&fPostProcIndex);
}

void TPZMatLoadCasesBase::Write(TPZStream& buf, int withclassid) const
{
    buf.Write(&fNumLoadCases);
    buf.Write(&fPostProcIndex);
}


// this method is your chance to verify if the material to which this
// BC interface applies is compatible with this boundary interface
// it is called in the method SetMaterial of class TPZBndCondBase
template<class TVar>
void TPZMatLoadCasesBC<TVar>::SetMaterialImpl(TPZMaterial *mat)
{
    fMatLoadCases = dynamic_cast<TPZMatLoadCases<TVar>*>(mat);
    if(!fMatLoadCases){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Invalid reference for creating BC.\nAborting...\n";
        DebugStop();
    }
    this->fNumLoadCases = fMatLoadCases->NumLoadCases();
}
template<class TVar>
void TPZMatLoadCasesBC<TVar>::SetBCRhsValVec(TPZVec<TPZVec<TVar>>& bcValVec)
{
    const auto &nLoadCases = this->fNumLoadCases;
    const auto vecSize = bcValVec.size();
    if(nLoadCases != vecSize){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR:\nNum load cases = "<<nLoadCases<<'\n';
        PZError<<"\nVec size = "<<vecSize<<'\n';
        PZError<<"Aborting...\n";
        DebugStop();
    }
    fBCRhsValVec = bcValVec;
}
template<class TVar>
const TPZVec<TVar>& TPZMatLoadCasesBC<TVar>::GetBCRhsVal(int i) const
{
    const auto valVecSize = fBCRhsValVec.size();
    const auto &nLoadCases = this->fNumLoadCases;
#ifdef PZDEBUG
    if(i <0 || (valVecSize && i >= valVecSize)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR:\nNum load cases = "<<nLoadCases<<'\n';
        PZError<<"i = "<<i<<'\n';
        PZError<<"Aborting...\n";
        DebugStop();
    }
#endif
    if (valVecSize == 0){
        auto tmp =
            dynamic_cast<const TPZBndCondT<TVar>*>(this);
        if(!tmp){
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"\nERROR: Invalid reference for creating BC.\nAborting...\n";
            DebugStop();
        }
        return tmp->Val2();        
    }else{
        return fBCRhsValVec[i];
    }
}


template class TPZMatLoadCases<STATE>;
template class TPZMatLoadCases<CSTATE>;
template class TPZMatLoadCasesBC<STATE>;
template class TPZMatLoadCasesBC<CSTATE>;

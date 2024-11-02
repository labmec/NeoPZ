#include "TPZMatCombinedSpaces.h"
#include "TPZBndCondT.h"
#include "TPZMatBase.h"
#include "TPZMaterialDataT.h"

#include "Hash/TPZHash.h"


int TPZMatCombinedSpaces::ClassId() const{
    return Hash("TPZMatCombinedSpaces");
}

/** @brief Returns the variable index associated with a given name */

int TPZMatCombinedSpaces::VariableIndex(const std::string &name) const
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" should be implemented in your material";
    PZError<<" for any sort of post processing of the FEM solution\n";
    PZError<<"Aborting...\n";
    DebugStop();
}
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */

int TPZMatCombinedSpaces::NSolutionVariables(int var) const
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" should be implemented in your material";
    PZError<<" for any sort of post processing of the FEM solution\n";
    PZError<<"Aborting...\n";
    DebugStop();
}

template<class TVar>
void TPZMatCombinedSpacesT<TVar>::FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const
{
    for(auto &data : datavec) data.SetAllRequirements(false);
}
template<class TVar>
void TPZMatCombinedSpacesT<TVar>::FillDataRequirements(std::map<int, TPZMaterialDataT<TVar>> &datavec) const
{
    for(auto &data : datavec) data.second.SetAllRequirements(false);
}

template<class TVar>
void TPZMatCombinedSpacesT<TVar>::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<TVar>> &datavec) const
{
    for(auto &data : datavec) data.SetAllRequirements(false);
}

template<class TVar>
void TPZMatCombinedSpacesT<TVar>::Contribute(
    const TPZVec<TPZMaterialDataT<TVar>> &datavec,
    REAL weight,TPZFMatrix<TVar> &ef)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->Contribute(datavec, weight, fakeek, ef);
}

template<class TVar>
void TPZMatCombinedSpacesT<TVar>::ContributeBC(
    const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
    TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeBC(datavec, weight, fakeek, ef, bc);
}


template<class TVar>
int TPZMatCombinedSpacesT<TVar>::IntegrationRuleOrder(const TPZVec<int>& elPMaxOrder) const
{
    const int maxOrder = [&elPMaxOrder=std::as_const(elPMaxOrder)](){
        int max = 0;
        for (auto ord : elPMaxOrder)
            if (ord > max) max = ord;
        return max;
    }();
    auto *tmp = dynamic_cast<const TPZMaterialT<TVar>*>(this);
    auto ffporder = tmp->ForcingFunctionPOrder();
    //Check BC forcing function order
    auto *tmpBC = dynamic_cast<const TPZBndCondT<TVar>*>(this);
    if (tmpBC){
        ffporder = tmpBC->ForcingFunctionBCPOrder();
    }

    const int intOrder = maxOrder < ffporder ?
        maxOrder + ffporder : 2 * maxOrder;

    return  intOrder;
}

template<class TVar>
int TPZMatCombinedSpacesT<TVar>::ClassId() const{
    return Hash("TPZMatCombinedSpacesT") ^
        ClassIdOrHash<TVar>() << 1 ^
        TPZMatCombinedSpaces::ClassId() << 2;
}

// this method is your chance to verify if the material to which this
// BC interface applies is compatible with this boundary interface
// it is called in the method SetMaterial of class TPZBndCondBase
template<class TVar>
void TPZMatCombinedSpacesBC<TVar>::SetMaterialImpl(TPZMaterial *mat)
{
    fMatCombinedSpaces = dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(mat);
    if(!fMatCombinedSpaces){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Invalid reference for creating BC.\nAborting...\n";
        DebugStop();
    }
}

template<class TVar>
void TPZMatCombinedSpacesBC<TVar>::
Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
           TPZFMatrix<TVar> &ek,
           TPZFMatrix<TVar> &ef)
{
    auto *tmp = dynamic_cast<TPZBndCondT<TVar>*>(this);
    fMatCombinedSpaces->ContributeBC(datavec,weight,ek,ef,*tmp);
}

template<class TVar>
void TPZMatCombinedSpacesBC<TVar>::Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                TPZFMatrix<TVar> &ef)
{
    auto *tmp = dynamic_cast<TPZBndCondT<TVar>*>(this);
    fMatCombinedSpaces->ContributeBC(datavec,weight,ef,*tmp);
}

/**
 * @brief Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered
 * as not necessary.
 */
template<class TVar>
void TPZMatCombinedSpacesBC<TVar>::FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const
{
    auto *tmp = dynamic_cast<const TPZBndCondT<TVar>*>(this);
    fMatCombinedSpaces->FillBoundaryConditionDataRequirements(tmp->Type(), datavec);
}

template<class TVar>
void TPZMatCombinedSpacesBC<TVar>::ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                  TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                  TPZBndCondT<TVar> &bc)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}


template<class TVar>
void TPZMatCombinedSpacesBC<TVar>::Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec, int var,
              TPZVec<TVar> &sol)
{
    auto *tmp = dynamic_cast<TPZBndCondT<TVar>*>(this);
    fMatCombinedSpaces->Solution(datavec,var,sol);

}

template<class TVar>
int TPZMatCombinedSpacesBC<TVar>::ClassId() const
{
    return Hash("TPZMatCombinedSpacesBC");
}

template class TPZMatCombinedSpacesT<STATE>;
template class TPZMatCombinedSpacesT<CSTATE>;
template class TPZMatCombinedSpacesBC<STATE>;
template class TPZMatCombinedSpacesBC<CSTATE>;

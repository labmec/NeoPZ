#include "TPZMatSingleSpace.h"
#include "TPZMatBase.h"
#include "TPZBndCondT.h"
#include "TPZMaterialDataT.h"

#include "Hash/TPZHash.h"


int TPZMatSingleSpace::ClassId() const{
    return Hash("TPZMatSingleSpace");
}

void TPZMatSingleSpace::FillDataRequirements(TPZMaterialData &data) const{
    data.SetAllRequirements(false);
}

void TPZMatSingleSpace::FillBoundaryConditionDataRequirements(int type,
                                                                 TPZMaterialData &data) const{
    data.SetAllRequirements(false);
}

void TPZMatSingleSpace::GetSolDimensions(uint64_t &u_len,
                                         uint64_t &du_row,
                                         uint64_t &du_col) const{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nGetSolDimensions and Solution methods";
    PZError<<" should be implemented in your material";
    PZError<<" for any sort of post processing of the FEM solution\n";
    PZError<<"Aborting...\n";
    DebugStop();
}

template<class TVar>
void TPZMatSingleSpaceT<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                    REAL weight,TPZFMatrix<TVar> &ef)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->Contribute(data, weight, fakeek, ef);
}

template<class TVar>
void TPZMatSingleSpaceT<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                                      TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeBC(data, weight, fakeek, ef, bc);
}

template<class TVar>
void TPZMatSingleSpaceT<TVar>::Solution(const TPZMaterialDataT<TVar> &data, int var,
         TPZVec<TVar> &sol){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nGetSolDimensions and Solution methods";
    PZError<<" should be implemented in your material";
    PZError<<" for any sort of post processing of the FEM solution\n";
    PZError<<"Aborting...\n";
    DebugStop();
}

   
/** @brief Returns the variable index associated with a given name */

int TPZMatSingleSpace::VariableIndex(const std::string &name) const
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

int TPZMatSingleSpace::NSolutionVariables(int var) const
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" should be implemented in your material";
    PZError<<" for any sort of post processing of the FEM solution\n";
    PZError<<"Aborting...\n";
    DebugStop();
}


// template<class TVar>
// void TPZMatSingleSpaceT<TVar>::SolutionBC(const TPZMaterialDataT<TVar> &data, int var,
//          TPZVec<TVar> &sol){
//     PZError<<__PRETTY_FUNCTION__;
//     PZError<<" should be implemented in your material\n"
//            <<" for any sort of post processing of the FEM solution\n"
//            <<" over the boundaries\n"
//            <<"Aborting..."<<std::endl;
//     DebugStop();
// }

template<class TVar>
int TPZMatSingleSpaceT<TVar>::IntegrationRuleOrder(const int elPMaxOrder) const
{
    auto *tmp = dynamic_cast<const TPZMaterialT<TVar>*>(this);
    auto ffporder = tmp->ForcingFunctionPOrder();
    int integrationorder = 2*elPMaxOrder;
    //if the forcing function is not set, fForcingFunctionPOrder=0
    if (elPMaxOrder < ffporder) {
        integrationorder =
            ffporder+elPMaxOrder;
    }
    return  integrationorder;
}

template<class TVar>
int TPZMatSingleSpaceT<TVar>::ClassId() const{
    return Hash("TPZMatSingleSpaceT") ^
        ClassIdOrHash<TVar>() << 1 ^
        TPZMatSingleSpace::ClassId() << 2;
}

template<class TVar>
void TPZMatSingleSpaceBC<TVar>::Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                    TPZFMatrix<TVar> &ek,
                    TPZFMatrix<TVar> &ef)
{
    auto *tmp = dynamic_cast<TPZBndCondT<TVar>*>(this);
    fMatSingleSpace->ContributeBC(data,weight,ek,ef,*tmp);
}

template<class TVar>
void TPZMatSingleSpaceBC<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                           REAL weight,
                                           TPZFMatrix<TVar> &ef)
{
    auto *tmp = dynamic_cast<TPZBndCondT<TVar>*>(this);
    fMatSingleSpace->ContributeBC(data,weight,ef,*tmp);
}


template<class TVar>
void TPZMatSingleSpaceBC<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data,
                                             REAL weight,
                                             TPZFMatrix<TVar> &ek,
                                             TPZFMatrix<TVar> &ef,
                                             TPZBndCondT<TVar> &bc)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}



template<class TVar>
void TPZMatSingleSpaceBC<TVar>::FillDataRequirements(TPZMaterialData &data) const{
    auto *tmp = dynamic_cast<const TPZBndCondT<TVar>*>(this);
    fMatSingleSpace->FillBoundaryConditionDataRequirements(tmp->Type(), data);
}   

template<class TVar>
int TPZMatSingleSpaceBC<TVar>::ClassId() const
{
    return Hash("TPZMatSingleSpaceBC");
}

// this method is your chance to verify if the material to which this
// BC interface applies is compatible with this boundary interface
// it is called in the method SetMaterial of class TPZBndCondBase
template<class TVar>
void TPZMatSingleSpaceBC<TVar>::SetMaterialImpl(TPZMaterial *mat)
{
    fMatSingleSpace = dynamic_cast<TPZMatSingleSpaceT<TVar>*>(mat);

    if (!fMatSingleSpace) {
      PZError << __PRETTY_FUNCTION__;
      PZError << "\nERROR: Invalid reference for creating BC.\nAborting...\n";
      DebugStop();
    }
}

template class TPZMatSingleSpaceT<STATE>;
template class TPZMatSingleSpaceT<CSTATE>;
template class TPZMatSingleSpaceBC<STATE>;
template class TPZMatSingleSpaceBC<CSTATE>;

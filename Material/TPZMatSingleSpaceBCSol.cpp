#include "TPZMatSingleSpaceBCSol.h"
#include "TPZMaterial.h"
#include "Hash/TPZHash.h"

    template<class TVar>
    int TPZMatSingleSpaceBCSol<TVar>::ClassId() const
    {
        return Hash("TPZMatSingleSpaceBCSol");
    }

    // template<class TVar>
    // void TPZMatSingleSpaceBCSol<TVar>::Read(TPZStream&,void*)
    // {

    // }
    // template<class TVar>
    // void TPZMatSingleSpaceBCSol<TVar>::Write(TPZStream&,int) const
    // {

    // }

    template<class TVar>
    int TPZMatSingleSpaceBCSolBC<TVar>::ClassId() const
    {
        return Hash("TPZBndCondPostProcessInterface");
    }

    template<class TVar>
    void TPZMatSingleSpaceBCSolBC<TVar>::Read(TPZStream&,void*)
    {

    }

    template<class TVar>
    void TPZMatSingleSpaceBCSolBC<TVar>::Write(TPZStream&,int) const
    {

    }

  template<class TVar>
  void TPZMatSingleSpaceBCSolBC<TVar>::SetMaterialImpl(TPZMaterial* mat){
    fMatBndCondPPInterface = dynamic_cast<TPZMatSingleSpaceBCSol<TVar>*>(mat);
    TPZMatSingleSpaceBC<TVar>::SetMaterialImpl(mat);
  }

template class TPZMatSingleSpaceBCSol<STATE>;
template class TPZMatSingleSpaceBCSol<CSTATE>;



template class TPZMatSingleSpaceBCSolBC<STATE>;
template class TPZMatSingleSpaceBCSolBC<CSTATE>;

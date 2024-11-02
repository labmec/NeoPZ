#ifndef TPZBNDCONDBASE_H
#define TPZBNDCONDBASE_H
#include "TPZMaterialT.h"
#include "TPZBndCondT.h"


/**
 * @brief This class ensures that the boundary condition material is compatible
 * with any given interfaces.
 * @ingroup matsinglespace
 */
template<class TVar, class...Interfaces>
class TPZBndCondBase :
    public TPZMaterialT<TVar>,
    public TPZBndCondT<TVar>,
    public Interfaces... {
 protected:
    //! Default constructor
    TPZBndCondBase() = default;
    bool IsMatImpl() final{return true;}
 public:
    /*no need to specify destructor, copy/move assignment/destructors.
     the compiler can do it for us*/

    //! Constructor from an existing material setting `fVal1` and `fVal2`
    TPZBndCondBase(TPZMaterial * material,
                   int id,int type,
                   const TPZFMatrix<TVar> &val1,
                   const TPZVec<TVar> &val2);
    
    //! Prints data associated with the material.
    void Print(std::ostream &out = std::cout) const override;


    [[nodiscard]] int ClassId() const override;
    void Read(TPZStream &buf, void*context) override;
    void Write(TPZStream &buf, int withclassid) const override;

    
    void SetMaterial(TPZMaterial *) final;
    
    [[nodiscard]] int Dimension() const final
    {return this->fMaterial->Dimension() - 1;}
    
    [[nodiscard]] int NStateVariables() const final
    {return this->fMaterial->NStateVariables();}



    
    [[nodiscard]] int Id() const override {
        return TPZMaterial::Id();
    }

    [[nodiscard]] TPZMaterial* NewMaterial() const override {
        return new TPZBndCondBase(*this);
    }
        
    virtual void Clone(std::map<int, TPZMaterial * >&matvec) override {
        auto bccopy = new TPZBndCondBase(*this);
#ifdef PZDEBUG
        if(matvec.find(this->Id()) != matvec.end()) DebugStop();
#endif
        matvec[this->Id()] = bccopy;
        auto newmat = matvec[this->Id()];
        TPZBndCondBase *newbnd = dynamic_cast<TPZBndCondBase *>(newmat);
#ifdef PZDEBUG
        if(!newbnd) DebugStop();
#endif
        TPZMaterial *prev = this->Material();
        int previd = prev->Id();
        auto newmat_it = matvec.find(previd);
        if(newmat_it == matvec.end()) DebugStop();
        newbnd->SetMaterial(newmat_it->second);
    }
};

template<class TVar, class...Interfaces>
TPZBndCondBase<TVar, Interfaces...>::TPZBndCondBase(TPZMaterial * material,
                                                    int id,int type,
                                                    const TPZFMatrix<TVar> &val1,
                                                    const TPZVec<TVar> &val2) :
    TPZMaterialT<TVar>(id), TPZBndCondT<TVar>(type,val1,val2)
{
    SetMaterial(material);
}

template<class TVar, class...Interfaces>
void TPZBndCondBase<TVar, Interfaces...>::SetMaterial(TPZMaterial *mat){
    auto *tmp =
        dynamic_cast<TPZMaterialT<TVar>*>(mat);
    if(!tmp){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nIncompatible material type. Aborting...\n";
        DebugStop();
    }
    this->fMaterial = tmp;
    // the next statement will call the method SetMaterialImpl(mat) for all interfaces
    // it is called a fold expression (from c++17 on)
    // https://en.cppreference.com/w/cpp/language/fold
    (Interfaces::SetMaterialImpl(mat),...);
}



template<class TVar, class...Interfaces>
int TPZBndCondBase<TVar, Interfaces...>::ClassId() const{
    constexpr int nInterfaces = sizeof...(Interfaces);
    int id = Hash("TPZBndCondBase") ^
        TPZMaterialT<TVar>::ClassId() << 1 ^
        TPZBndCondT<TVar>::ClassId() << 2;
    
//    if constexpr (nInterfaces>0){
//        TPZVec<int> classIds = {Interfaces::ClassId()...};
//        for(auto count =0;count< nInterfaces;count++){
//            id = id ^ (classIds[count] << (count+3));
//        }
//    }
    return id;
}

//! Prints data associated with the material.
template<class TVar, class...Interfaces>
void TPZBndCondBase<TVar, Interfaces...>::Print(std::ostream &out) const {
    TPZMaterialT<TVar>::Print(out);
    TPZBndCondT<TVar>::Print(out);
    /* The following will perform calls to all the Read methods of the interfaces.
       This is a c++17 addition called fold expressions.*/
//    (Interfaces::Print(out),...);
}



template<class TVar, class...Interfaces>
void TPZBndCondBase<TVar, Interfaces...>::Read(
    TPZStream& buf, void* context)
{
    TPZMaterialT<TVar>::Read(buf,context);
    TPZBndCondT<TVar>::Read(buf,context);
    /* The following will perform calls to all the Read methods of the interfaces. 
       This is a c++17 addition called fold expressions.*/
//    (Interfaces::Read(buf,context),...);
}

template<class TVar, class...Interfaces>
void TPZBndCondBase<TVar, Interfaces...>::Write(
    TPZStream& buf, int withclassid) const
{
    TPZMaterialT<TVar>::Write(buf,withclassid);
    TPZBndCondT<TVar>::Write(buf,withclassid);
    /* The following will perform calls to all the Read methods of the interfaces. 
       This is a c++17 addition called fold expressions.*/
//    (Interfaces::Write(buf,withclassid),...);
}
#endif

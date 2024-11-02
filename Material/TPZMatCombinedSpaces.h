/**
 * @file TPZMatCombinedSpaces.h
 * @brief Header file for abstract class TPZMatCombinedSpacesT.\n
 * It implements interface for implementing the weak statement of the differential equation within the PZ environment for multiphysics materials.
 */

#ifndef TPZMATERIALCOMBINEDSPACES_H
#define TPZMATERIALCOMBINEDSPACES_H

#include "TPZSavable.h"
#include "TPZMatTypes.h"

class TPZMaterial;
class TPZBndCond;

template<class T>
class TPZMaterialDataT;
template<class T>
class TPZFMatrix;
template<class T>
class TPZBndCondT;
template<class T>
class TPZVec;



/**
 * @ingroup matcombinedspaces
 * @brief This class is a fundamental part of the TPZMaterial group.
 * It contains the type-agnostic interfaces that multiphysics  materials
 * should implement.
 */
class TPZMatCombinedSpaces : public virtual TPZSavable{
public:
    //! Default constructor
    TPZMatCombinedSpaces() = default;
    /** @brief Gets the order of the integration rule necessary to
        integrate an element with polynomial order `p` taking the 
        forcing function order into account.
        @param[in] elPMaxOrder maximum element order for each space
        @return adjusted integration rule order*/
    [[nodiscard]]virtual int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const = 0;
    
    [[nodiscard]] int ClassId() const override;

        /** @name PostProcess
     * @{
     */
    
    /** @brief Returns the variable index associated with a given name */
    [[nodiscard]] virtual int VariableIndex(const std::string &name) const;
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    [[nodiscard]] virtual int NSolutionVariables(int var) const;

    /**
     * @}
     */

};

template<class TVar>
class TPZMatCombinedSpacesBC;

/**
 * @ingroup matcombinedspaces
 * @brief This class is a fundamental part of the TPZMaterial group.
 * It contains the interfaces that multiphysics  materials
 * should implement.
 * It is noteworthy to observe that this definition does not depend
 * on the definition of the approximation spaces.
 * @note Materials should derive from TPZMatBase class.
 * @tparam TVar State variable type (`CSTATE`/`STATE`)
 */
template<class TVar>
class TPZMatCombinedSpacesT : public TPZMatCombinedSpaces{
public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatCombinedSpacesBC<TVar>;
    //! Default constructor
    TPZMatCombinedSpacesT() = default;
    /*no need to specify destructor, copy/move assignment/destructors.
     the compiler can do it for us*/

    /** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered
     * as not necessary.
     */
	virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const;
    /**
     * @brief Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered 
     * as not necessary.
     */
    virtual void FillDataRequirements(std::map<int, TPZMaterialDataT<TVar>> &datavec) const;

    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<TVar> > &datavec) const;

    /** @name Contribute */
    /** @{ */
    /**
     * @brief It computes a contribution to the stiffness matrix 
     * and load vector at one integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                            REAL weight,TPZFMatrix<TVar> &ek,
                            TPZFMatrix<TVar> &ef) = 0;
    /**
     * @brief It computes a contribution to the residual vector
     * at one integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ef is the residual vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                            REAL weight,TPZFMatrix<TVar> &ef);
    /** @name ContributeBC
        @ingroup Contribute*/
    /**@{*/
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one BC integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                              REAL weight, TPZFMatrix<TVar> &ek,
                              TPZFMatrix<TVar> &ef,
                              TPZBndCondT<TVar> &bc) = 0;
    /**
     * @brief It computes a contribution to the stiffness matrix 
     * and load vector at one BC integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                              REAL weight, TPZFMatrix<TVar> &ef,
                              TPZBndCondT<TVar> &bc);
    /**@}*/
    /**@}*/
    /** @brief Returns the solution associated with a given index
        based on the finite element approximation.
        @param[in] datavec Stores all the input data.
        @param[in] var Index of the queried solution
        @param[out] sol FEM Solution at the integration point
    */
    virtual void Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                          int var, TPZVec<TVar> &sol) = 0;

    [[nodiscard]] int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;
    
    [[nodiscard]] int ClassId() const override;

};

class TPZMaterial;

/**
   @brief Boundary condition interface for TPZMatCombinedSpacesT.
*/
template<class TVar>
class TPZMatCombinedSpacesBC : public TPZMatCombinedSpacesT<TVar>{
public:
    //! Default constructor
    TPZMatCombinedSpacesBC() = default;
    /** @brief Forward the call to TPZMatCombinedSpacesT::ContributeBC.*/
    void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                    TPZFMatrix<TVar> &ek,
                    TPZFMatrix<TVar> &ef) override;
    /** @brief Forward the call to TPZMatCombinedSpacesT::ContributeBC.*/
    void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                    TPZFMatrix<TVar> &ef) override;
    /** @brief This method should never be called. Throws.*/
    void ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    
    /** @brief This method passes the call to the associated material.*/
    virtual void Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec, int var,
                          TPZVec<TVar> &sol) override;
    /**
     * @brief Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as not necessary.
     */
    virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const override;

    [[nodiscard]] int ClassId() const override;
protected:
    /** @brief Pointer to material which created this BC. */
    TPZMatCombinedSpacesT<TVar>* fMatCombinedSpaces{nullptr};
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial *mat);    
};


extern template class TPZMatCombinedSpacesT<STATE>;
extern template class TPZMatCombinedSpacesT<CSTATE>;
extern template class TPZMatCombinedSpacesBC<STATE>;
extern template class TPZMatCombinedSpacesBC<CSTATE>;
#endif

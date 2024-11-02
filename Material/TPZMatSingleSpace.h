/**
 * @file TPZMatSingleSpace.h
 * @brief Header file for abstract class TPZMatSingleSpaceT.\n
 * It implements interface for implementing the weak statement of the differential equation within the PZ environment.
 */

#ifndef TPZMATERIALSINGLESPACE_H
#define TPZMATERIALSINGLESPACE_H

#include "TPZSavable.h"
#include "TPZMatTypes.h"

class TPZBndCond;
class TPZMaterialData;
template<class T>
class TPZMaterialDataT;
template<class T>
class TPZFMatrix;
template<class T>
class TPZBndCondT;
template<class T>
class TPZVec;


/**
 * @ingroup matsinglespace
 * @brief This class is a fundamental part of the TPZMaterial group.
 * It contains the type-agonstic interfaces that simple (non-multiphysics) 
 * materials should implement.
 */
class TPZMatSingleSpace : public virtual TPZSavable {
public:
    //! Default constructor
    TPZMatSingleSpace() = default;
    
    /** 
	 * \brief Fill the TPZMaterialData parameter.
     * The needed parameters for the calculations of the 
     * Contribute methods should be set.By default, 
     * all requirements are considered to be false.
     * \param The TPZMaterialData instance to be filled
     */
    virtual void FillDataRequirements(TPZMaterialData &data) const;

    /**
     * \brief Fill the TPZMaterialData parameter for boundary condition.
     * The needed parameters for the calculations at the boundaries
     * should be set. By default, all requirements are considered to be false.
     * \param type Boundary condition type.
     * \param data The TPZMaterialData instance to be filled
     */
    virtual void FillBoundaryConditionDataRequirements([[maybe_unused]]int type,
                                                      [[maybe_unused]]TPZMaterialData &data) const;
    
    /** @brief Gets the order of the integration rule necessary to
        integrate an element with polynomial order `p` taking the 
        forcing function order into account.
        @param[in] elPMaxOrder maximum element order
        @return adjusted integration rule order*/
    [[nodiscard]] virtual int IntegrationRuleOrder(const int elPMaxOrder) const = 0;

    /** 
       @brief Get the dimensions of the solution for each state variable.
       This will be used for initializing the corresponding TPZMaterialData.
       @param [out] u_len solution dimension.
       @param [out] du_row number of rows of the derivative
       @param [out] du_col number of columns of the derivative
    */
    virtual void GetSolDimensions(uint64_t &u_len,
                                  uint64_t &du_row,
                                  uint64_t &du_col) const;
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
class TPZMatSingleSpaceBC;
/**
 * @ingroup matsinglespace
 * @brief This class is a fundamental part of the TPZMaterial group.
 * It contains the interfaces that simple (non-multiphysics) materials
 * should implement.
 * It is noteworthy to observe that this definition does not depend
 * on the definition of the approximation space.
 * @note Materials should derive from TPZMatBase class.
 * @tparam TVar State variable type (`CSTATE`/`STATE`)
 */
template<class TVar>
class TPZMatSingleSpaceT : public TPZMatSingleSpace{
    friend class TPZBndCondT<TVar>;
public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC=TPZMatSingleSpaceBC<TVar>;
    //! Default constructor
    TPZMatSingleSpaceT() = default;
    /*no need to specify destructor, copy/move assignment/destructors.
     the compiler can do it for us*/
    /** @name Contribute*/
    /** @{ */
    /**
     * @brief It computes a contribution to the stiffness matrix 
     * and load vector at one integration point.
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    virtual void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                            TPZFMatrix<TVar> &ek,
                            TPZFMatrix<TVar> &ef) = 0;
    /**
     * @brief It computes a contribution to the residual vector
     * at one integration point.
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ef is the residual vector
     */
    virtual void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                            TPZFMatrix<TVar> &ef);
    /** @name ContributeBC
     @ingroup Contribute*/

    /**@{*/
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one BC integration point.
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                              TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                              TPZBndCondT<TVar> &bc) = 0;
    /**
     * @brief It computes a contribution to the stiffness matrix 
     * and load vector at one BC integration point.
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                              TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc);
    /**@}*/
    /**@}*/

    /** @brief Returns the solution associated with a given index
        based on the finite element approximation.
        This method should be implemented if any computations 
        on the solution are to be done.
        @param[in] data Stores all the input data.
        @param[in] var Index of the queried solution
        @param[out] sol FEM Solution at the integration point
    */
    virtual void Solution(const TPZMaterialDataT<TVar> &data, int var,
                          TPZVec<TVar> &sol);

    [[nodiscard]] int IntegrationRuleOrder(const int elPMaxOrder) const override;
    
    [[nodiscard]] int ClassId() const override;
};


class TPZMaterial;
/**
   @brief Boundary condition interface for TPZMatSingleSpaceT.
*/
template <class TVar>
class TPZMatSingleSpaceBC : public TPZMatSingleSpaceT<TVar>
{
public:
    //!Default constructor
    TPZMatSingleSpaceBC() = default;
    /** @brief Forward the call to TPZMatSingleSpaceT::ContributeBC.*/
    void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                    TPZFMatrix<TVar> &ek,
                    TPZFMatrix<TVar> &ef) override;
    /** @brief Forward the call to TPZMatSingleSpaceT::ContributeBC.*/
    void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                    TPZFMatrix<TVar> &ef) override;
    /** @brief This method should never be called. Throws.*/
    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief Forward the call to TPZMatSingleSpaceT::FillBoundaryConditionDataRequirements*/
    void FillDataRequirements(TPZMaterialData &data) const override;
    
    int ClassId() const override;
protected:
    /** @brief Pointer to material which created this BC. */
	TPZMatSingleSpaceT<TVar> * fMatSingleSpace{nullptr};
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial *mat);
};

extern template class TPZMatSingleSpaceT<STATE>;
extern template class TPZMatSingleSpaceT<CSTATE>;
extern template class TPZMatSingleSpaceBC<STATE>;
extern template class TPZMatSingleSpaceBC<CSTATE>;
#endif

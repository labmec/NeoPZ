/**
 * @file TPZBCPostProcessing.h
 * @brief Header file for defining a boundary condition post processing interface.\n
 * It defines the interface for post processing data associated with a boundary condition
 */

#ifndef TPZBCPOSTPROCESSINGINTERFACE_H
#define TPZBCPOSTPROCESSINGINTERFACE_H

#include "TPZSavable.h"
#include "TPZMatTypes.h"
#include "TPZMatSingleSpace.h"

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
 * @ingroup matsinglespaces
 * @brief This class is a fundamental part of the TPZMaterial group.
 * It contains the type-agnostic interfaces that singe space materials
 * should implement.
 */
template<class TVar>
class TPZMatSingleSpaceBCSolBC;

template<class TVar>
class TPZMatSingleSpaceBCSol : public TPZMatSingleSpaceT<TVar>
{
public:
    using TInterfaceBC = TPZMatSingleSpaceBCSolBC<TVar>;

public:
//declare all the virtual methods and etc as usual
    //@{
    //! Read and Write methods
    virtual int ClassId() const override;
    //@}
      /** @brief Returns the variable index associated with a given name (for boundary)*/
    [[nodiscard]] virtual int VariableIndexBC(const std::string &name) const = 0;
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var for a boundary material. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    [[nodiscard]] virtual int NSolutionVariablesBC(int var) const = 0;

    /** @brief Returns the solution associated with a given index
        on boundary elements.
        This method should be implemented if any computations 
        on the solution are to be done.
        @param[in] data Stores all the input data.
        @param[in] var Index of the queried solution
        @param[out] sol FEM Solution at the integration point
    */
    virtual void SolutionBC(const TPZMaterialDataT<TVar> &data, int var,
                          TPZVec<TVar> &sol) = 0;    

};

template<class TVar>
class TPZMatSingleSpaceBCSolBC : public TPZMatSingleSpaceBC<TVar>
{

public:

//declare all the virtual methods and etc as usual
    //@{
    //! Read and Write methods
    virtual int ClassId() const override;
	/** @brief Writes this object to the TPZStream buffer. Include the classid if `withclassid = true` */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads an objects from the TPZStream buffer. */
	virtual void Read(TPZStream &buf, void *context) override;

    /** @brief This method should never be called. Throws.*/
    void Solution(const TPZMaterialDataT<TVar> &data, int var,
                  TPZVec<TVar> &sol) override
    {
        fMatBndCondPPInterface->SolutionBC(data,var,sol);
    }
	/** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override
    {
        return fMatBndCondPPInterface->VariableIndexBC(name);
    }
	int NSolutionVariables(int var) const override
    {
        return fMatBndCondPPInterface->NSolutionVariablesBC(var);
    }
  

protected:
  TPZMatSingleSpaceBCSol<TVar> *fMatBndCondPPInterface{nullptr};
  void SetMaterialImpl(TPZMaterial* mat);

};
#endif
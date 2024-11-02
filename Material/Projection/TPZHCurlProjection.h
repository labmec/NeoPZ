/**
 * @file TPZHCurlProjection.h
 * @brief Header file for TPZHCurlProjection.
 * It implements a projection of a given solution using the HCurl norm
 */

#ifndef TPZHCURLPROJECTION_H
#define TPZHCURLPROJECTION_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpaceBCSol.h"
#include "TPZMatErrorSingleSpace.h"

/**
 * @brief Implements a projection of a given solution on a HCurl approximation space.
 * This material uses a HCurl approximation space in a geometric domain of
 * dimension defined by the user (2 or 3D) and projects a given solution using
 * the HCurl norm.
 * The solution to be projected is set by the Forcing Function of the 
 * TPZMaterialT.
 * It is expected to be a TPZVec<TVar> of size 3+ dim curl, as the solution is always in
 * a 3D space and the last positions are the curl.
 * For the error analysis, the solution is set through the ExactSol of
 * TPZMatErrorSingleSpace.
 */
template<class TVar=STATE>
class  TPZHCurlProjection :
    public TPZMatBase<TVar,TPZMatSingleSpaceBCSol<TVar>,
                      TPZMatErrorSingleSpace<TVar>  >
{
	using TBase = TPZMatBase<TVar,TPZMatSingleSpaceBCSol<TVar>,TPZMatErrorSingleSpace<TVar>>;
public:
    //! Default constructor
    TPZHCurlProjection();
    /** @brief Creates the TPZHCurlProjection material.
     * @param id material identifier
     * @param dim dimension
     */
    TPZHCurlProjection(int id, int dim);
    
    std::string Name() const override { return "TPZHCurlProjection"; }

    /** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 , ECurl = 2};
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {fScale = scale;}
    
    [[nodiscard]] REAL ScaleFactor() const
    {return fScale;}

	inline int Dimension() const  override { return this->fDim; }
    
    /** @brief Sets problem dimension */
    inline void SetDimension(int dim) {
        this->fDim = dim;
        if(this->fDim==1) fCurlDim = 1;
        else fCurlDim = 2*dim - 3;//1 for 2D and 2 for 3D
    }
	
    int NStateVariables() const override { return 1; }
    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief To create another material of the same type */
	TPZMaterial * NewMaterial() const override;
	
	/** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override;

  int VariableIndexBC(const std::string &name) const override
  {return VariableIndex(name);}
	
	int NSolutionVariables(int var) const override;
  
  int NSolutionVariablesBC(int var) const override;

  void Solution(const TPZMaterialDataT<TVar> &data,
                int var, TPZVec<TVar> &solOut) override;
  void SolutionBC(const TPZMaterialDataT<TVar> &data,
                int var, TPZVec<TVar> &solOut) override;
  void GetSolDimensions(uint64_t &u_len,
                        uint64_t &du_row,
                        uint64_t &du_col) const override;

    void Errors(const TPZMaterialDataT<TVar> &data,
                TPZVec<REAL> &errors) override;
    
    [[nodiscard]] int IntegrationRuleOrder(const int elPMaxOrder) const override;
    virtual int ClassId() const override;
protected:
    /** @brief Problem dimension */
	int fDim{-1};
    int fCurlDim{-1};
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
};

extern template class TPZHCurlProjection<STATE>;
extern template class TPZHCurlProjection<CSTATE>;
#endif


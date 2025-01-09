/**
 * @file TPZL2ProjectionCS.h
 * @brief Contains the TPZL2ProjectionCS class which implements an L2 projection of a given scalar solution
 */

#ifndef TPZL2PROJECTIONCS_H
#define TPZL2PROJECTIONCS_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "Projection/TPZL2Projection.h"
/**
 * @file TPZL2ProjectionCS.h
 * @brief Contains the TPZL2ProjectionCS class 
 * which implements an L2 projection of a given solution on a scalar approximation space 
 * to be used in FEM formulations with combined spaces.
 */

template<class TVar=STATE>
class TPZL2ProjectionCS :
    public TPZMatBase<TVar,TPZMatCombinedSpacesT<TVar>,TPZMatErrorCombinedSpaces<TVar>>{
	using TBase = TPZMatBase<TVar,TPZMatCombinedSpacesT<TVar>,TPZMatErrorCombinedSpaces<TVar>>;	
public:
    //! Default constructor
    TPZL2ProjectionCS() = default;
    /**
	 * @brief Class constructor 
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 */
	TPZL2ProjectionCS(int id, int dim, int nstate=1);
	/**
	 * @brief Class constructor setting a default constant solution
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 * @param sol constant solution. Ignored if a forcing function is set.
	 */
	TPZL2ProjectionCS(int id, int dim, int nstate, const TPZVec<TVar> &sol);
        
        TPZL2ProjectionCS(const TPZL2ProjectionCS &copy) = default;
        
        TPZL2ProjectionCS &operator=(const TPZL2ProjectionCS &copy) = default;

    std::string Name() const override { return "TPZL2ProjectionCS"; }
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 , EDerivative = 2};
	
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {fScale = scale;}

    void SetSol(const TPZVec<TVar> &sol)
    {
        fSol = sol;
        if(fSol.size()!=fNStateVars){
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"\nn state variables: "<<fNStateVars;
            PZError<<"\nn solutions: "<<fSol.size();
            PZError<<"\nAborting...\n";
            DebugStop();
        }
    }
    
    [[nodiscard]] REAL ScaleFactor() const
    {return fScale;}

        /// modifying the behavior of needing the solution
        void SetNeedSol(int need) {
            fNeedsSol = need;
        }
        
        
	int Dimension() const  override { return this->fDim; }
    
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }
	
    int NStateVariables() const override { return this->fNStateVars; }
    void Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief To create another material of the same type */
	TPZMaterial * NewMaterial() const override;
	
	/** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override;
	
	int NSolutionVariables(int var) const override;

    void Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                  int var, TPZVec<TVar> &solOut) override;
    
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const;

    void Errors(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                TPZVec<REAL> &errors) override;
    
    virtual int ClassId() const override;
        
        /**
         * @brief Fill material data parameter with necessary requirements for the
         * Contribute method. Here, in base class, all requirements are considered
         * as not necessary.
         */
        virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const override ;

protected:
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Number of state variables */
	int fNStateVars{1};
	
	/** @brief Constant solution vector. Ignored if forcing function is set. */
	TPZVec<TVar> fSol;
    
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
        
    /// Variable indicating whether the material will compute the residual or not
    int fNeedsSol = 0;
};


extern template class TPZL2ProjectionCS<STATE>;
extern template class TPZL2ProjectionCS<CSTATE>;

#endif

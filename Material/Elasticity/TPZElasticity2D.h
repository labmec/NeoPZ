/**
 * @file
 * @brief Contains the TPZElasticity2D class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef TPZELASTICITY2D_H
#define TPZELASTICITY2D_H


#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZMatLoadCases.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZElasticity2D : public TPZMatBase<STATE,
                                          TPZMatSingleSpaceT<STATE>,
                                          TPZMatErrorSingleSpace<STATE>,
                                          TPZMatLoadCases<STATE>>
{
public:
    using TBase = TPZMatBase<STATE,
                             TPZMatSingleSpaceT<STATE>,
                             TPZMatErrorSingleSpace<STATE>,
                             TPZMatLoadCases<STATE>>;
public :
    using ElasticityFunctionType = std::function<void (const TPZVec<REAL>&,
                                                       TPZVec<STATE> &,
                                                       TPZFMatrix<STATE>&)>;
	/** 
	 * @brief Creates an elastic material with:
	 * @param id material id
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
	 */
	TPZElasticity2D(int id, STATE E, STATE nu, STATE fx, STATE fy, int planestress = 1);
    
    TPZElasticity2D(int id);

    TPZElasticity2D(const TPZElasticity2D &copy) : TBase(copy), fE_def(copy.fE_def), fnu_def(copy.fnu_def),
    fElasticity(copy.fElasticity), ff(copy.ff), fEover21PlusNu_def(copy.fEover21PlusNu_def),
    fEover1MinNu2_def(copy.fEover1MinNu2_def), fPreStressXX(copy.fPreStressXX),
    fPreStressYY(copy.fPreStressYY), fPreStressXY(copy.fPreStressXY), fPreStressZZ(copy.fPreStressZZ),
    fPlaneStress(copy.fPlaneStress)
    {
        
    }

    /** @name Elasticity */
    /** @brief Set elasticity parameters */
    void SetElasticity(STATE E, STATE nu)
    {
        fE_def	= E;  // Young modulus
        fnu_def	= nu;   // poisson coefficient

    }
    
    /// Set a variable elasticity and poisson coefficient
    void SetElasticityFunction(ElasticityFunctionType func)
    {
        fElasticity = func;
    }
    
    /// Set the material configuration to plane strain
    void SetPlaneStrain()
    {
        fPlaneStress = 0;
    }
    
    /// Set the material configuration to plane stress
    void SetPlaneStress()
    {
        fPlaneStress = 1;
    }
    
    /** @brief Set forcing function */
    void SetBodyForce(STATE fx, STATE fy)
    {
        ff[0] = fx;
        ff[1] = fy;
        ff[2] = 0.;
    }
    
    STATE GetLambda(STATE E, STATE nu) const
    {
        STATE lambda = (nu*E)/((1.+nu)*(1.-2.*nu));
        return lambda;
    }
    
    STATE GetMU(STATE E, STATE nu) const
    {
        STATE mu = E/(2.*(1.+nu));
        return mu;
    }
/** @brief Returns the elasticity modulus E */
	STATE E() {return fE_def;}
	
	/** @brief Returns the poison coefficient modulus E */
	STATE Nu() {return fnu_def;}
	
	/** @brief Set PresStress Tensor */
	void SetPreStress(STATE Sigxx, STATE Sigyy, STATE Sigxy, STATE Sigzz);
       /**@}*/
    

	/** @brief Returns the model dimension */
	int Dimension() const  override { return 2;}
	
	/** @brief Returns the number of state variables associated with the material */
	 int NStateVariables() const override;
	
	/** @brief Print the material data*/
	void Print(std::ostream & out = std::cout) const override;
	
	/** @brief Returns the material name*/
	std::string Name()  const override { return "TPZElasticity2D"; }
		
	/** @name Contribute methods */
	/** @{ */
	
	/** @brief Calculates the element stiffness matrix */
	void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight,
                    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    /** @brief Applies the element boundary conditions */
	void ContributeBC(const TPZMaterialDataT<STATE> &data,STATE weight,
                      TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    void FillDataRequirements(TPZMaterialData &data) const override;

    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;
        
	
	/** @} */

    /** @name Solution */
	/** @{ */
    
	/** @brief Returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override;
	
	/** 
	 * @brief Returns the number of variables associated with the variable indexed by var.
	 */
	int NSolutionVariables(int var) const override;

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void Solution(const TPZMaterialDataT<STATE> &data, int var,
                  TPZVec<STATE> &Solout) override;
    
	/** @} */

    /** @name Errors */
	/** @{ */
    /** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
    int NEvalErrors() const override {return 6;}

    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override
    {
        u_len=2;
        du_row=2;
        du_col=2;
    };

    /** @} */
    
    /** @brief Creates a new material from the current object   ??*/
	TPZMaterial * NewMaterial()  const override;
    /**
       @name ReadWrite
       @{*/
    //!Read and Write methods
    int ClassId() const override;

	void Read(TPZStream &buf, void *context) override;
	
	void Write(TPZStream &buf, int withclassid) const override;
    /**@}*/
	/** 
	 * @brief Computes the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
protected:
    /** @brief Default constructor */
	TPZElasticity2D();
    
     
    void ContributeVecShape(const TPZMaterialDataT<STATE> &data, STATE weight,
                            TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    void ContributeVecShapeBC(const TPZMaterialDataT<STATE> &data, STATE weight,
                              TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                              TPZBndCondT<STATE> &bc);

    /** @name Errors */
	/** @{ */
	void Errors(const TPZMaterialDataT<STATE> &data,
                TPZVec<STATE> &values) override;
    /** @} */
    
	/** @brief Elasticity modulus */
	STATE fE_def;
	
	/** @brief Poison coeficient */
	STATE fnu_def;
	
    /** Elasticity function */
    ElasticityFunctionType fElasticity;
    
	/** @brief Forcing vector */
	TPZManVector<STATE,3> ff;
	
	/** @brief \f$ G = E/2(1-nu) \f$ */
	STATE fEover21PlusNu_def;
	
	/** @brief \f$ E/(1-nu) \f$ */
	STATE fEover1MinNu2_def;
	
	/** @brief Pre Stress Tensor - Sigma XX */
	STATE fPreStressXX;
	
	/** @brief Pre Stress Tensor - Sigma YY */
	STATE fPreStressYY;
	
	/** @brief Pre Stress Tensor - Sigma XY */
	STATE fPreStressXY;
    
    /** @brief Pre Stress Tensor - Sigma ZZ */
    STATE fPreStressZZ;
	
	/** @brief Uses plain stress */
	int fPlaneStress;
    
};

#endif

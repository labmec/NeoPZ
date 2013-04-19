//
//  pzporoelasticmf2d.h
//  PZ
//
//  Created by Agnaldo Farias on 7/20/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

/**
 * @file
 * @brief Contains the TPZPoroElasticMF2d class which implements the poroelastic 2d problem for simulation multi-phyisics,
 for the variables of displacement, flow and pressure.
 */

#ifndef PZ_pzporoelasticmf2d_h
#define PZ_pzporoelasticmf2d_h


#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzvec.h"

#include <iostream>


/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
/**
 **@ingroup equacao da elasticidade
 * \f$  div(T(u))  + fXf2 = 0  ==> Int{Grad(v).T(u)}dx - Int{v.gN}ds  = Int{ff.v}dx  \f$ (Eq. 1) 
 *
 *\f$ T(u) =  tr(E(u) − alpha*p*I)lambda*I + 2*nu*(E(u) − alpha*p*I)\f$
 *
 *@ingroup equacao da pressao
 * \f$ -1/visc*div(gradu k)  = 0 ==> k/visc*Int{Grad(u)Grad(v)}dx - Int{k/visv*Grad(u).n v}ds  = 0   (Eq. 2)  \f$ 
 *
 */

class TPZPoroElasticMF2d : public TPZDiscontinuousGalerkin{
	
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
	
	/** @brief Elasticity modulus */
	REAL fE;
	
	/** @brief Poison coeficient */
	REAL fnu;
	
	/** @brief constants poroelastic Biot*/
	REAL falpha; //parameter poroelastic Biot-Willis [dimensionless]
	REAL fSe; //or 1/M poroelastic storage coefficient at constant volume [dimensionless]
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Permeability of the rock and fluid viscosity*/
	REAL fk; 
	REAL fvisc;
	
	/** @brief Uses plain stress 
     * @note \f$fPlaneStress = 1\f$ => Plain stress state 
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state 
     */
	int fPlaneStress;
	
	/// timestep [s]
	REAL fTimeStep;
    
    REAL fTimeValue;
    
	int fmatId;
	
	/** @brief State: n ou n+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	static EState gState;
	
public:
	TPZPoroElasticMF2d();
	
	TPZPoroElasticMF2d(int matid, int dim);
	
	virtual ~TPZPoroElasticMF2d();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZPoroElasticMF2d"; }
	
	int Dimension() {return fDim;}
	
	virtual int NStateVariables();
	
	void SetLastState(){ gState = ELastState; }
	void SetCurrentState(){ gState = ECurrentState; }
    
	/** @brief Parameters of rock and fluid: */
	void SetParameters(REAL perm, REAL visc)
	{
		fk = perm;
		fvisc = visc;
	}
	
	/** 
	 * @brief Set parameters of elastic material:
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 */
	void SetParameters(REAL E, REAL nu,  REAL fx, REAL fy)
	{
		fE = E;
		fnu = nu;
		ff[0] = fx;
		ff[1] = fy;
	}
	
	/**
	 * @brief Set falpha parameter
	 * @param alpha : constant poroelastic Biot [dimensionless]
	 * @param Se : Coeficiente poroelastico de armazenamento a volume constante [adimensional]
	 */
	void SetBiotParameters(REAL alpha, REAL Se)
	{
		falpha = alpha;
		fSe = Se; 
	}
	
	/**
	 * @brief Set plane problem  
	 * planestress = 1 => Plain stress state 
	 * planestress != 1 => Plain Strain state 
	 */
	void SetfPlaneProblem(int planestress)
	{
		fPlaneStress = planestress;
	}
	
	/// Set the timestep
	void SetTimeStep(REAL delt)
	{
		fTimeStep = delt;
	}
    
    /// Set the timestep
	void SetTimeValue(REAL TimeValue)
	{
		fTimeValue = TimeValue;
	}	
    
    int MatId()
    {
        return fmatId;
    }
	
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /**
     * @brief Applies to Dirichlet boundary condition for the elasticity equation (mechanical problem)
     */
    void ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Neumann boundary condition for elasticity equation (mechanical problem)
     */
    void ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Mixed boundary condition for elasticity equation (mechanical problem)
     */
    void ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /**
     *@brief Dirichlet with free boundary in the y-direction in the equation of elasticity.  
     */
    void ApplyDirichletFreeY_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /**
     *@brief Neumann with free boundary in the y-direction in the equation of elasticity.  
     */
    void ApplyNeumannFreeX_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc);
    
    
    /*
     * @brief Applies to Dirichlet boundary condition for mixed problem (pressure and flux)
     * In the mixed formulation, the contribution of the Dirichlet boundary condition for the pressure appears in the flow equation and not in the equation of pressure
     */
    void ApplyDirichlet_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Neumann boundary condition for mixed problem (pressure and flux)
     */
    void ApplyNeumann_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Mixed boundary condition for mixed problem (pressure and flux)
     */
    void ApplyMixed_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
    
	//public:
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                       REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
	
	/** @name Contribute methods
	 * @{
	 */
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef) {
		DebugStop();
	}
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc){
		DebugStop();
	}
	
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
};

#endif

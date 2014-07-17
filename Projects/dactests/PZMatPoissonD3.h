//
//  PZMatPoissonD3.h
//  PZ
//
//  Created by Douglas Castro on 5/23/14.
//
//

#ifndef __PZ__PZMatPoissonD3__
#define __PZ__PZMatPoissonD3__

#include <iostream>
#include "pzdiscgal.h"
#include "pzmaterial.h"


/** Material para problema de Poisson 3D */
/**  Div.(-fK(x,y,z) Grad u) = fF    em \Omega
 *   u  = uD                       em \partial\Omega
 *   du/dn = g                     em \partial\Omega
 */

/**
 *   q = -fK(x,y,z) Grad u          em \Omega
 *   Div.(q) = fF                    em \Omega
 *   u  = uD                       em \partial\Omega
 *   du/dn = g                     em \partial\Omega
 */


class TPZMatPoissonD3 : public TPZDiscontinuousGalerkin {
    
protected:
	/** Material Id */
    int fMatId;
    
    /** Valor da funcao de carga */
    REAL fCarga; //fF
    
    /** Dimensao do dominio */
    int fDim;
    
    /** Coeficiente que multiplica o gradiente */
    REAL fK;
	
public:
    
	TPZMatPoissonD3();
	
	TPZMatPoissonD3(int matid, int dim);
    
	virtual ~TPZMatPoissonD3();
    
    /** @brief copy constructor */
    TPZMatPoissonD3(const TPZMatPoissonD3 &copy);
    
    TPZMatPoissonD3 &operator=(const TPZMatPoissonD3 &copy);
	
	virtual std::string Name() { return "TPZMatPoissonD3"; }
	
	int Dimension() const {return fDim;}
    
    int MatId()
    {
        return fMatId;
    }
	
	virtual int NStateVariables();
    
//    void SetForcesPressure(REAL fxfPQ);
//    
//    void SetForcesSaturation(REAL fxfS);
//    
//    void SetPermeability(REAL perm);
//    
//    void SetViscosity(REAL visc);
//    
//    void SetPorosity(REAL poros);
//    
//    void GetPermeability(REAL &perm);
//    
//	void SetConvectionDirection(TPZVec<REAL> convdir);
//	
//	void GetConvectionDirection(TPZVec<REAL> &convdir);
//    
//    void SetLastState(){
//        gState = ELastState;
//    }
//    
//	void SetCurrentState(){
//        gState = ECurrentState;
//    }
//    
//    void SetPressureEqFilter(){
//        fPressureEquationFilter = true;
//    }
//    void SetFalsePressureEqFilter(){
//        fPressureEquationFilter = false;
//    }
//    
//	void SetTimeStep(REAL delt){
//		fTimeStep = delt;
//	}
//    
//	void SetTimeValue(REAL TimeValue){
//		fTimeValue = TimeValue;
//	}
//    
//    void GetTimeValue(REAL &TimeValue){
//		TimeValue = fTimeValue;
//	}
    
    void Print(std::ostream &out);
    
	/** @name Contribute methods
	 * @{
	 */
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since June 2, 2014
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since June 2, 2014
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
		DebugStop();
	}
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since June 2, 2014
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since June 2, 2014
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
		DebugStop();
	}
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param datavec [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since June 2, 2014
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since June 2, 2014
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    void         ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param datavec [in]
	 * @param dataleft [in]
     * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 2, 2014
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 2, 2014
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {DebugStop();}
    
	
	//virtual int VariableIndex(const std::string &name);
	
	//virtual int NSolutionVariables(int var);
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
	//virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
    
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                                 TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &/*flux*/,
                TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
};



#endif /* defined(__PZ__PZMatPoissonD3__) */

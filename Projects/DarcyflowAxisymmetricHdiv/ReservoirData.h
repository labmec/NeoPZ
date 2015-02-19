#ifndef TPZProblemDATAH
#define TPZProblemDATAH
/*
 *  ReservoirData.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

class ReservoirData {
	
public:

    /** @brief State: Stiffness or Mass Matrix Calculations */
    enum EState { ELastState = 0, ECurrentState = 1 };
    
    /**
     * @ingroup Characteristic Parameters
     * @brief Define characteristic parameters for Darcy linear flow.
     * @since December 08, 2014
     */
    
    /** @brief Characteristic length - m */
    REAL fLref;
    
    /** @brief Characteristic Permeability - m2 */
    REAL fKref;
    
    /** @brief Characteristic Pressure - Pa */
    REAL fPref;
    
    /** @brief Characteristic Density - kg/m3 */
    REAL fRhoref;
    
    /** @brief Characteristic viscosity - Pa s */
    REAL fEtaref;
    
    /** @brief Density at P of reference - kg/m3 */
    REAL fRhoRef;
    
    /** @brief Porosity at P of reference - */
    REAL fPhiRef;
    
    /** @brief absolute permeability */
    TPZFMatrix<REAL> fKab;
	
	ReservoirData();
	
	~ReservoirData();

	/**
	 * @brief \f$ Rock porosity. \f$ Phi = Phi( P ) \f$
	 * @param P fluid pressure
	 */	
	void Porosity(REAL P, REAL &poros, REAL &dPorosDp);

	/** 
	 * @brief \f$ Oil density RhoOil = RhoOil( P ) \f$
	 * @param P fluid pressure
	 */
	void RhoOil(REAL P, REAL &Rho, REAL &dRhoDpo);

	/**
	 * @brief Oil viscosity. \f$ OilViscosity = ViscOil( P ) \f$
	 * @param P fluid pressure
	 */
	void OilViscosity(REAL P, REAL &Viscosity, REAL &dViscosityDpo);

    /** @brief Set the characteristic length - m */
    void SetLref(REAL Lref) {fLref = Lref; }
    
    /** @brief Characteristic length - m */
    REAL Lref() {return fLref; }

    /** @brief Set the characteristic Permeability - m2 */
    void SetKref(REAL Kref) {fKref = Kref;}
    
    /** @brief Characteristic Permeability - m2 */
    REAL Kref() {return fKref;}

    /** @brief Set the characteristic Pressure - Pa */
    void SetPref(REAL Pref) {fPref = Pref;}
    
    /** @brief Characteristic Pressure - Pa */
    REAL Pref() {return fPref;}

    /** @brief Set the characteristic Density - kg/m3 */
    void Rhoref(REAL Rhoref) {fRhoref = Rhoref;}
    
    /** @brief Characteristic Density - kg/m3 */
    REAL Rhoref() {return fRhoref;}

    /** @brief Set the characteristic viscosity - Pa s */
    void SetEtaref(REAL Etaref) {fEtaref = Etaref;}
    
    /** @brief Characteristic viscosity - Pa s */
    REAL Etaref() {return fEtaref;}

    /** @brief Porosity at P of reference - */
    void SetPhiRef(REAL PhiRef) {fPhiRef = PhiRef;}
    
    /** @brief Porosity at P of reference - */
    REAL PhiRef() {return fPhiRef;}
	
};


#endif
/*
 *  TPZMohrCoulomb.h
 *  FEMPZ
 *
 *  Created by Diogo Cecilio on 5/4/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


/* Generated by Together */// $Id: TPZMohrCoulomb.h,v 1.2 2010-06-11 22:12:14 diogo Exp $

#ifndef TPZMOHRCOULOMB_H
#define TPZMOHRCOULOMB_H

#include "pzlog.h"
#include "TPZPlasticStep.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzvec_extras.h"
#include "TPZPlasticStepID.h"

#ifdef LOG4CXX_PLASTICITY
static LoggerPtr loggerMohrCoulomb(Logger::getLogger("MCC"));
#endif

#define MOHRCOULOMBPARENT TPZPlasticStep<TPZYCMohrCoulomb, TPZThermoForceA, TPZElasticResponse>


class TPZMohrCoulomb : public MOHRCOULOMBPARENT  {
	
public:
	
	enum {NYield = TPZYCMohrCoulomb::NYield};
	
public:
	
    TPZMohrCoulomb():MOHRCOULOMBPARENT()
    {
		fMaterialTensionSign  = 1; // internally in this material tension is negative
		fInterfaceTensionSign =  1; // by default
    }
	
    TPZMohrCoulomb(const TPZMohrCoulomb & source):MOHRCOULOMBPARENT(source)
    {
    }
	
    TPZMohrCoulomb & operator=(const TPZMohrCoulomb & source)
    {
		MOHRCOULOMBPARENT::operator=(source);		
		return *this;
    }
	
	virtual const char * Name() const
	{
		return "TPZMohrCoulomb";	
	}
	
    static void ConventionalConcrete(TPZMohrCoulomb & material)
	{
    	REAL pi = M_PI;
	    REAL cohesion = 11.2033; //yield- coesao inicial correspondeno a fck igual 32 Mpa
	    REAL phi =  20./180. * pi; //phi=20
	    REAL hardening = 1000.; //Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao
	    REAL young = 20000.;
	    REAL poisson = 0.2;
        material.fYC.SetUp(phi);
		material.fTFA.SetUp(cohesion, hardening);
		material.fER.SetEngineeringData(young, poisson);
	}
    
    static void TaludeMaterial(TPZMohrCoulomb & material)
	{
    	REAL pi = M_PI;
        REAL cohesion = 50.; //yield- coesao inicialem KPa
	    REAL phi =  20./180. * pi; //phi=20
	    REAL hardening = 10.; //Modulo de hardening da coesao equivante 0.01 Mpa a cada 0.1% de deformacao
	    REAL young = 20000.;//E em KPa
	    REAL poisson = 0.49;
        material.fYC.SetUp(phi);
		material.fTFA.SetUp(cohesion, hardening);
		material.fER.SetEngineeringData(young, poisson);
	}
	
	void SetUp(REAL & cohesion, REAL & phi, REAL & hardening, REAL &young, REAL &poisson)
	{
		MOHRCOULOMBPARENT::fYC.SetUp(phi);
		MOHRCOULOMBPARENT::fTFA.SetUp(cohesion, hardening);
		MOHRCOULOMBPARENT::fER.SetEngineeringData(young, poisson);
	}
	
    virtual void SetUp(const TPZTensor<REAL> & epsTotal) {
        MOHRCOULOMBPARENT::SetUp(epsTotal);
    }

	virtual void Print(std::ostream & out) const
	{
		out << "\n" << this->Name();
		out << "\n Base Class Data:\n";
		MOHRCOULOMBPARENT::Print(out);		
	}
	
	public:
int ClassId() const override;


    void Write(TPZStream &buf, int withclassid) const override{
        MOHRCOULOMBPARENT::Write(buf, withclassid);

        buf.Write(&fYC.fPhi, 1);

        REAL lambda = fER.Lambda();
        REAL mu = fER.Mu();
        buf.Write(&lambda, 1);
        buf.Write(&mu, 1);

        buf.Write(&fTFA.fSigmaYield0, 1);
        buf.Write(&fTFA.fK, 1);

        buf.Write(&fResTol, 1);
        buf.Write(&fIntegrTol, 1);
        buf.Write(&fMaxNewton, 1);
        buf.Write(&fMinLambda, 1);

        buf.Write(&fN.m_eps_t.fData[0], 6);
        buf.Write(&fN.m_eps_p.fData[0], 6);
        buf.Write(&fN.m_hardening, 1);

        // fPlasticMem does not need to be stored
    }
    
    void Read(TPZStream& buf, void* context) override {
        MOHRCOULOMBPARENT::Read(buf, context);

        buf.Read(&fYC.fPhi, 1);

        REAL lambda = fER.Lambda();
        REAL mu = fER.Mu();
        buf.Read(&lambda, 1);
        buf.Read(&mu, 1);

        buf.Read(&fTFA.fSigmaYield0, 1);
        buf.Read(&fTFA.fK, 1);

        buf.Read(&fResTol, 1);
        buf.Read(&fIntegrTol, 1);
        buf.Read(&fMaxNewton, 1);
        buf.Read(&fMinLambda, 1);

        buf.Read(&fN.m_eps_t.fData[0], 6);
        buf.Read(&fN.m_eps_p.fData[0], 6);
        buf.Read(&fN.m_hardening, 1);

        fPlasticMem.Resize(0);
    }
    
		
public:
    
    virtual int GetNYield() const {
        return as_integer(NYield);
    }

};


#endif //TPZMohrCoulomb_H

#ifndef MATHYPERELASTICHPP
#define MATHYPERELASTICHPP

#include "pzmaterial.h"
#include "pzfmatrix.h"




class TPZMatHyperElastic : public TPZMaterial {

   REAL fK2[3][3],fK3[3][3],fK4[3][3],fK6[3][3],fK7[3][3],fK8[3][3],fXf[3];
   REAL fL1[3][3],fL2[3][3],fL3[3][3],fL4[3][3],fL5[3][3],fL6[3][3],fL7[3][3],fL8[3][3],fL9[3][3],fGradtrC[3][3];
	REAL fE1[3],fE5[3],fE9[3],fGradDetF[3][3];
	REAL fLambda,fNu,fE,fMu;
   REAL fCoef1,fCoef2,fCoef3;

public :

TPZMatHyperElastic(int nummat,REAL e,REAL mu,REAL nu=-1.,REAL lambda=-1.,REAL coef1=-1.,REAL coef2=-1.,REAL coef3=-1.);

virtual ~TPZMatHyperElastic();

void SetMaterial(TPZFMatrix &xfin){
   fXf[0] = xfin(0,0);
   fXf[1] = xfin(1,0);
   fXf[2] = xfin(2,0);
}

int Dimension() { return 3;}

int NStateVariables();

virtual void Print(ostream & out);

char *Name() { return "TPZMatHyperElastic"; }

virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv ,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);

virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,double weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

virtual int VariableIndex(char *name);

virtual int NSolutionVariables(int var);

virtual int NFluxes(){ return 9;}

virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

/**compute the value of the flux function to be used by ZZ error estimator*/
virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);

void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				  TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
		        TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric
};

#endif

//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOLINEARH
#define TPZGEOLINEARH


#include "pzvec.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;

class TPZGeoLinear {

public:
	enum {NNodes = 2, NSides = 3};

static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

static void Jacobian(TPZFMatrix nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

static TPZIntPoints *CreateSideIntegrationRule(int side, int order);

static int NSubElements();
};
#endif


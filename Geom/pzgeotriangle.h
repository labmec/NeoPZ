// TPZGeoQuad.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOTRIANGLEH
#define TPZGEOTRIANGLEH

#include "pzvec.h"

class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;

class TPZGeoTriangle  
{
public:

	enum {NNodes = 3, NSides = 7};

  /** implementation of two-dimensional bilinear interpolation*/
static  void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);

  /**Computes the jacobian*/
static  void Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

  /**Computes the geometric location*/
static  void X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);

  /**
   * Method which creates a geometric boundary condition 
   * element based on the current geometric element, 
   * a side and a boundary condition number
   */
static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);

  /**
   * Create an integration rule 
   * @param order order of the integration rule to be created
   * @param side side to create integration rule
   */
static TPZIntPoints * CreateSideIntegrationRule(int side, int order);


};

#endif 

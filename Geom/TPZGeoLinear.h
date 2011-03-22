// -*- c++ -*-
// $ Id: $
//HEADER FILE FOR CLASS TPZGeoCube

#ifndef TPZGEOLINEARH
#define TPZGEOLINEARH

#include "pznoderep.h"

#include "pzvec.h"
#include "pzeltype.h"
#include "pzfmatrix.h"
#include "tpzline.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {

/// implements the geometry of a one dimensional linear element
  class TPZGeoLinear : public TPZNodeRep<2, pztopology::TPZLine> {

public:
	enum {NNodes = 2};

  /**
  * Constructor with list of nodes
   */
 TPZGeoLinear(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes, pztopology::TPZLine>(nodeindexes)
 {
 }
  
  /**
  * Empty constructor
   */
 TPZGeoLinear() : TPZNodeRep<NNodes, pztopology::TPZLine>()
 {
 }
  
  /**
  * Constructor with node map
   */
 TPZGeoLinear(const TPZGeoLinear &cp,
                std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp,gl2lcNdMap)
 {
 }
  
  /**
  * Copy constructor
   */
 TPZGeoLinear(const TPZGeoLinear &cp) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
 {
 }

  /**
  * Copy constructor
   */
 TPZGeoLinear(const TPZGeoLinear &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
 {
 }


/**
 * returns the type name of the element
 */
static std::string TypeName() { return "Linear";} 

static void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);

/**
 * returns the projection of a given point from "NSide - 1" side to "side".
 */
static bool MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide);

static void Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi);

static void Jacobian(TPZFMatrix &nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,
				TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);

static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);

	  
  public:
	  /**
	   * Creates a geometric element according to the type of the father element
	   */
	  static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										TPZVec<int>& nodeindexes,
										int matid,
										int& index);
	  
	  


};

inline void TPZGeoLinear::Jacobian(TPZFMatrix &coord,TPZVec<REAL> &param,TPZFMatrix &jacobian,
							TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {


//   REAL mod1 = (coord(0,1)-coord(0,0))*0.5;
//   jacobian(0,0) = mod1;
//   detjac = mod1;
//   jacinv(0,0) = 1./detjac;
//   axes(0,0) = 1.;

  //VERSAO FUNCIONAL
  jacobian.Resize(1,1); axes.Resize(1,3); jacinv.Resize(1,1);
  int ic;
  REAL v1[3] = {0.};
  int nrow = coord.Rows();
  REAL mod1 = 0.;
  for(ic=0; ic<nrow; ic++) {
    v1[ic] = (coord(ic,1)-coord(ic,0))*0.5;
    mod1 += v1[ic]*v1[ic];
  }
  mod1 = sqrt(mod1);
  jacobian(0,0) = mod1;
  detjac = mod1;

  if(IsZero(detjac)){
    detjac = ZeroTolerance();
  }

  jacinv(0,0) = 1./detjac;

  for(ic=0; ic<3; ic++) {
    axes(0,ic) = v1[ic]/detjac;
  }

  // VERSAO ORIGINAL
//   TPZFNMatrix<9> phi(NNodes, pztopology::TPZLine,1);
//   TPZFNMatrix<18> dphi(1,NNodes, pztopology::TPZLine);
//   Shape(param,phi,dphi);

//   int ic;
//   TPZManVector<REAL,3> v1(3,0.);
//   REAL mod1 = 0.;

//   for(int i=0; i < NNodes, pztopology::TPZLine; i++) {
//     for(ic = 0; ic < 3; ic++) {
//       v1[ic] += coord(ic,i)*dphi(0,i);
//     }
//   }

//   for(ic=0; ic<3; ic++) {
//     mod1 += v1[ic]*v1[ic];
//   }
//   mod1 = sqrt(mod1);
//   jacobian(0,0) = mod1;
//   detjac = mod1;
//   jacinv(0,0) = 1./detjac;


  //  axes.Zero();
//   for(ic=0; ic<3; ic++) {
//     axes(0,ic) = v1[ic]/mod1;
//   }
}

inline void TPZGeoLinear::X(TPZFMatrix &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){

  int ic;
  REAL xi = loc[0];
  int nrow = coord.Rows();
  for(ic=0; ic<nrow; ic++) result[ic] = coord(ic,0)*(1.-xi)*0.5+coord(ic,1)*(1.+xi)*0.5;

//   TPZFNMatrix<9> phi(NNodes, pztopology::TPZLine,1);
//   TPZFNMatrix<18> dphi(1,NNodes, pztopology::TPZLine);
//   Shape(loc,phi,dphi);
//   int in,ic;
//   for(in=0; in<3; in++) result[in] = 0.;
//   for(in = 0; in < NNodes, pztopology::TPZLine; in++) {
//     for(ic=0; ic<3 ; ic++) {
//       result[ic] += coord(ic,in)*phi(in,0);
//     }
//   }

}

inline bool TPZGeoLinear::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide) {
     TPZTransform Transf = pztopology::TPZLine::SideToSideTransform(TPZGeoLinear::NSides - 1, side);
     SidePar.Resize(SideDimension(side));
     Transf.Apply(InternalPar,SidePar);

     JacToSide = Transf.Mult();
	   return true;
}

};
#endif


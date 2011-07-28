//$Id: pzgeoelbc.cpp,v 1.5 2006-10-17 01:38:03 phil Exp $

#include "pzgeoelbc.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzgmesh.h"

TPZGeoElBC::TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh):fCreatedElement(NULL){
	if (!el){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " - TPZGeoEl *el is NULL\n";
		return;
	}
	this->fCreatedElement = el->CreateBCGeoEl( side, matid );
}

TPZGeoElBC::TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh):fCreatedElement(NULL){
	TPZGeoEl * el = elside.Element();
	const int side = elside.Side();
	if (!el || side == -1){
		if (!el) PZError << "Error at " << __PRETTY_FUNCTION__ << " - TPZGeoEl *elside.Element() is NULL\n";
		if (side == -1) PZError << "Error at " << __PRETTY_FUNCTION__ << " - int elside.Side() is -1\n";
		return;
	}
	this->fCreatedElement = el->CreateBCGeoEl( side, matid );
}

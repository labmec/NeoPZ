/* Generated by Together */
#include "pzcartsys.h"

/**
 *  Default empty constructor
 */
TPZCartsys::TPZCartsys():TPZCosys(),fOrigin(3,0.){
  int i,j;
  fOrigin.Resize(3,0.);
  //	fOrigin->Redim(3);
  for (i=0 ; i<3 ; ++i)	fOrigin[i] = 0.;
  fReference = NULL;
  fTr.Redim(3,3);
  fTr.Zero();
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      fTr(i,j)=(i == j) ? 1.0 : 0.0;
}

/**Constructor - 
 * create one object from a reference coordinate system
 * and a origin and euller angles
 */
TPZCartsys::TPZCartsys(int num, TPZCartsys * ref, TPZVec<REAL> * org, TPZVec<REAL> *angles):TPZCosys(num,ref),fOrigin(3,0.){
  int i,j;
  if (org) {
    fOrigin = *org;
  }
  if (angles) CalcRotMatrix(*angles);
  else {
    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
	fTr(i,j) = (i == j) ? 1.0 : 0.0;
  }
}

/**Constructor - 
 * create one object from a reference coordinate system
 * and a origin and a rotation matrix
 */
TPZCartsys::TPZCartsys(int num, TPZCartsys * ref, TPZVec<REAL> * org, TPZFMatrix *angles):TPZCosys(num,ref),fOrigin(3,0.){
  int i,j;
  if (org) fOrigin = *org;
  if (angles) fTr = *angles;
  else 
    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
	fTr(i,j) = (i == j) ? 1.0 : 0.0;
}

/**Destructor - 
 *  
 */
TPZCartsys::~TPZCartsys(){
}

void TPZCartsys::GetAxes(TPZFMatrix &axes){
  axes=fTr;
}

void TPZCartsys::SetAxes(TPZVec<REAL> &x, TPZVec<REAL> &z){

  x[0] -= fOrigin[0]; x[1] -= fOrigin[1]; x[2] -= fOrigin[2];
  Normalise(x);
  fTr(0,0) = x[0]; fTr(0,1) = x[1]; fTr(0,2) = x[2];
  
  z[0] -= fOrigin[0]; z[1] -= fOrigin[1]; z[2] -= fOrigin[2];
  Normalise(z);
  fTr(2,0) = z[0]; fTr(2,1) = z[1]; fTr(2,2) = z[2];
  
  TPZVec<REAL> y(3,0.);
  GetNormal(z,x,y);
  fTr(1,0) = y[0]; fTr(1,1) = y[1]; fTr(1,2) = y[2];
}


void TPZCartsys::SetAxes(TPZFMatrix &RotMat){
  fTr = RotMat;
}

void TPZCartsys::SetOrigin(TPZVec<REAL> *org){
  fOrigin = *org;
}

void TPZCartsys::TransformGradient(TPZVec<REAL> &X, TPZFMatrix &GradX, TPZVec<REAL> &x, TPZFMatrix &Gradx, TPZCosys *dest)  {
  int i,j;
  Gradx.Zero();
  for (i=0;i<3;i++){
    x=X;
    for (j=0;j<3;j++)
      Gradx(i,i)+=fTr(i,j)*GradX(j,j);
  }
  if (dest != fReference) {
    TPZFMatrix gradin(GradX);
    FromReference(x);
    TransformGradient(x,gradin,x,Gradx,dest);
  }
}

void TPZCartsys::FromReference (TPZVec<REAL> &point) {

    if (!fReference) return;
    point[0] -= fOrigin[0];
    point[1] -= fOrigin[1];
    point[2] -= fOrigin[2];

    TPZVec<REAL> newp(3,0.);
    newp[0] = fTr(0,0)*point[0] + fTr(1,0)*point[1] + fTr(2,0)*point[2];
    newp[1] = fTr(0,1)*point[0] + fTr(1,1)*point[1] + fTr(2,1)*point[2];
    newp[2] = fTr(0,2)*point[0] + fTr(1,2)*point[1] + fTr(2,2)*point[2];

    point[0] = newp[0];
    point[1] = newp[1];
    point[2] = newp[2];
  
}


void TPZCartsys::ToReference (TPZVec<REAL> &point) {
  
    TPZVec<REAL> newp(3,0.);
    newp[0] = fTr(0,0)*point[0] + fTr(0,1)*point[1] + fTr(0,2)*point[2];
    newp[1] = fTr(1,0)*point[0] + fTr(1,1)*point[1] + fTr(1,2)*point[2];
    newp[2] = fTr(2,0)*point[0] + fTr(2,1)*point[1] + fTr(2,2)*point[2];

    point[0] = newp[0] + fOrigin[0];
    point[1] = newp[1] + fOrigin[1];
    point[2] = newp[2] + fOrigin[2];
	if (!fReference) return;
}


void TPZCartsys::CalcRotMatrix(TPZVec<REAL> &angle){
	fTr(0,0)=cos(angle[1])*cos(angle[2]);
	fTr(0,1)=cos(angle[1])*sin(angle[2]);
	fTr(0,2)=sin(angle[1]);
	fTr(1,0)=-(cos(angle[2])*sin(angle[0])*sin(angle[1])) - cos(angle[0])*sin(angle[2]);
	fTr(1,1)=cos(angle[0])*cos(angle[2]) - sin(angle[0])*sin(angle[1])*sin(angle[2]);
	fTr(1,2)=cos(angle[1])*sin(angle[0]);
	fTr(2,0)=-(cos(angle[0])*cos(angle[2])*sin(angle[1])) + sin(angle[0])*sin(angle[2]);
	fTr(2,1)=-(cos(angle[2])*sin(angle[0])) - cos(angle[0])*sin(angle[1])*sin(angle[2]);
	fTr(2,2)=cos(angle[0])*cos(angle[1]);
}


void TPZCartsys::Normalise(TPZVec<REAL> &p){
	 
	REAL norm = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    norm = sqrt(norm);
    if (norm >= 1.0E-10) {
        p[0] /= norm;
        p[1] /= norm;
        p[2] /= norm;
    }
}
	
void TPZCartsys::GetNormal(TPZVec<REAL> &vec1, TPZVec<REAL> &vec2, TPZVec<REAL> &norm){
	norm[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
	norm[1] = vec2[0] * vec1[2] - vec1[0] * vec2[2];
	norm[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];

	Normalise(norm);
}

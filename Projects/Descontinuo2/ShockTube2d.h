#ifndef SHOCKTUBE2DHPP
#define SHOCKTUBE2DHPP

#include "pzreal.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "pzgmesh.h"
#include "pzflowcmesh.h"
#include "pzeulerconslaw.h"

// Creates a mesh for the simple shock problem

void STMeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms);
TPZGeoMesh * CreateSTGeoMesh(TPZVec< TPZVec< REAL > > & nodes,
                             TPZVec< TPZVec< int > > & elms,
                             MElementType ElType, int matId,
                             TPZVec<TPZGeoEl *> & gEls,
                             int nSubdiv);

// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh * STCompMesh(REAL CFL, REAL delta,
                             int degree, int nSubdiv,
                             TPZArtDiffType DiffType,
                             TPZTimeDiscr Diff_TD,
                             TPZTimeDiscr ConvVol_TD,
                             TPZTimeDiscr ConvFace_TD);

#endif
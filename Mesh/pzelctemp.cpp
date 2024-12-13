/**
 * @file
 * @brief Contains the implementation of the TPZIntelGen methods.
 */

#include "pzelctemp.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "pzcmesh.h"

#include "TPZMatSingleSpace.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzintelgen");
#endif

template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel) : TPZRegisterClassId(&TPZIntelGen::ClassId),
TPZInterpolatedElement(mesh,gel){

	//  RemoveSideRestraintsII(EInsert);
	gel->SetReference(this);
  int matid = gel->MaterialId();
#ifdef PZDEBUG
  if (mesh.FindMaterial(matid) == 0) {
      std::cout << "Creating element " << __PRETTY_FUNCTION__ << " matid " << matid << std::endl;
      std::cout << "Material not found, bailing out \n";
    DebugStop();
  }
#endif
  
   bool isMeshBuilt = false;
   AdjustIntegrationRule(isMeshBuilt,fPreferredOrder);
}

template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, int nocreate) : TPZRegisterClassId(&TPZIntelGen::ClassId),
TPZInterpolatedElement(mesh,gel)
{
	fPreferredOrder = -1;
}

template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, const TPZIntelGen<TSHAPE> &copy) : TPZRegisterClassId(&TPZIntelGen::ClassId),
TPZInterpolatedElement(mesh,copy), fIntRule(copy.fIntRule) {
	fPreferredOrder = copy.fPreferredOrder;
}


template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh,
								 const TPZIntelGen<TSHAPE> &copy,
								 std::map<int64_t,int64_t> & gl2lcConMap,
								 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZIntelGen::ClassId), 
TPZInterpolatedElement(mesh,copy,gl2lcElMap), fIntRule(copy.fIntRule)
{
	
	fPreferredOrder = copy.fPreferredOrder;
	int i;
	// for(i=0;i<TSHAPE::NSides;i++)
	// {
	// 	int64_t lcIdx = -1;
	// 	int64_t glIdx = copy.fConnectIndexes[i];
	// 	if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
	// 	else
	// 	{
	// 		std::stringstream sout;
	// 		sout << "ERROR in : " << __PRETTY_FUNCTION__
	// 		<< " trying to clone the connect index: " << glIdx
	// 		<< " wich is not in mapped connect indexes!";
	// 		LOGPZ_ERROR(logger, sout.str().c_str());
	// 		fConnectIndexes[i] = -1;
	// 		return;
	// 	}
	// 	fConnectIndexes[i] = lcIdx;
	// }
  //TODO:NATHANFRAN
}


template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen() :
TPZRegisterClassId(&TPZIntelGen::ClassId), 
TPZInterpolatedElement(), fIntRule() {
	fPreferredOrder = -1;
}

template<class TSHAPE>
MElementType TPZIntelGen<TSHAPE>::Type() {
	return TSHAPE::Type();
}


template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetIntegrationRule(int ord) {
	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	fIntRule.SetOrder(order);
}

/**Sets the interpolation order for the interior of the element*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetInterpolationOrder(int order) {
	fPreferredOrder = order;
}

//TODO:NATHANFRAN ASKPHIL
/**return the preferred order of the polynomial along side iside*/
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::PreferredSideOrder(int side) {
	if(side < TSHAPE::NCornerNodes) return 0;
	if(side<TSHAPE::NSides) {
		int order =fPreferredOrder;
		return AdjustPreferredSideOrder(side,order);
	}
	PZError << "TPZIntelgen::PreferredSideOrder called for side = " << side << "\n";
	return 0;
	
}

template<class TSHAPE>
int64_t TPZIntelGen<TSHAPE>::ConnectIndex(int con) const{
	
#ifndef PZNODEBUG
	if(con<0 || con>= NConnects()) {
		std::cout << "TPZIntelgen::ConnectIndex wrong parameter con " << con <<
		" NSides " << TSHAPE::NSides << " NConnects " << NConnects() << std::endl;
		DebugStop();
	}
	
#endif

	return this->ConnectVec()[con];
}



/**Sets the preferred interpolation order along a side
 This method only updates the datastructure of the element
 In order to change the interpolation order of an element, use the method PRefine
 */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetPreferredOrder(int order)
{
	fPreferredOrder = order;
}


/**returns the actual interpolation order of the polynomial for a connect*/
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::ConnectOrder(int connect) const {
	return Connect(connect).Order();
}

/** Returns the transformation which transform a point from the side to the interior of the element */
template<class TSHAPE>
TPZTransform<> TPZIntelGen<TSHAPE>::TransformSideToElement(int side) {
	return TSHAPE::TransformSideToElement(side);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZInterpolatedElement::Write(buf,withclassid);
	TPZManVector<int,3> order(3,0);
	fIntRule.GetOrder(order);
	buf.Write(order);
	buf.Write(&fPreferredOrder,1);
	int classid = this->ClassId();
	buf.Write ( &classid, 1 );
}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZInterpolatedElement::Read(buf,context);
	TPZManVector<int,3> order;
	buf.Read(order);
	fIntRule.SetOrder(order);
	buf.Read(&fPreferredOrder,1);
	int classid = -1;
	buf.Read( &classid, 1 );
	if ( classid != this->ClassId())
	{
		std::stringstream sout;
		sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " for an package of id = " << classid;
		LOGPZ_ERROR ( logger, sout.str().c_str() );
	}
}

using namespace pzshape;

#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapepiramHdiv.h"
template<>
void TPZIntelGen<pzshape::TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension /* && Material()->Id() > 0 */) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}

#include "TPZRefCube.h"

#include "TPZRefLinear.h"
#include "pzrefquad.h"

#include "pzgeoquad.h"

#include "pzreftriangle.h"
#include "pzgeotriangle.h"

#include "pzrefprism.h"
#include "pzgeoprism.h"

#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"

#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"



#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"
#include "tpzpyramid.h"

template class TPZIntelGen<TPZShapeTriang>;
template class TPZIntelGen<TPZShapePoint>;
template class TPZIntelGen<TPZShapeLinear>;
template class TPZIntelGen<TPZShapeQuad>;
template class TPZIntelGen<TPZShapeTetra>;
template class TPZIntelGen<TPZShapePrism>;
template class TPZIntelGen<TPZShapePiram>;
template class TPZIntelGen<TPZShapeCube>;


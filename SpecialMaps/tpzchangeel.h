/**
 * @file
 * @brief Contains the TPZChangeEl class. It is a special map.
 */

#ifndef TPZCHANGEEL_H
#define TPZCHANGEEL_H

#include "pzgeotriangle.h"
#include "tpzgeoblend.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * @ingroup geometry
 * @since 2007
 * @brief Special map. It is util to convert a linear element for quadratic element, but the same topology. \ref geometry "Geometry"
 */
/** Also turns slide middle nodes of an quadratic element to the quarterpoint with respect to a indicated side. */
class TPZChangeEl {
public:

    TPZChangeEl();
   ~TPZChangeEl();

    /** @brief Turns an linear geoelement to quadratic */
    static TPZGeoEl * ChangeToQuadratic(TPZGeoMesh *Mesh, int64_t ElemIndex);

    /** @brief Turns a regular element into a geoblend */
    static TPZGeoEl * ChangeToGeoBlend(TPZGeoMesh *Mesh, int64_t ElemIndex);

    /** @brief Turns a regular linear element into a TPZArc3D*/
    static TPZGeoEl * ChangeToArc3D(TPZGeoMesh *mesh, const int64_t ElemIndex,
                                    const TPZVec<REAL> &xcenter, const REAL radius);
    /** @brief Turns a regular 2D element into a TPZCylinderMap(using rotation matrix)*/
    static TPZGeoEl * ChangeToCylinder(TPZGeoMesh *mesh, const int64_t ElemIndex,
                                       const TPZVec<REAL> &xcenter,
                                       const TPZFMatrix<REAL> &rotmatrix,
                                       const REAL radius);
    /** @brief Turns a regular 2D element into a TPZCylinderMap(using cylinder axis)*/
    static TPZGeoEl * ChangeToCylinder(TPZGeoMesh *mesh, const int64_t ElemIndex,
                                       const TPZVec<REAL> &xcenter,
                                       const TPZVec<REAL> &axis,
                                       const REAL radius);
    
    /** @brief Slide middle nodes of an quadratic geoelement to the quarterpoint with respect to a given side */
    static TPZGeoEl * ChangeToQuarterPoint(TPZGeoMesh *Mesh, int64_t ElemIndex, int targetSide);

    /**
	 * @brief Return if a given point x is near to some node of a given geo element
	 * @param gel [in] given geo element
	 * @param x [in] given point
	 * @param meshNode [out] id of node that is in the x range
	 * @param tol [in] x range radius
	 */
	static bool NearestNode(TPZGeoEl * gel, TPZVec<REAL> &x, int64_t &meshNode, double tol);
    /**
	 * @brief Return the id of the node into the geometric mesh nearest a given point x
	 * @param gmesh [in] given geometric mesh
	 * @param x [in] given point
	 * @param tol [in] x range radius
	 */
    static int64_t NearestNode(TPZGeoMesh * gmesh, TPZVec<REAL> &x, double tol);
    static bool CreateMiddleNodeAtEdge(TPZGeoMesh *Mesh, int64_t ElemIndex, int edge, int64_t &middleNodeId);

    /**
       @brief Stores the first neighbour for each side of a given element
       @note If, for a given side, the neighbour is itself, an empty TPZGeoElSide is placed in its position
       @param gel [in] given geometric element
       @param neighs [out] vector of neighbours
     */
    static void StoreNeighbours(TPZGeoEl* gel, TPZVec<TPZGeoElSide> &neighs);
    /**
       @brief Restores neighbourhood information, checking for empty TPZGeoElSides in the vector
       @param gel [in/out] given geometric element
       @param neighs [in] vector of neighbours
     */
    static void RestoreNeighbours(TPZGeoEl* gel, TPZVec<TPZGeoElSide> &neighs);

    /**
       @brief Change ordering of element's nodes and adjust connectivity
       @param gel [in] given geometric element
       @param ordering [in] element nodes in the new ordering
       @note NOT YET TESTED!!!!
     */
    static void ChangeNodeOrdering(TPZGeoEl *gel, const TPZVec<int64_t> & mapped_nodes);
};

#endif

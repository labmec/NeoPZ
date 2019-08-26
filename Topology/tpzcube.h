/**
 * @file
 * @brief Contains the TPZCube class which defines the topology of the hexahedron element.
 */

#ifndef PZTOPOLOGYTPZCUBE_H
#define PZTOPOLOGYTPZCUBE_H


#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"
#include "pznumeric.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZIntPoints;
class TPZIntCube3D;
class TPZGraphElQ3dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of the hexahedron element. \ref topology "Topology"
	 * Sides 0 to 7 are vertices, sides 8 to 19 are lines, 20 to 25 are quadrilaterals 
	 * and side 26 is the hexahedra (cube).
	 */
	class TPZCube : public TPZSavable {
	public:
		
		/** @brief enumerate for topological characteristics */
		enum {NSides = 27, NCornerNodes = 8, Dimension = 3, NFaces = 6};
		
                int ClassId() const override;
                void Read(TPZStream &buf, void *context) override;
                void Write(TPZStream &buf, int withclassid) const override;
                
		/** @brief Default constructor */
        TPZCube() : TPZRegisterClassId(&TPZCube::ClassId) {
		}
		
		/** @brief Default destructor */
		virtual ~TPZCube() {
		}
		
		/** @name About sides of the topological element
		 * @{ */
		
		/** @brief Returns the dimension of the side */
		static int SideDimension(int side);
		
		/** @brief Get all sides with lower dimension on side */	
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		/** @brief Get all sides with lower dimension but equal to DimTarget on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side Smaller dimension side
		 * @param high Vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);
		
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NSideNodes(int side);
		/** @brief Returns the local node number of the node "node" along side "side" */
		static int SideNodeLocId(int side, int node);

		/** @brief Returns number of connects of the element (27)  ??? */
		static int NumSides();
		/** @brief Returns the number of connects for a set dimension */
		static int NumSides(int dimension);
		
		/** @brief Returns the number of connectivities associated with a side */
		static int NContainedSides(int side);
		/** @brief Returns the local connect number of the connect "c" along side "side" */
		static int ContainedSideLocId(int side, int c);

        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        /**
         * This method calculates the influence (a.k.a. the blend function) of the side side regarding an
         * interior point qsi. It is used by the TPZGeoBlend class.
         * @param side the index of the side
         * @param xi coordinates of the interior point
         * @param correctionFactor influence (0 <= correctionFactor <= 1)
         * * @param corrFactorDxi derivative of the correctionFactor in respect to xi
         */
        template<class T>
        static void CalcSideInfluence(const int &side, const TPZVec<T> &xi, T &correctionFactor,
                                      TPZVec<T> &corrFactorDxi);
		/** @} */

		/** @name About points at the parametric spaces
		 * @{ */
		
		/** @brief Returns the barycentric coordinates in the master element space of the original element */
		static void CenterPoint(int side, TPZVec<REAL> &center);
		
		/** @brief Verifies if the parametric point pt is in the element parametric domain
		 */
		static bool IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol = 1e-6L);
        #ifdef _AUTODIFF
        template<typename T,
                typename std::enable_if<std::is_same<T,Fad<REAL>>::value>::type* = nullptr>
        static bool IsInParametricDomain(const TPZVec<T> &pt, REAL tol){
            TPZVec<REAL> qsiReal(pt.size(),-1);
            for(int i = 0; i < qsiReal.size(); i++) qsiReal[i] = pt[i].val();
            return IsInParametricDomain(qsiReal,tol);
        }
        #endif
        
        /** @brief Generates a random point in the master domain */
        static void RandomPoint(TPZVec<REAL> &pt);
        

        template<class T>
        static bool MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide);
        
        static void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
		
		/** @} */
		
		/** @name About type of the topological element
		 * @{ */
		
		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static MElementType Type();// { return ECube;}
		
		/** @brief Returns the type of the element side as specified in file pzeltype.h */
		static MElementType Type(int side);
		
		/** @} */
		
		/** @name About Transformations
		 * @{ */
				
		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom to the side sideto
		 * @param sidefrom Side where the point resides
		 * @param sideto Side whose closure contains sidefrom
		 * @see the class TPZTransform
		 */
		static TPZTransform<> SideToSideTransform(int sidefrom, int sideto);

		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=26)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformSideToElement(int side);
		/**
		 * @brief Returns the transformation which projects a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=26)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformElementToSide(int side);
		
		/**
		 * @brief Method which identifies the transformation based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(TPZVec<int64_t> &id);
		
		/**
		 * @brief Method which identifies the transformation of a side based on the IDs of the corner nodes
		 * @param side Index of side
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */	
		static int GetTransformId(int side, TPZVec<int64_t> &id);

		/** @} */
		
		/** @name Methods related over numeric integration
		 * @{ */
		
		/**
		 * @brief Create an integration rule over side
		 * @param side Side to create integration rule
		 * @param order Order of the integration rule to be created
		 */
		static TPZIntPoints *CreateSideIntegrationRule(int side, int order);
		
		/** @brief Typedef to numerical integration rule */
		typedef TPZIntCube3D IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphElQ3dd GraphElType;

		/** @} */
		
		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather);

		/** @brief Volume of the master element*/
		static REAL RefElVolume(){return 8.0;}
        
        /* Given side and gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        
        /// Compute the directions of the HDiv vectors
        static void ComputeDirections(TPZFMatrix<REAL> &gradx, REAL detjac, TPZFMatrix<REAL> &directions);
        
        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao);
        
        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors);
        
        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();
	
	protected:
		/** @name Data structure which defines the hexahedral transformations */
		/** @{ */
		
		/** @brief Nodes over quadrilateral sides (2d - faces). */
		static int FaceNodes[6][4];
		
		/** @brief Nodes over lines sides (1d) */
		static int SideNodes[12][2];
		
		/** @brief Ids of the shape face */
		static int ShapeFaceId[6][2];
		
		/** @} */
        
		
	};
	
}

#endif

/**
 * @file
 * @brief Contains declaration of TPZElementMatrix struct which defines the interface for associating an element matrix with the coeficients of its contribution in the global stiffness matrix
 */

#ifndef ELMATHPP
#define ELMATHPP

#include "pzblock.h"
#include "pzconnect.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzstack.h"


class TPZBaseMatrix;



/**
 * @brief This class associates an element matrix with the coeficients of its contribution in the global stiffness matrix. \ref interpolation "Aproximation space"
 * @ingroup interpolation
 */
/**
 This class groups all information associated with an element stiffness matrix so that it can be used independent of the element object itself
 Objects of this class provide storage as well for the constrained stiffness matrix, i.e. the stiffness matrix from which the constrained connects have been eliminated
 In future versions, the computation of the contraints will be incorporated in a method of this class
 */
struct TPZElementMatrix {
    enum MType{Unknown = 0, EF = 1, EK = 2};
	/// Reset the data structure
	void Reset(TPZCompMesh *mesh = NULL, MType type=Unknown)
	{
      fMesh = mesh;
      fType = type;
      fConnect.Resize(0);
      fConstrConnect.Resize(0);
      fDestinationIndex.Resize(0);
      fSourceIndex.Resize(0);
	}
	
	TPZElementMatrix(TPZCompMesh *mesh, MType type) :
        fType(type), fMesh(mesh), fConnect(),
        fConstrConnect(),
        fDestinationIndex(), fSourceIndex()
    {
    }

    TPZElementMatrix() :
        fType(Unknown), fMesh(NULL), fConnect(), fConstrConnect(), 
        fDestinationIndex(), fSourceIndex()
    {}
    TPZElementMatrix &operator=(const TPZElementMatrix &copy);
	
    TPZElementMatrix(const TPZElementMatrix &copy);

	virtual ~TPZElementMatrix(){}
	
	/** @brief Returns the number of nodes of TElementMatrix*/
	int NConnects(){
		return fConnect.NElements();
	}
	
	/** @brief Returns the pointer to the ith node of the element*/
	int64_t ConnectIndex(int i)const{
		return fConnect[i];
	}
	
	virtual void Print(std::ostream &out) = 0;
	
	virtual void SetMatrixSize(short NumBli, short NumBlj, short BlSizei, short BlSizej) = 0;
	
	virtual void SetMatrixMinSize(short NumBli, short NumBlj, short BlMinSizei, short BlMinSizej) = 0;
	
	void ComputeDestinationIndices();
    
    /** @brief permute the order of the connects */
    virtual void PermuteGather(TPZVec<int64_t> &permute) = 0;
	
	
	/** @brief Returns true if the element has at least one dependent node. Returns false otherwise */
	bool HasDependency();
	
	/** @brief Apply the constraints applied to the nodes by transforming the tangent matrix and right hand side */
	virtual void ApplyConstraints() = 0;
    
    
    /// Compute the dependency order of the connects, considering the one shape restraints
    void BuildDependencyOrder(TPZVec<int64_t> &connectlist, TPZVec<int> &DependenceOrder, TPZCompMesh &mesh);

    virtual TPZBaseMatrix &Matrix() = 0;
    virtual TPZBaseMatrix &ConstrMatrix() = 0;
    virtual TPZBlock &Block() = 0;
    virtual TPZBlock &ConstrBlock() = 0;

    
	MType fType;
	
	TPZCompMesh * fMesh;
	/** @brief Vector of pointers to TPZConnect objects*/
	TPZStack<int64_t> fConnect;
	/** @brief Vector of all nodes connected to the element*/
	TPZStack<int64_t> fConstrConnect;
	
	TPZManVector<int64_t> fDestinationIndex, fSourceIndex;
};

#endif

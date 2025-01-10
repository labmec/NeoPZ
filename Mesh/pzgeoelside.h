/**
 * @file
 * @brief Contains declaration of TPZGeoElSide class which represents an element and its side, and TPZGeoElSideIndex class which represents an TPZGeoElSide index.
 */

#ifndef PZGEOELSIDEH
#define PZGEOELSIDEH

/*******       TPZGeoElSide       *******/


#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"
#include "tpzintpoints.h"
#include <set>

#include "fadType.h"

class TPZCompElSide;

class TPZGeoElSide;
class TPZGeoMesh;
class TPZGeoEl;

/**
 * @ingroup geometry
 * @brief Utility class which represents an element index with its side. \ref geometry "Geometry"
 */
class TPZGeoElSideIndex : public TPZSavable{
private:
	int64_t fGeoElIndex{-1};
	int fSide{-1};
	
public:
	/** @brief Destructor. */
	~TPZGeoElSideIndex();
	/** @brief Simple constructor. */
	TPZGeoElSideIndex();
	/** @brief Constructor with geometric element referenced and corresponding side. */
	TPZGeoElSideIndex(TPZGeoEl *gel,int side);
	
	TPZGeoElSideIndex(int64_t gelindex,int side);
	
	TPZGeoElSideIndex(const TPZGeoElSide &side);
	
	TPZGeoElSideIndex(const TPZGeoElSideIndex &cp);
	/** @brief To clone current object */
	TPZGeoElSideIndex * Clone();
	/** @brief Redefines operator = attribuition to TPZGeoElSideIndex object. */
	TPZGeoElSideIndex &operator= (const TPZGeoElSideIndex &A );
	
    /** @brief cast to bool to indicate whether it is an initialized side */
    operator bool() const
    {
        return fGeoElIndex == -1 || fSide == -1;
    }
    
    bool operator==(const TPZGeoElSideIndex &copy)
    {
        return (fGeoElIndex == copy.fGeoElIndex && fSide == copy.fSide);
    }
    
	int Side() const;
	
	void SetSide(int side);
	
	TPZGeoEl *Element(const TPZGeoMesh *mesh) const;
	
	void SetElement(TPZGeoEl* geoel);
	
	int64_t ElementIndex() const;
	
	void SetElementIndex(int64_t i);
	        int ClassId() const override;
        void Read(TPZStream &buf, void *context) override;
        void Write(TPZStream &buf, int withclassid) const override;
};

/**
 * @brief Utility class which represents an element with its side. \ref geometry Geometry
 * @ingroup geometry
 */
/** This class is often used to manipulate neighbouring information between elements */
class TPZGeoElSide : public TPZSavable {
	
	TPZGeoEl *fGeoEl{nullptr};
	int fSide{-1};
public:
    
    /// return the number of element/side pairs which compose the current set of points
	int NSubElements();
    
    /// build the list of element/side pairs which compose the current set of points
	void GetSubElements2(TPZStack<TPZGeoElSide> &subelements);
    
    /// returns the father/side pair which contains this/side and is strictly larger than this/side
	TPZGeoElSide StrictFather() const;
    
    /// returns the father/side pair which contains the set of points associated with this/side
	TPZGeoElSide Father2() const;
    
    /// return the element/side pair which contains this/side and has a computational element associated
	TPZCompElSide LowerLevelCompElementList2(int onlyinterpolated);
	
	/** @brief Checks whether other is a relative (son or ancestor) of this */
	bool IsRelative(const TPZGeoElSide &other) const;
	
	/** @brief Checks whether other is an ancestor of this */
	bool IsAncestor(const TPZGeoElSide &other) const;
    
    /** @brief Checks whether other is a neighbour of the current element */
    bool IsNeighbour(const TPZGeoElSide &other) const
    {
        return NeighbourExists(other);
    }
    
    /** @brief X coordinate of a point loc of the side */
	void X(TPZVec< REAL > &loc, TPZVec< REAL > &result) const;
    
    /** @brief Parametric coordinate of a point loc of the side and return parametric element point */
    void QsiElement(TPZVec< REAL > &qsi_side, TPZVec< REAL > &qsi_element) const;
    
    /** @brief X coordinate of a point loc of the side */
    void GradX(TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx) const;

    bool ResetBlendConnectivity(const int64_t &index);
	
    /** @brief X coordinate of a point loc of the side */
    void X(TPZVec< Fad<REAL> > &loc, TPZVec< Fad<REAL> > &result) const;
    
    /** @brief X coordinate of a point loc of the side */
    void GradX(TPZVec< Fad<REAL> > &loc, TPZFMatrix< Fad<REAL> > &gradx) const;

    /** @brief Jacobian associated with the side of the element */
	void Jacobian(TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;
    
    /** @brief Area associated with the side */
    REAL Area() const;
	
	/** @returns Total number of neighbours through this side*/
	int NNeighbours() const;

	/** @brief Get number of neighbours of a given dimension */
	int NNeighbours(int dimfilter) const;

	/** @brief Get number of neighbours filtered by dimension and/or material id
	 * @param dimfilter: only return elements of this dimension. Ignore filter if set to < 0;
	 * @param matids: only return elements of a material id contained in this set. Ignore filter if passing empty set;
	*/
	int NNeighbours(int dimfilter, const std::set<int>& matids) const;
	
	/** @brief Returns the number of neighbours, excluding the given element (thisElem) */
	int NNeighboursButThisElem(TPZGeoEl *thisElem) const;
	
	TPZGeoElSide(){ fGeoEl = 0; fSide  = -1;}
	
	TPZGeoElSide(TPZGeoEl *gel,int side){  fGeoEl = gel; fSide = side;}
	
	/** @brief This constructor set an TPZGeoElSide based in the cornerNodes of an side of gel */
	/** If the cornerNodes are not consistent, the TPZGeoElSide created is NULL */
	TPZGeoElSide(TPZGeoEl *gel, std::set<int64_t> &sideCornerNodes);
	
	TPZGeoElSide(const TPZGeoElSideIndex &index, const TPZGeoMesh * mesh){
		this->fSide = index.Side();
		this->fGeoEl = index.Element(mesh);
		if(fGeoEl == 0 && index.ElementIndex() != -1)
        {
		    DebugStop();
        }
	}
    
    TPZGeoElSide(int zero) : fGeoEl(0), fSide(-1)
    {
        
    }
    
    TPZGeoElSide(TPZGeoEl *gel);
	
	TPZGeoEl *Element()const{return fGeoEl;}
	
	void SetElement(TPZGeoEl* geoel){ fGeoEl = geoel; }
    
    /** @brief print geometric characteristics of the element/side */
    void Print(std::ostream &out) const;
	
	int Side() const {return fSide;}
	
	void SetSide(int side) { fSide = side; }
	
	bool IsLinearMapping() const;
	
	int Exists() const {return (fGeoEl != 0 && fSide > -1);}
	
    /** @brief return the coordinates of the center in master element space (associated with the side) */
	void CenterPoint(TPZVec<REAL> &center) const;
    
    /** @brief return the coordinates of the center of the side in real space */
	void CenterX(TPZVec<REAL> &Xcenter) const;
    
    /** @brief compute the normal to the point from left to right neighbour */
    void Normal(TPZVec<REAL> &point, TPZGeoEl *left, TPZGeoEl *right, TPZVec<REAL> &normal) const;
    
    /** @brief compute the normal to the point */
    void Normal(TPZVec<REAL> &point, TPZVec<REAL> &normal) const;
    
    /** @brief Returns the number of sides in which the current side can be decomposed */
    int NSides() const;
	
	TPZGeoElSide Neighbour() const;//return neighbour of the side fSide
	
	/** @brief Returns the set of neighbours which can directly be accessed by the datastructure */
	void AllNeighbours(TPZStack<TPZGeoElSide> &allneigh);
	
	/** @brief Returns the set of neighbours as computed by the intersection of neighbours along vertices */
	void ComputeNeighbours(TPZStack<TPZGeoElSide> &compneigh);
	
	int64_t Id();
    
    /** @brief the dimension associated with the element/side */
	int Dimension() const;
	
	int operator==(const TPZGeoElSide &other) const {
		return fGeoEl == other.fGeoEl && fSide == other.fSide;
	}
	int operator!=(const TPZGeoElSide &other) const {
		return fGeoEl != other.fGeoEl || fSide != other.fSide;
	}
	
	int operator<(const TPZGeoElSide &other) const {
		return (fGeoEl < other.fGeoEl || (fGeoEl == other.fGeoEl && fSide < other.fSide));
	}
	
	int operator>(const TPZGeoElSide &other) const {
		return (fGeoEl > other.fGeoEl || (fGeoEl == other.fGeoEl && fSide > other.fSide));
	}

	/** @brief Next neighbour operator as post-increment */
	TPZGeoElSide operator++(int){
		TPZGeoElSide pre = *this;
		*this = this->Neighbour();
		return pre;
	}
	/** @brief Next neighbour operator as pre-increment */
	TPZGeoElSide& operator++(){
		*this = this->Neighbour();
		return *this;
	}
    /** @brief Next neighbour operator as post-increment */
    TPZGeoElSide operator--(int){
        TPZGeoElSide pre = *this;
        TPZGeoElSide neigh = *this;
        while(neigh.Neighbour() != *this) neigh = neigh.Neighbour();
        *this = neigh;
        return pre;
    }
    /** @brief Next neighbour operator as pre-increment */
    TPZGeoElSide& operator--(){
        TPZGeoElSide neigh = *this;
        while(neigh.Neighbour() != *this) neigh = neigh.Neighbour();
        *this = neigh;
        return *this;
    }

    /** @brief The conversion to bool indicates whether the object has an associated element */
    operator bool() const
    {
        return fGeoEl != 0;
    }
	
	/** @brief Accumulates the transformations from the current element/side to the neighbour/side
	 * @note Third improved version */
	void SideTransform3(TPZGeoElSide neighbour,TPZTransform<> &t);
	
	/// @brief Insert neighbour in the connectivity loop
	/// @param neighbour will be my direct neighbour
	void SetConnectivity(const TPZGeoElSide &neighbour) const;
    
	/** @brief This method inserts the element/side and all lowerdimension sides into the connectivity loop */
	void InsertConnectivity(TPZGeoElSide &neighbour);
	
    /** @brief This method inserts the element/side and all lowerdimension sides into the connectivity loop
                neighbour : permuted element neighbour
                mapsides : indicates the permutation of sides between neighbours
     */
    void InsertConnectivity(TPZGeoElSide &neighbour, const TPZVec<int> &mapsides);

    /// Remove the element/side from the connectivity loop
	// the neighbour of the element/side will be undefined
	void RemoveConnectivity();
	
	static void BuildConnectivities(TPZVec<TPZGeoElSide> &elvec, TPZVec<TPZGeoElSide> &neighvec);
	
	/** @brief Fill in the data structure for the neighbouring information*/
	void SetNeighbour(const TPZGeoElSide &neighbour) const;
	
	TPZTransform<REAL> NeighbourSideTransform(const TPZGeoElSide &neighbour);
	
    /// Compute if the neighbour along a face has normal pointing outward
    bool IsNeighbourCounterClockWise(TPZGeoElSide &neighbour);
    
	/** 
	 * @brief Compute the transformation between the master element space of one side
	 * of an element to the master element space of a higher dimension side
	 */
	TPZTransform<REAL> SideToSideTransform(TPZGeoElSide &higherdimensionside);
    
    /// return the lowest level direct ancestor
    TPZGeoElSide LowestFatherSide();
    
    /// return the TPZGeoElSide element that contains the current element/side
    TPZGeoElSide LowerLevelSide() const;
    
    /**
     * @brief [deprecated] use YoungestChildren
     */
    virtual void GetAllSiblings(TPZStack<TPZGeoElSide> &unrefinedSons);
    /**
     * @brief This method will return all children at the bottom of the refinement tree of the element. i.e. all children that have no subelements
     */
    virtual void YoungestChildren(TPZStack<TPZGeoElSide> &unrefinedSons);

	/** @brief Returns a pointer to the elementside referenced by the geometric elementside*/
	TPZCompElSide Reference() const;
	/** @brief Return 1 if the element has subelements along side*/
	int HasSubElement();
	
	/* @brief Return the number of nodes for a particular side*/
	int NSideNodes() const;
	
	/** @brief Returns the index of the nodenum node of side*/
	int64_t SideNodeIndex(int nodenum) const;
	
	/** @brief Returns the index of the local nodenum node of side*/
	int64_t SideNodeLocIndex(int nodenum) const;
    
	
	/** @brief Returns 1 if neighbour is a neighbour of the element along side*/
	int NeighbourExists(const TPZGeoElSide &neighbour) const;
    
    /**
      *      verify if a neighbour with the given material id exists
     */
    TPZGeoElSide HasNeighbour(int materialid) const;

    /**
     *      verify if a neighbour with the given material id exists
    */
    TPZGeoElSide HasNeighbour(const std::set<int> &matIDs) const;
    
    /** verifiy if a larger (lower level) neighbour exists with the given material id
     */
    TPZGeoElSide HasLowerLevelNeighbour(int materialid) const;
    
    /** verifiy if a larger (lower level) neighbour exists with the given material id
     */
    TPZGeoElSide HasLowerLevelNeighbour(const std::set<int> &matids) const;
    

    /** @brief Will return all elements of equal or higher level than than the current element */
	void EqualorHigherCompElementList2(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates);
    
	/** @brief Pushes all connected computational elements which have higher dimension than the current element/side \n
	 * if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack \n
     * if removeduplicates == 1 no elements which are direct neighbours will be put on the stack */
	void HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated);
	
	/** @brief Returns all connected computational elements which have level higher to the current element \n
     * if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack \n
     * if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void HigherLevelCompElementList2(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated, int removeduplicates);
	
	/** @brief Returns all connected computational elements to the current element \n
	 * if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack \n
     * if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void ConnectedCompElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates);
	
	/** @brief Returns all connected computational elements which have equal level to the current element */
	/** This method will not put this on the stack \n
     * if onlyinterpolated == 1 only elements TPZInterpolatedElement will be put on the stack \n
     * if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,	int onlyinterpolated, int removeduplicates);
    
    /** @brief Will return all elements of equal or higher level than than the current element \n
     * if onlymultiphysicelement == 1 only elements TPZMultiphysicsElement will be put on the stack \n
     * if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void EqualorHigherCompElementList3(TPZStack<TPZCompElSide> &celside, int onlymultiphysicelement, int removeduplicates);
    
    /** @brief Returns all connected computational elements which have equal level to the current element */
	/** This method will not put this on the stack \n
     * if onlymultiphysicelement == 1 only elements TPZMultiphysicsElement will be put on the stack \n
     * if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void EqualLevelCompElementList3(TPZStack<TPZCompElSide> &elsidevec,	int onlymultiphysicelement, int removeduplicates);
    
    /** @brief Returns all connected computational elements which have level higher to the current element \n
     * if onlymultiphysicelement == 1 only elements TPZMultiphysicsElement will be put on the stack \n
     * if removeduplicates == 1 no elements which are direct neighbours will be put on the stack*/
	void HigherLevelCompElementList3(TPZStack<TPZCompElSide> &elsidevec, int onlymultiphysicelement, int removeduplicates);


    /* @brief Creates an integration rule for the topology of this side */
	TPZIntPoints * CreateIntegrationRule(int order);

    int GelLocIndex(int index) const;
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;
};

/** @brief Overload operator << to print geometric element side data */
std::ostream  &operator << (std::ostream & out,const TPZGeoElSide &geoside);

inline void TPZGeoElSide::AllNeighbours(TPZStack<TPZGeoElSide> &allneigh) {
	TPZGeoElSide neigh = Neighbour();
#ifndef PZNODEBUG
	if(! Exists() || ! neigh.Exists()) 
    {
		std::cout << "TPZGeoElSide AllNeighbours inconsistent\n";
		return;
    }
#endif
	while(neigh != *this)
    {
		allneigh.Push(neigh);
		neigh = neigh.Neighbour();
    }
}


/** Implementing TPZGeoElSideIndex methods */

inline TPZGeoElSideIndex::~TPZGeoElSideIndex(){
	//nothing to be done here
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(){
    this->fSide = -1;
    this->fGeoElIndex = -1;
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(int64_t gelindex,int side){  
    this->fGeoElIndex = gelindex;
    this->fSide = side;
}

inline TPZGeoElSideIndex::TPZGeoElSideIndex(const TPZGeoElSideIndex &cp){
    this->operator =(cp);
}

inline TPZGeoElSideIndex * TPZGeoElSideIndex::Clone(){
    return new TPZGeoElSideIndex(*this);
}

inline TPZGeoElSideIndex & TPZGeoElSideIndex::operator= (const TPZGeoElSideIndex &A ){
    this->fGeoElIndex = A.fGeoElIndex;
    this->fSide = A.fSide;
    return *this;
}

inline int TPZGeoElSideIndex::Side() const{
    if (this->fSide == -1 || this->fGeoElIndex == -1){
		return -1;
    }
    return this->fSide;
}

inline void TPZGeoElSideIndex::SetSide(int side){
    this->fSide = side;
}  

inline int64_t TPZGeoElSideIndex::ElementIndex() const{
    return this->fGeoElIndex;
}

inline void TPZGeoElSideIndex::SetElementIndex(int64_t i){
    this->fGeoElIndex = i;
}

#endif

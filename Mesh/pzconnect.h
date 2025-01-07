/**
 * @file
 * @brief Contains declaration of TPZConnect class which represents a set of shape functions associated with a computational element
 */

#ifndef  PZCONNECTH
#define  PZCONNECTH

#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>
#include <set>


class TPZBndCond;
class TPZCompMesh;
class TPZBlock;

/** 
 * @brief Represents a set of shape functions associated with a computational element/side. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/**
 * This class keeps track of information associated with an element/side such as order of interpolation
 * sequence number in the vector of blocks of equations \n
 * Objects of this class also contain the information necessary for constraints between shapefunctions
 */
class TPZConnect : public TPZSavable {
public:
    enum EConnectType {ENone = 0, EPressure = 1, ECondensed = 2};
	/** @brief Node block number */
	int64_t		fSequenceNumber;
	
	/** @brief Number of element connected */
	int		fNElConnected;
	    
    /** @brief Flag to signal the type of connect 
     * ECondensed : the connect will be condensed before assembly 
     * EPressure : a pressure connect which should be assembled after the flux connects
     * 
     */
    union
    {
        unsigned int fFlags;
        struct
        {
            /** @brief Interpolation order of the associated shape functions */
            unsigned char fOrder;
            /** @brief Number of state variables associated with each shape function */
            unsigned char fNState;
            /** @brief Whether the connnect is associated with a lagrange multiplier, indicating the level 
             * level n multipliers need to be numbered after the multipliers of level n-1
            */
            unsigned char fLagrangeMultiplier;
            /** @brief Whether the equations associated with the connect should be/are condensed */
            bool fIsCondensed;
        } fCompose;
    };
    /** @brief Number of shape functions associated with the connect */
    unsigned int fNShape;

	
public:
	/** @brief Structure to reference dependency */
	class TPZDependBase : public TPZSavable {
        public :
		int64_t			fDepConnectIndex;
		TPZDependBase		*fNext;
		
		TPZDependBase();
		
		TPZDependBase(const TPZDependBase &copy);
		TPZDependBase(int64_t connectindex);
		
		virtual ~TPZDependBase();
		TPZDependBase *HasDepend(int64_t DepConnectIndex);
		TPZDependBase *RemoveDepend(TPZDependBase *Ptr);
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;

    virtual int DepRows() = 0;
    virtual int DepCols() = 0;
    virtual TPZDependBase* Clone() = 0;
    virtual TPZDependBase* CreateEmptyInstance() = 0;
    virtual void PrintMat(std::ostream &out) = 0;
		/**
		 * @brief Copy a depend data structure to a clone depend in a clone mesh
		 * @param orig original depend to be copied
		 * @param gl2lcIdx global to local indexes map
		 */
		virtual void CopyFrom(TPZDependBase *orig , std::map<int64_t,int64_t>& gl2lcIdx);
	};

  template<class TVar>
  class TPZDepend : public TPZDependBase{
  public:
    TPZFNMatrix<50,TVar> fDepMatrix;
    TPZDepend(int64_t DepConnectIndex,TPZFMatrix<TVar> &depmat,int64_t ipos,int64_t jpos, int isize, int jsize);
    TPZDepend() = default;
    TPZDepend(const TPZDepend &copy) = default;
    ~TPZDepend() = default;
    
    int DepRows() override {return fDepMatrix.Rows();}
    int DepCols() override {return fDepMatrix.Cols();}

    TPZFMatrix<TVar> const &GetDepMatrix(){
      return fDepMatrix;
    }
    

    TPZDepend<TVar>* CreateEmptyInstance() override{
      return new TPZDepend<TVar>;
    }
    
    TPZDepend<TVar>* Clone() override{
      auto *newdep = new TPZDepend<TVar>(*this);
      return newdep;
    }

    void PrintMat(std::ostream &out) override{
      fDepMatrix.Print(" ",out);
    }
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;
    void CopyFrom(TPZDependBase *orig , std::map<int64_t,int64_t>& gl2lcIdx) override;
  };
private:
    /// Dependency list of the connect (NULL if the connect is not constrained)
	TPZDependBase *fDependList;
	
public:
	/** @brief Default constructor */
	TPZConnect();
    
    /** @brief copy constructor */
    TPZConnect(const TPZConnect &copy);
	/** @brief Default destructor */
	~TPZConnect();
	
	TPZConnect &operator=(const TPZConnect &con);
        
        virtual int ClassId() const override{
            return Hash("TPZConnect");
        }
	
    /** @brief Reset the data of the connect */
    void Reset()
    {
        if (fDependList) {
            DebugStop();
            delete fDependList;
            fDependList = 0;
        }
        SetSequenceNumber(-1);
        SetNState(0);
        SetOrder(0,-1);
        SetNShape(0);
        ResetElConnected();
        SetCondensed(false);
        SetLagrangeMultiplier(0);
    }
	/**
	 * @brief Number of degrees of freedom associated with the object
	 * @note In order to compute NDof, the Connect object needs to know the mesh
	 */
	int NDof(TPZCompMesh &mesh);
    
	/**
	 * @brief Number of degrees of freedom associated with the object
	 * @note In order to compute NDof, the Connect object needs to know the mesh
	 */
	int NDof() const
    {
        return fNShape*fCompose.fNState;
    }
    
    /**
     * @brief Number of state variables associated with the connect
     */
    unsigned char NState() const
    {
        return fCompose.fNState;
    }
    
    unsigned int NShape() const
    {
        return fNShape;
    }

	/** @brief Returns the Sequence number of the connect object */
	/** If the \f$ sequence number == -1 \f$ this means that the node is unused */
	int64_t SequenceNumber() const
    {
        return fSequenceNumber;
    }

	/** @brief Set the sequence number for the global system of equations of the connect object*/
	/** If the argument \f$i==-1\f$ this means that the node is out of use*/
	void SetSequenceNumber(int64_t i) {fSequenceNumber = i;}
	
	/** @brief Set the order of the shapefunction associated with the connect */
	void SetOrder(int order, int64_t index = -1) {
#ifdef PZDEBUG
        if(order < 0 || order > 255)
        {
            DebugStop();
        }
//        int64_t crit[] = {5004,5005,4978,5170,5516,5515,5544,5555,4997,4973,4966,4965,5520,5519,5546,5556,5523,5549,5557,5558,4962,5527,5552,5559,5560,5561,5562};
//        static std::set<int64_t> critical;
//        static int once = 1;
//        if (once)
//        {
//            once = 0;
//            for (int i=0; i<27; i++) {
//                critical.insert(crit[i]);
//            }
//        }
//        if (critical.find(index) != critical.end()) {
//            std::cout << "Connect index " << index << " set to order " << order << std::endl;
//            if(order == 5)
//            {
//                DebugStop();
//            }
//        }
#endif
		fCompose.fOrder = order;
	}
    
    /** @brief Set the number of state variables */
    void SetNState(int nstate)
    {
        if(nstate<0 || nstate > 255)
        {
            DebugStop();
        }
        fCompose.fNState = nstate;
    }
    /** @brief Set the number of shape functions associated with the connect */
    void SetNShape(int nshape)
    {
#ifdef PZDEBUG
        if(nshape < 0 )
        {
            DebugStop();
        }
#endif
		fNShape = nshape;
    }
	
	/** @brief Access function to return the order associated with the connect */
	unsigned char Order() const {
		return fCompose.fOrder;
	}
    
    /** @brief Access method to return the indication whether the connect is condensed or not
     */
    bool IsCondensed() const
    {
        return fCompose.fIsCondensed;
    }
    
    /** @brief Access method to return the indication whether the connect is associated with a pressure lagrange multiplier
     */
    unsigned char LagrangeMultiplier() const
    {
        return fCompose.fLagrangeMultiplier;
    }
    
    /** @brief Set the connect as a pressure connect or not */
    void SetLagrangeMultiplier(unsigned char mult)
    {
        fCompose.fLagrangeMultiplier = mult;
    }

    /** @brief Set the connect as a condensed connect or not */
    void SetCondensed(bool flag)
    {
        fCompose.fIsCondensed = flag;
    }

	/** @brief Print the information for the connect element */
	/**
	 * The mesh argument allows the object to identify the number of variables
	 * associated with it and the solution
	 */
	void Print(const TPZCompMesh &mesh, std::ostream & out = std::cout);
	
	/** @brief Also print the center point of the side associated to the connect */
	void Print(TPZCompMesh &mesh, TPZVec<REAL> &cp, std::ostream & out = std::cout);
	
	/** @brief Initialize with zero fNElConnected */
	void ResetElConnected() { fNElConnected = 0; }
	/** @brief Increment fNElConnected */
	void IncrementElConnected() { fNElConnected++; }
	/** @brief Decrement fNElConnected */
	void DecrementElConnected() { fNElConnected--; }
	/** @brief Returns fNElConnected */
	int NElConnected() const { return fNElConnected; }
	
	/**
	 * @brief Add dependency between connects
	 * @param myindex [in] index of this connect
	 * @param dependindex [in] index of the connect this will depend upon
	 * @param depmat [in] dependency matrix which defines the relation between the connects
	 * @param ipos, jpos, isize, jsize are parameters which define the submatrix within depmat which is to be used
	 */
  template<class TVar>
	TPZDepend<TVar> *AddDependency(int64_t myindex, int64_t dependindex,TPZFMatrix<TVar> &depmat,int64_t ipos,int64_t jpos, int isize, int jsize);
	
	/**
	 * @brief Remove dependency between connects if exist
	 * @param myindex [in] index of this connect
	 * @param dependindex [in] index of the connect was exist dependency
	 */
	void RemoveDepend(int64_t myindex, int64_t dependindex);
	
	/** @brief Deletes all dependency information */
	void RemoveDepend();

    /// the maximum number of recursive dependencies
	int DependencyDepth(TPZCompMesh &mesh);
    
    /// size of the dependency list
    int NumDepend() const;
    
	/** @brief Returns whether exist dependecy information */
	int HasDependency()const { return fDependList != 0; }
	
	int CheckDependency(int nshape, TPZCompMesh *mesh, int nstate);
	
	TPZDependBase *FirstDepend() { return fDependList; }
	
	/** @brief Adds itself and the connects from which it depends to the list*/
	void AddToList(int64_t myindex, TPZCompMesh &mesh, TPZStack<int64_t> &connectlist);
	
	/** @brief Adds itself and the connects from which it depends to the list*/
	void AddToList(int64_t myindex, TPZCompMesh &mesh, std::set<int64_t> &connectlist);
	
	void SetDependenceOrder(int64_t myindex, TPZCompMesh &mesh, int CurrentOrder,TPZVec<int64_t> &connectlist,TPZVec<int> &DependenceOrder);
	
	void ExpandShape(int64_t cind, TPZVec<int64_t> &connectlist, TPZVec<int> &blocksize, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	/** @brief Saves the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
	
	/**
	 * @brief Copy a connect data structure from an original connect to a new connect mapping their indexes
	 * @param orig original connect to be copied
	 * @param gl2lcIdx global to local indexes map
	 */
	void CopyFrom(TPZConnect &orig,std::map<int64_t,int64_t> & gl2lcIdx);
	
	/**
	 * @brief Builds the list of all connectivities related to ConnectIndex including the
	 * connects pointed to by dependent connects
	 * @note Note : this method does not reset the stack to zero. The calling
	 * method should do this
	 * @param index [in] index of the current connect
	 * @param indepconnectlist [out] set which contains the indices of independent connects
	 * @param depconnectlist [out] set which contains the indices of dependent connects
	 * @param mesh [in]
	 */
	void BuildConnectList(int64_t index, std::set<int64_t> &indepconnectlist, std::set<int64_t> &depconnectlist, TPZCompMesh &mesh);
	
	/**
	 * @brief Builds the list of all connectivities related to ConnectIndex including the
	 * connects pointed to by dependent connects
	 * @note Note : this method does not reset the stack to zero. The calling
	 * method should do this
	 * @param connectlist [out] stack to receive the list
	 * @param ConnectIndex [in]
	 * @param mesh [in]
	 */
	static void BuildConnectList(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &ConnectIndex, TPZCompMesh &mesh);
	
	/**
	 * @brief Builds the list of all connectivities related to ConnectIndex including the
	 * connects pointed to by dependent connects
	 * @note Note : this method does not reset the set to zero. The calling
	 * method should do this
	 * @param connectlist [out] set to receive the list
	 * @param additional [in] connects which should be added along with their dependencies
	 * @param mesh [in]
	 */
	static void BuildConnectList(std::set<int64_t> &connectlist, std::set<int64_t> &additional, TPZCompMesh &mesh);
	
	/**
	 * @brief This method builds the vector DependenceOrder which indicates in which
	 * order constrained nodes need to be processed
	 */
	/** connectlist need to be computed by BuildConnectList */
	static void BuildDependencyOrder(TPZVec<int64_t> &connectlist, TPZVec<int> &DependenceOrder, TPZCompMesh &mesh);
	
};

/** @ingroup CompElement */
/** @brief Overload operator << to write node connect data */
inline std::ostream & operator<<(std::ostream &out,TPZConnect &con)
{
	out << "seq num: " << con.SequenceNumber()
	<< " nel con: " << con.NElConnected()
	<< " order: " << con.Order()
    << " is condensed " << con.IsCondensed()
	<< " hasdepend: " << con.HasDependency();
	return out;
}

#endif


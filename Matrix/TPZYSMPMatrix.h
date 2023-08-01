/**
 * @file
 * @brief Contains the TPZFYsmpMatrix class which implements a non symmetric sparse matrix. \n
 * Purpose: Defines operations on non-symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */

#ifndef YSMPMATH
#define YSMPMATH

#include "pz_config.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

template<class TVar>
class TPZVerySparseMatrix;
/**
 * @brief Implements a non symmetric sparse matrix (Yale Sparse Matrix Storage). \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Defines operations on general sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */
template<class TVar>
class TPZFYsmpMatrix : public TPZMatrix<TVar> {
	    
	public :
	
	/** @brief An auxiliary structure to hold the data of the subset \n of equations used to multiply in a multi-threaded environment */
	/**
	 In future versions this structure should be defined in a derived class
	 */
	struct TPZMThread {
		const TPZFYsmpMatrix<TVar> *target;
		int64_t fFirsteq;
		int64_t fLasteq;
		const TPZFMatrix<TVar> *fX;
		TPZFMatrix<TVar> *fZ;
		TVar fAlpha;
		int fOpt;
	};
	
private:
	
	static void * ExecuteMT(void *entrydata);
	
public:
    
  TPZFYsmpMatrix();
  TPZFYsmpMatrix(const int64_t rows,const int64_t cols );
	TPZFYsmpMatrix(const TPZFYsmpMatrix<TVar>&) = default;
  TPZFYsmpMatrix(TPZFYsmpMatrix<TVar>&&) = default;
  
	TPZFYsmpMatrix(const TPZVerySparseMatrix<TVar> &cp);
	
	TPZFYsmpMatrix &operator=(const TPZFYsmpMatrix<TVar> &copy) = default;
  TPZFYsmpMatrix &operator=(TPZFYsmpMatrix<TVar> &&copy) = default;
	
	TPZFYsmpMatrix &operator=(const TPZVerySparseMatrix<TVar> &cp);
  inline TPZFYsmpMatrix<TVar>*NewMatrix() const override {return new TPZFYsmpMatrix<TVar>{};}
	CLONEDEF(TPZFYsmpMatrix)
	
	virtual ~TPZFYsmpMatrix();	

  /** @brief Creates a copy from another TPZFYsmpMatrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override
  {                                                           
    auto *from = dynamic_cast<const TPZFYsmpMatrix<TVar> *>(mat);                
    if (from) {                                               
      *this = *from;                                          
    }                                                         
    else                                                      
      {                                                       
        PZError<<__PRETTY_FUNCTION__;                         
        PZError<<"\nERROR: Called with incompatible type\n."; 
        PZError<<"Aborting...\n";                             
        DebugStop();                                          
      }                                                       
  }

  /** @brief Fill matrix storage with randomic values */
  /** This method use GetVal and PutVal which are implemented by each type matrices */
  void AutoFill(int64_t nrow, int64_t ncol, SymProp symmetric) override;
    
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const TVar GetVal(const int64_t row,const int64_t col ) const override;
	
	int64_t NumTerms()
	{
		return fIA[this->Rows()];
	}
	
	int PutVal(const int64_t row, const int64_t col, const TVar &Value) override;
  /**@name Arithmetics */
  /** @{ */
  TPZFYsmpMatrix<TVar> operator+(const TPZFYsmpMatrix<TVar> & A) const;
  TPZFYsmpMatrix<TVar> operator-(const TPZFYsmpMatrix<TVar> & A) const;
  TPZFYsmpMatrix<TVar> operator*(const TVar alpha) const;
  TPZFYsmpMatrix<TVar> &operator+=(const TPZFYsmpMatrix<TVar> &A );
  TPZFYsmpMatrix<TVar> &operator-=(const TPZFYsmpMatrix<TVar> &A );
  TPZMatrix<TVar> &operator*=(const TVar val) override;
  
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0., const int opt = 0) const override;
	
	virtual void MultAddMT(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						   const TVar alpha=1.,const TVar beta = 0., const int opt = 0);

  TVar RowTimesVector(const int row, const TPZFMatrix<TVar> &v) const override;
	/** @} */
	// this method is not working! default version does work
	// virtual int GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
	// 				   const int64_t colSize, TPZFMatrix<TVar> & A ) const override;
	


  
  
	/** @brief Pass the data to the class. */
	virtual void SetData( int64_t *IA, int64_t *JA, TVar *A );
    
    /** @brief Pass the data to the class. */
    virtual void SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A );
  virtual void SetData(TPZVec<int64_t> &&IA, TPZVec<int64_t> &&JA, TPZVec<TVar> &&A);
  /** @brief Get the data from the class*/
  virtual void GetData(TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A);
  /** @brief Get the mem location of data from the class*/
  virtual void GetData(int64_t* &IA, int64_t* &JA, TVar* &A);

  //! Creates a sparse matrix containing only the eqs given in indices
  void GetSubSparseMatrix(const TPZVec<int64_t> &indices, TPZVec<int64_t> &ia,
                          TPZVec<int64_t> &ja, TPZVec<TVar> &aa);
	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const override;
	
	/**
	 * @name Solvers
	 * @brief Linear system solvers. \n
	 */
	 /** For symmetric decompositions lower triangular matrix is used. \n
	 * Solves a system A*X = B returning X in B
	 */  
	//@{
	/**
	 * @brief Solves the linear system using Jacobi method. \n
	 * @param numiterations The number of interations for the process.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param scratch Available manipulation area on memory.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveJacobi(int64_t & numiterations, const TPZFMatrix<TVar> & F, TPZFMatrix<TVar> & result,
							 TPZFMatrix<TVar> * residual, TPZFMatrix<TVar> & scratch, REAL & tol, const int FromCurrent = 0)  override;
	
	void SolveSOR(int64_t &numiterations, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
				  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,
				  const REAL overrelax, REAL &tol,
				  const int FromCurrent = 0,const int direction = 1 )  override;
	// @}
	
	/**
	 * @brief Add a contribution of a stiffness matrix
	 * putting it on destination indexes position
	 */
	virtual void AddKelOld(
						   TPZFMatrix<TVar> & elmat //! Member stiffness matrix beeing added
						   , TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
						   );    
	
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex) override;
	
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex) override;
	
	void MultiplyDummy(TPZFYsmpMatrix<TVar> & B, TPZFYsmpMatrix<TVar> & Res);

    void GetRowIndices(const int64_t i, TPZVec<int64_t>& indices) const override;
	virtual int Zero() override;
	
    /// this is a class that doesn't implement direct decompostion
        /** @brief decompose the system of equations acording to the decomposition
         * scheme */
        virtual int Decompose(const DecomposeType dt) override {
            std::cout << "ERROR! Trying to Decompose base class TPZFYsmpMatrix. To do so, you should use one of its child classes." << std::endl;
            DebugStop();
            return 0;
        }
        /**
         * @brief Solves the linear system using Direct methods
         * @param F The right hand side of the system and where the solution is stored.
         * @param dt Indicates type of decomposition
         */
        virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) override
        {
            std::cout << "ERROR! Trying to SolveDirect on base class TPZFYsmpMatrix. To do so, you should use one of its child classes." << std::endl;
            DebugStop();
            return 0;
        }
        virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const override{
            std::cout << "ERROR! Trying to SolveDirect on base class TPZFYsmpMatrix. To do so, you should use one of its child classes." << std::endl;
            DebugStop();
            return 0;
        }

	
	
  int ClassId() const override;

  /**
   * @brief Unpacks the object structure from a stream of bytes
   * @param buf The buffer containing the object in a packed form
   * @param context
   */
  void Read(TPZStream &buf, void *context) override;

  /**
   * @brief Packs the object structure in a stream of bytes
   * @param buf Buffer which will receive the bytes
   * @param withclassid
   */
  void Write(TPZStream &buf, int withclassid) const override;


	
	void ComputeDiagonal();

private:
	/*
	 * @brief Perform row update of the sparse matrix
	 */
	void RowLUUpdate(int64_t sourcerow, int64_t destrow);
	
protected:
  void CheckTypeCompatibility(const TPZMatrix<TVar>*aPtr,
                              const TPZMatrix<TVar>*bPtr) const override;
  inline TVar *&Elem() override
  {
    return fA.begin();
  }
  inline const TVar *Elem() const override
  {
    return fA.begin();
  }
  inline int64_t Size() const override
  {
    return fA.size();
  }
	TPZVec<int64_t>  fIA;
	TPZVec<int64_t>  fJA;
	TPZVec<TVar> fA;
	
	TPZVec<TVar> fDiag;
	
	int   fSymmetric;
	
protected:
	
	/**
	 * @brief Implements a initialization method for the sparse structure. It sets the initial value for the fIA and fJA.
	 */ 
	/**
	 * -fIA will contain the initial positions for all the equations
	 * -fJA will contain (-1) on all its positions
	 * -fA will contain 0 on all its value 
	 */
	void InitializeData();
};


template<class TVar>
inline void TPZFYsmpMatrix<TVar>::SetData( int64_t *IA, int64_t *JA, TVar *A ) {
	// Pass the data to the class.
    int nel = this->Rows()+1;
    fIA.resize(nel);
//    memccpy(&fIA[0], IA, nel, sizeof(int64_t));
    memcpy(&fIA[0], IA, nel*sizeof(int64_t));
    int64_t nval = fIA[nel-1];
    fJA.resize(nval);
//    memccpy(&fJA[0], JA, nval, sizeof(int64_t));
    memcpy(&fJA[0], JA, nval*sizeof(int64_t));
    fA.resize(nval);
//    memccpy(&fA[0], A, nval, sizeof(TVar));
    memcpy(&fA[0], A, nval*sizeof(TVar));
	ComputeDiagonal();
}

/** @brief Pass the data to the class. */
template<class TVar>
inline void TPZFYsmpMatrix<TVar>::SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A ){
    
    if (IA.size() != this->Rows() + 1 ) {
        DebugStop();
    }
    
    if (JA.size() != IA[this->Rows()]) {
        DebugStop();
    }
    
    fIA = IA;
    fJA = JA;
    fA = A;
    ComputeDiagonal();
}

template<class TVar>
inline void TPZFYsmpMatrix<TVar>::GetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A ){
    IA = fIA;
    JA = fJA;
    A = fA;
}

template<class TVar>
inline void TPZFYsmpMatrix<TVar>::GetData( int64_t* &IA, int64_t* &JA, TVar* &A ){
    IA = fIA.begin();
    JA = fJA.begin();
    A = fA.begin();
}

template<class TVar>
inline void TPZFYsmpMatrix<TVar>::SetData( TPZVec<int64_t> &&IA, TPZVec<int64_t> &&JA, TPZVec<TVar> &&A ){
  if (IA.size() != this->Rows() + 1 ) {
    DebugStop();
  }
    
  if (JA.size() != IA[this->Rows()]) {
    DebugStop();
  }
  fIA = std::move(IA);
  fJA = std::move(JA);
  fA = std::move(A);
}

#endif

/**
 * @file
 * @brief Contains TPZMatrix<TVar>class, root matrix class.
 */

#ifndef _TMATRIXHH_
#define _TMATRIXHH_

#include "pzbasematrix.h"
#include "pzreal.h"


template<class TVar>
class TPZVec;
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZFMatrixRef;
template<class TVar>
class TPZAutoPointer;

class TPZStream;
/** @brief To create clone matrix */
#define CLONEDEF(A) virtual TPZMatrix<TVar>*Clone() const override { return new A(*this); }


class TPZSolver;

/** @brief Root matrix class (abstract). \ref matrix "Matrix" */
/** Abstract class TPZMatrix<TVar>which defines interface of derived matrix classes. */
template<class TVar=STATE>
class TPZMatrix: public TPZBaseMatrix

{
public:
    /// Allows to use TPZMatrix::Type as synonimous to TVar
  using Type = TVar;
	/** @brief Default constructor */
	TPZMatrix() = default;
	/**@brief Copy constructor */
	TPZMatrix<TVar>(const TPZMatrix<TVar>&cp) : TPZRegisterClassId(&TPZMatrix<TVar>::ClassId), TPZBaseMatrix(cp)
	{
	}
    
    /**
     * @brief Constructor
     * @param row Number of rows
     * @param col Number of cols
     */
    inline  TPZMatrix<TVar>(const int64_t row,const int64_t col ) : TPZRegisterClassId(&TPZMatrix<TVar>::ClassId),
        TPZBaseMatrix(row,col)
    {
        
    }

  /** @brief Move constructor */
  TPZMatrix<TVar>(TPZMatrix<TVar> &&cp) = default;
	/** @brief Simple destructor */
	virtual ~TPZMatrix() = default;

  /** @brief Copy assignment operator*/
  TPZMatrix<TVar> &operator=(const TPZMatrix<TVar> &);
  /** @brief Move assignment operator*/
  TPZMatrix<TVar> &operator=(TPZMatrix<TVar> &&);
  /** @brief To create a matrix of the same type */
  virtual TPZMatrix<TVar>*NewMatrix() const = 0;
  /** @brief To create clone matrix */
  virtual TPZMatrix<TVar>*Clone() const = 0;

	/**
	 * @brief Returns the approximate size of the memory footprint (amount
	 * of memory required to store this object).
	 */
    int64_t MemoryFootprint() const override{
	  std::cout << __PRETTY_FUNCTION__ 
		    << ": Please, implement me! (class = " << ClassId() 
	            << std::endl;
	  //<< ") (class name = " << ClassName() << ")" << std::endl;
	  //DebugStop();
	  return 0;
	}

  /** @{ */
  /** @brief Reference to a nx1 TPZFMatrix associated with the contiguous memory area 
   of the matrix.
  This method is useful for arithmetic operations. Derived class should implement
  methods `Size()` and `Elem()`*/
  TPZFMatrixRef<TVar> Storage();
  const TPZFMatrix<TVar> Storage() const;
  /** @} */
  template<class TVar2>
  void CopyFromDiffPrecision(TPZMatrix<TVar2> &copy)
  {
    fDecomposed = copy.IsDecomposed();
    fDefPositive = copy.IsDefPositive();
    fRow = copy.Rows();
    fCol = copy.Cols();

  }

  /** @brief Creates a copy from a given matrix of arbitrary storage format. 
      Every implementation should check for type compatibility */
  virtual void CopyFrom(const TPZMatrix<TVar> *mat) = 0;
  
	/** @brief Fill matrix storage with randomic values */
	/** This method use GetVal and PutVal which are implemented by each type matrices */
	void AutoFill(int64_t nrow, int64_t ncol, SymProp sym) override;
	
	/** @brief Checks if current matrix value is symmetric */
	SymProp VerifySymmetry(REAL tol = ZeroTolerance()) const override;
	
	/**
     * @brief Put value with bounds checking
     * @param row Row number.
     * @param col Column number.
     * @param value Value being put.
	 */
	virtual int Put(const int64_t row,const int64_t col,const TVar & value );
	/**
	 * @brief Get value with bound checking
     * @param row Row number.
     * @param col Column number.
	 */
	virtual const TVar Get(const int64_t row,const int64_t col ) const;
	
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
     * @param row Row number.
	 * @param col Column number.
	 */
	TVar &operator() (const int64_t row,const int64_t col );

    /**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
     * @param row Row number.
	 * @param col Column number.
	 */
	TVar operator() (const int64_t row,const int64_t col ) const;
  /**
   * @brief The operators check on the bounds if the DEBUG variable is defined
   * @param row Row number.
   * @param col Column number.
   */
  TVar &at(const std::pair<int64_t,int64_t> &rowcol )
  {
    return operator()(rowcol.first,rowcol.second);
  }
  const TVar at(const std::pair<int64_t,int64_t> &rowcol ) const
  {
    return Get(rowcol.first,rowcol.second);
  }
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
   * @param row Row number.
   * @param col Column number.
	 */
	virtual TVar &s(const int64_t row, const int64_t col);
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
     * @param row Row number.
	 */
	TVar &operator()(const int64_t row);
	
	/** @brief Put values without bounds checking \n
	 *  This method is faster than "Put" if DEBUG is defined.
	 */
	virtual int PutVal(const int64_t /*row*/,const int64_t /*col*/,const TVar & val )
    {
        if(val != ((TVar)(0.))) DebugStop();
        return 0;
    }
	/** @brief Get matrix entry without bound checking.
	 */
    virtual const TVar GetVal(const int64_t /*row*/, const int64_t /*col*/ ) const;
	
	/** @name Algebraic
	 *  @brief Implements algebraic operations with matrices
	 * @{
	 */
	
	/**
	 * @brief It mutiplies itself by TPZMatrix<TVar>A putting the result in res
	 * @param A TPZMatrix<TVar>object to multiplied to
	 * @param res TPZFMatrix<TVar>containing the result
	 * @param opt Indicates if is Transpose or not
	 */
	virtual void Multiply(const TPZFMatrix<TVar>& A,TPZFMatrix<TVar>& res, int opt = 0) const;
  /**
	 * @brief It mutiplies itself by a scalar alpha putting the result in res
	 * @param alpha scalar to be multiplied with
	 * @param res TPZFMatrix<TVar>containing the result
	 */
  void MultiplyByScalar(const TVar alpha,TPZMatrix<TVar>& res) const;
	/**
	 * @brief It adds itself to TPZMatrix<TVar>A putting the result in res
	 * @param A TPZMatrix<TVar>to added to current matrix
	 * @param res Contains the result
	 */
	void Add(const TPZMatrix<TVar>& A,TPZMatrix<TVar>& res) const;
  /** @brief It substracts A from storing the result in result */
	void Subtract(const TPZMatrix<TVar>& A,TPZMatrix<TVar>& result) const;

  virtual TPZMatrix<TVar> &operator*=(const TVar val);

  TPZFMatrix<TVar> operator*(const TPZFMatrix<TVar> &B ) const;
  //! Unifies input parameters checks of MultAdd
  void MultAddChecks(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
                     const TVar alpha=1., const TVar beta = 0., const int opt = 0) const;
	/**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 */
	virtual void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
						 const TVar alpha=1., const TVar beta = 0., const int opt = 0) const;

    virtual TVar RowTimesVector(const int row, const TPZFMatrix<TVar> &v) const;

  /**
   * @brief It computes this += alpha*(A * B), where A or B can be transposed.
   * @param i Is the row of (this) where the first element of the matrices product should be added 
   * @param j Is the column of (this) where the first element of the matrices product should be added 
   * @param A Is A on the above operation
   * @param transpA 0: A not transpose, 1: A tranpose, everything else: A conj transpose
   * @param B Is B on the above operation
   * @param transpB 0: B not transpose, 1: B tranpose, everything else: B conj transpose
   * @param alpha Is alpha on the above operation
   */
	virtual void AddContribution(int64_t i, int64_t j, const TPZFMatrix<TVar> & A, int transpA, const TPZFMatrix<TVar>& B, 
						 		 int transpB, const TVar alpha = 1.0);
	
	/** @brief Computes res = rhs - this * x */
	virtual void Residual(const TPZFMatrix<TVar>& x,const TPZFMatrix<TVar>& rhs, TPZFMatrix<TVar>& res ) ;
	
	/** @brief Converts the matrix in an identity matrix*/
	virtual void Identity();
    /** @brief Converts the matrix in a diagonal matrix*/
    virtual void Diagonal(TVar val);
	/** @brief It makes *T the transpose of current matrix. */
	virtual void Transpose(TPZMatrix<TVar>*const T, bool conj = false) const;
	
	/**
	 * @brief It makes Inv =[this].
	 * IMPORTANT OBSERVATION --> The original matrix (calling object) no is more equal. 
	 * It containts the some decomposition (LU or Cholesky or ...)
	 */
	int Inverse(TPZFMatrix<TVar>&Inv, DecomposeType dec);
	
	/**
	 * @brief Computes the matrix norm of this
	 * @param p interpolation order
	 * @param numiter is used by 2-norm calculation in the SolveEigenvaluesJacobi method required to compute the maximum eigenvalue
	 * @param tol - same of numiter
	 * @note It is available p-norm = 1, 2 and infinity.\n
	 * p=1 is the maximum absolute column sum norm \n
	 * p=2 is the spectral norm wich is the square root of the maximum eigenvalue of Tranpose[this].this \n
	 * p=infinity is the maximum absolute row sum norm - p infinity is implemented with p = 0 \n
	 * These operations are defined on the website of the Mathematica software: \n
	 * http://mathworld.wolfram.com/MatrixNorm.html \n
	 * Be careful when choosing 2-norm. It has a high computational cost.
	 */
	TVar MatrixNorm(int p, int64_t numiter = 2000000, REAL tol = 1.e-10) const; // -<
	
	/** @brief Computes the matrix condition number of this */
	/**
	 * It is available p-norm = 1, 2 and infinity.
	 * p=1 is the maximum absolute column sum norm
	 * p=2 is the spectral norm wich is the square root of the maximum eigenvalue of Tranpose[this].this
	 * p=infinity is the maximum absolute row sum norm - p infinity is implemented with p = 0
	 * These operations are defined on the website of the Mathematica software:
	 * http://mathworld.wolfram.com/MatrixNorm.html
	 * All norms require the computation of the inverse matrix.
	 * It has a high computational cost and a high memory requirement.
	 */
	TVar ConditionNumber(int p, int64_t numiter = 2000000, REAL tol = 1.e-10);
	
	/** @} */
	
	/**
	 * @name Generic
	 * @brief Performs generic operations on matrices
	 * @{
	 */
	
	/** @brief Input operation*/
	virtual void Input(std::istream & in = std::cin );
	
	/** @brief Input operation*/
	template <class TT> friend std::istream & operator >> (std::istream& in,TPZMatrix<TT>& A) ;
	
    virtual void Print(std::ostream &out) const 
    {
        std::string name = __PRETTY_FUNCTION__;
        Print(name.c_str(),out);
    }
	/** @brief It prints the matrix data in a MatrixFormat Rows X Cols */
	void Print(const char *name, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const override;
	
	/**
	 * @brief Redimensions a matriz keeping the previous values
	 * @param newRows Specifies the new number of rows in matrix
	 * @param newCols Specifies the new number of Columns in matrix
	 */
	int Resize(const int64_t newRows, const int64_t newCols ) override{
		fRow = newRows;
		fCol = newCols;
		return 0;
	}
	
	/**
	 * @brief Redimensions the matrix reinitializing it with zero
	 * @param newRows Specifies the new number of rows in matrix.
	 * @param newCols Specifies the new number of Columns in matrix.
	 */
	int Redim(const int64_t newRows, const int64_t newCols ) override {
		fRow = newRows;
		fCol = newCols;
		return 0;
	}
	
	/** @brief Zeroes the matrix */
	int Zero() override{
		std::cout << "WARNING! TPZMatrix<TVar>::Zero is called\n";
        DebugStop();
		return 0; }
	
	/** @} */
	
	/**
	 * @name SubMatrices
	 * Operations with SUB MATRICES.
	 * @{
	 */
	
	/**
	 * @brief It puts submatrix Source on actual matrix structure.
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param Source The matrix to be inserted
	 */
	virtual int PutSub( const int64_t sRow, const int64_t sCol, const TPZFMatrix<TVar>& Source );
	
	/**
	 * @brief Gets submatrix storing it on Target.
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param rowSize Specifies the amount of rows from sRow
	 * @param colSize Specifies the amount of columns from sCol
	 * @param Target The matrix to be aquired.
	 */
	virtual int GetSub( const int64_t sRow, const int64_t sCol, const int64_t rowSize,
					   const int64_t colSize, TPZFMatrix<TVar>& Target ) const;
	
	/**
	 * @brief It adds Source matrix on current matrix from position (sRow, sCol)
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param Source The matrix to be added
	 */
	virtual int AddSub(const int64_t sRow, const int64_t sCol, const TPZFMatrix<TVar>& Source );
	
	/**
	 * @brief Inserts a submatrix from current object on matrix *Target with no \n
	 * redimentioning
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param rowSize Specifies the amount of rows from sRow
	 * @param colSize Specifies the amount of columns from sCol
	 * @param pRow Specifies final row on current object
	 * @param pCol Specifies final column on current object.
	 * @param Target The matrix to be inserted.
	 */
	virtual int InsertSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,const int64_t colSize,
						  const int64_t pRow,const int64_t pCol, TPZMatrix<TVar>* Target ) const;
	
	/**
	 * @brief Adds a submatrix from current object in *Target
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param rowSize Specifies the amount of rows from sRow
	 * @param colSize Specifies the amount of columns from sCol
	 * @param pRow Specifies starting row on pA
	 * @param pCol Specifies starting column on pA.
	 * @param pA The matrix to be added.
	 */
	virtual int AddSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
					   const int64_t colSize,const int64_t pRow,const int64_t pCol, TPZMatrix<TVar>* pA ) const;
	
	/** @} */
	
	/**
	 * @brief Add a contribution of a stiffness matrix
	 * @param elmat Element matrix to be contributed
	 * @param destinationindex Contains destine indexes on current matrix
	 */
	virtual  void AddKel(TPZFMatrix<TVar>&elmat, TPZVec<int64_t> &destinationindex);

	/**
	 * @brief Add a contribution of a stiffness matrix
	 * @param elmat Element matrix to be contributed
	 * @param sourceindex Contains source indexes on current matrix
	 * @param destinationindex Contains destine indexes on current matrix
	 */
	virtual  void AddKel(TPZFMatrix<TVar>&elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex);

    /**
     * @brief Similar to Addkel, but uses atomic add instead of regular sum
     */
    virtual  void AddKelAtomic(TPZFMatrix<TVar>&elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex){
        DebugStop();
    };

	/**
	 * @name Inquire
	 * @brief Returns information of the current object
	 * @{
	 */
	
	/** @brief Updates the values of the matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > /* mat*/)    
	{
		std::cout << "TPZMatrix<TVar>::UdateFrom is not implemented\n";
        DebugStop();
	}
	
	/**
	 * @name Solvers
	 * @brief Linear system solvers. \n
	 * @{
	 */
	/** 
	 * For symmetric decompositions lower triangular matrix is used. \n
	 * Solves a system A*X = B returning X in B
	 */
	
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
	virtual void SolveJacobi(int64_t & numiterations, const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
							 TPZFMatrix<TVar>* residual, TPZFMatrix<TVar>& scratch, REAL & tol, const int FromCurrent = 0);
	
	/**
	 * @brief Solves the linear system using Successive Over Relaxation method (Gauss Seidel). \n
	 * @param numiterations The number of interations for the process.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param scratch Available manipulation area on memory.
	 * @param overrelax The over relaxation parameter
	 * @param tol The tolerance value..
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 * @param direction Indicates interaction direction, from first to last (default 1) or from last to first (-1)
	 */
	virtual void SolveSOR(int64_t & numiterations, const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						  TPZFMatrix<TVar>* residual,TPZFMatrix<TVar>& scratch,const REAL overrelax, REAL & tol,
						  const int FromCurrent = 0,const int direction = 1) ;
	/**
	 * @brief Solves the linear system using Symmetric Successive Over Relaxation method (Gauss Seidel). \n
	 * @param numiterations The number of interations for the process.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param scratch Available manipulation area on memory.
	 * @param overrelax The over relaxation parameter
	 * @param tol The tolerance value..
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveSSOR(int64_t & numiterations,const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						   TPZFMatrix<TVar>* residual, TPZFMatrix<TVar>& scratch, const REAL overrelax, REAL & tol,
						   const int FromCurrent = 0) ;
	
	/**
	 * @brief Solves the linear system using Conjugate Gradient method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent.
	 */
	virtual void SolveCG(int64_t & numiterations, TPZSolver & preconditioner,
						 const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						 TPZFMatrix<TVar>* residual, REAL & tol,
						 const int FromCurrent = 0) ;
	/**
	 * @brief Solves the linear system using Bi-Conjugate Gradient method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param tol The tolerance value.
	 */
	virtual void SolveBICG(int64_t & numiterations, TPZSolver & preconditioner,
						   const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						   REAL & tol) ;
	
	/**
	 * @brief Solves the linear system using Bi-Conjugate Gradient stabilized method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent.
	 */
	virtual void SolveBICGStab(int64_t & numiterations, TPZSolver & preconditioner,
							   const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
							   TPZFMatrix<TVar>* residual, REAL & tol,
							   const int FromCurrent = 0) ;
	
	/**
	 * @brief Solves the linear system using Generalized Minimal Residual (GMRES) method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param H The right hand side of the system
	 * @param numvectors The number of vectors involved
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveGMRES(int64_t & numiterations, TPZSolver & preconditioner,
							TPZFMatrix<TVar>& H, int & numvectors,
							const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
							TPZFMatrix<TVar>* residual, REAL & tol,const int FromCurrent) ;
	
	/**
	 * @brief Solves the linear system using IR method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveIR(int64_t & numiterations, TPZSolver & preconditioner,
						 const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						 TPZFMatrix<TVar>* residual, REAL & tol,
						 const int FromCurrent = 0);
	
	/** @brief Transforms this matrix in a diagonal matrix, where the diagonal values are its eigenvalues.
	 * This method is efficient only for small matrices.
	 * @param numiterations The number of interations for the process.
	 * @param tol The tolerance value.
	 * @param Sort diagonal values from big to small
	 * @return Returns true if tolerance is achieved or false otherwise.
	 */
	virtual bool SolveEigenvaluesJacobi(int64_t &numiterations, REAL &tol, TPZVec<TVar> * Sort = 0);
	
	/** @brief Compute Eigenvalues and Eigenvectors of this matrix. \n
	 * This method is efficient only for small matrices.
	 * @param numiterations The number of interations for the process.
	 * @param tol The tolerance value.
	 * @param Eigenvalues ordered from big to small
	 * @param Eigenvectors: each row represent one eigenvector. It is in same order of eigenvalues.
	 * @return Returns true if tolerance is achieved or false otherwise.
	 */
	virtual bool SolveEigensystemJacobi(int64_t &numiterations, REAL & tol, TPZVec<TVar> & Eigenvalues, TPZFMatrix<TVar>& Eigenvectors) const;
	
	/**
	 * @brief Solves the linear system using Direct methods
	 * @param F The right hand side of the system and where the solution is stored.
	 * @param dt Indicates type of decomposition
	 */
	virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) = 0;
    virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const = 0;

    
/** @brief Retorna o valor mais proximo a "val" (exceto valores no intervalo -tol <= val <= +tol) contido no vetor Vec */
	static TVar ReturnNearestValue(TVar val, TPZVec<TVar> &Vec, TVar tol);
	
	
	/**
	 * @name TPZSavable
	 * @brief Methods which would make TPZMatrix<TVar>compliant with TPZSavable
	 * @{
	 */
	
    int ClassId() const override;
	
	/** @} */
	
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSavable *copy, bool override = false) override;
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSavable *copy, bool override = false) const override;
	
	/** @brief Extract the block indicated by the indices from the matrix */
	virtual void GetSub(const TPZVec<int64_t> &indices,TPZFMatrix<TVar>&block) const;
	
	/** @brief Compare values of this to B, with a precision tolerance tol. */
    bool CompareValues(TPZMatrix<TVar>&M, TVar tol);
	
protected:
  /** @brief Checks compatibility of matrices before Add/Subtract operations*/
  virtual void CheckTypeCompatibility(const TPZMatrix<TVar>*A,
                                      const TPZMatrix<TVar>*B)const;
  /** @brief Number of entries storaged in the Matrix*/
  virtual int64_t Size() const {DebugStop(); return -1;}
  /** @{ */
  /** @brief Pointer to the beginning of the storage of the matrix*/
  virtual TVar* &Elem() { DebugStop(); static TVar* t{nullptr}; return t; }
  virtual const TVar* Elem() const { DebugStop(); return nullptr;}
  /** @} */
	/**
	 * @brief Is an auxiliar method used by MultiplyAdd
	 * @see MultAdd
	 */
	void PrepareZ(const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,const TVar beta,const int opt) const;
	
	

    TVar GetRandomVal() const;
protected:
	/**
	 * @brief Returns error messages and aborts executioin
	 * @param msg First message.
	 * @param msg2 Second message.
	 */
	static int Error(const char *msg ,const char *msg2 = 0);
	/** @brief It clears data structure. */
	virtual int Clear() { return 0; }
	
	/** @brief Swaps contents of a in b and b in a */
	static void Swap(int64_t *a, int64_t *b);
    /** @brief value considered as zero*/
	static TVar gZero;
};



/** @brief Initializing value to static variable */
template <class TVar>
TVar TPZMatrix<TVar>::gZero = TVar(0);
/** @} */

/** @brief Overload << operator to print entries of the matrix ***/
template<class TVar>
std::ostream & operator<<(std::ostream& out, const TPZMatrix<TVar> & A);

/******** Inline ********/


template<class TVar>
inline void TPZMatrix<TVar>::Residual(const TPZFMatrix<TVar>& x,const TPZFMatrix<TVar>& rhs, TPZFMatrix<TVar>& res )  {
	DebugStop();
}

template<>
inline void
TPZMatrix<float>::Residual(const TPZFMatrix<float> &x, const TPZFMatrix<float> &rhs, TPZFMatrix<float> &res) {
    MultAdd(x, rhs, res, ((float) -1.0), ((float) 1.0));
}

template<>
inline void TPZMatrix<double>::Residual(const TPZFMatrix<double>& x,const TPZFMatrix<double>& rhs, TPZFMatrix<double>& res )  {
	MultAdd( x, rhs, res, ((double)-1.0), ((double)1.0) );
}

template<>
inline void
TPZMatrix<std::complex<float>>::Residual(const TPZFMatrix<std::complex<float>> &x, const TPZFMatrix<std::complex<float>> &rhs,
                                    TPZFMatrix<std::complex<float>> &res) {
    MultAdd(x, rhs, res, ((std::complex<float>) -1.0), ((std::complex<float>) 1.0));
}

template<>
inline void
TPZMatrix<std::complex<double>>::Residual(const TPZFMatrix<std::complex<double>> &x, const TPZFMatrix<std::complex<double>> &rhs,
                                     TPZFMatrix<std::complex<double>> &res) {
    MultAdd(x, rhs, res, ((std::complex<double>) -1.0), ((std::complex<double>) 1.0));
}

/***********/
/*** Put ***/

template<class TVar>
inline int TPZMatrix<TVar>::Put(const int64_t row,const int64_t col,const TVar & value ) {
	// bound checking
#ifdef PZDEBUG
	if ( row >= Rows() || col >= Cols() || row <0 || col < 0 ) {
		std::stringstream sout;
		sout << "TPZMatrix<TVar>::Put" << " Index out of range row = " << row << " col = " << col << " Rows() " << Rows() << " Cols() " << Cols() << std::endl;
		Error(sout.str().c_str());
	}
#endif
	
	return( PutVal( row, col, value ) );
}



/***********/
/*** Get ***/

template<class TVar>
inline const TVar TPZMatrix<TVar>::Get(const int64_t row, const int64_t col ) const {
	// bound checking
#ifdef PZDEBUG
	if ( (row >= Rows()) || (col >= Cols()) || row <0 || col <0 ) {
		Error("TPZMatrix::Get", "Index out of range");
        DebugStop();
	}
#endif
	return( GetVal( row, col ) );
}

template<class TVar>
inline TVar &TPZMatrix<TVar>::operator()(const int64_t row, const int64_t col) {
	// bound checking
#ifndef PZNODEBUG
	if ( (row >= Rows()) || (col >= Cols()) || row <0 || col<0 ) {
		Error("TPZMatrix<TVar>::Operator()","Index out of range");
        DebugStop();
	}
#endif
	return s(row,col);
}

template<class TVar>
inline TVar TPZMatrix<TVar>::operator()(const int64_t row, const int64_t col) const{
	// bound checking
#ifndef PZNODEBUG
	if ( (row >= Rows()) || (col >= Cols()) || row <0 || col<0 ) {
		Error("TPZMatrix<TVar>::Operator()","Index out of range");
        DebugStop();
	}
#endif
	return GetVal(row,col);
}

template<class TVar>
inline TVar &TPZMatrix<TVar>::s(const int64_t row, const int64_t col) {
    DebugStop();
    throw "TPZMatrix<TVar>::s not implemented\n";
    
}

template<class TVar>
inline TVar &TPZMatrix<TVar>::operator()(const int64_t row) {
	return operator()(row,0);
}
//***Solve LU ***/



template<class TVar>
inline void
TPZMatrix<TVar>::Swap( int64_t *a, int64_t *b )
{
	int64_t aux = *a;
	*a = *b;
	*b = aux;
}
#endif

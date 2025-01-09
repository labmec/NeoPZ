/**
 * @file
 * @brief Contains TPZBaseMatrix class, root matrix class.
 * Since this class is not a template class, its will only
 * define methods that take no other matrices.
 */

#ifndef TPZBASEMATRIX_H
#define TPZBASEMATRIX_H

#include "TPZSavable.h"
#include "pzreal.h"

/** @brief Defines output format */
enum MatrixOutputFormat {EFormatted, EInputFormat, EMathematicaInput, EMatlabNonZeros, EMatrixMarket, ECSV, EFixedColumn};

class TPZStream;
/** \addtogroup matrix
 * @{
 */
/**
 * @enum DecomposeType
 * @brief Defines decomposition type for any matrix classes
 * @param ENoDecompose Not decomposed
 * @param ELU Decomposed using LU method
 * @param ECholesky Decomposed using Cholesky method
 * @param ELDLt Decomposed using LDLt method
 */
enum DecomposeType { ENoDecompose, ELU, ELUPivot, ECholesky, ELDLt };

/**
 * @enum SymProp
 * @brief Defined symmetry property of a given matrix
 * @param NonSym Non symmetric
 * @param Sym Symmetric
 * @param Herm Hermitian 
 */
enum class SymProp{NonSym, Sym, Herm};

template<class TVar>
class TPZVec;

/** @brief Root matrix class (abstract). \ref matrix "Matrix"
 * Abstract class TPZBaseMatrix which is agnostic with respect to
 * the arithmetic type being used. */
class TPZBaseMatrix : public TPZSavable{
public:
  /** @brief Default constructor */
  TPZBaseMatrix() : TPZRegisterClassId(&TPZBaseMatrix::ClassId) {
    fDecomposed = ENoDecompose;
    fDefPositive = 0;
    fRow = 0;
    fCol = 0;
  }
    TPZBaseMatrix(int64_t row, int64_t col) {
      fDecomposed = ENoDecompose;
      fDefPositive = 0;
      fRow = row;
      fCol = col;
    }
  /** @brief Copy constructor*/
  TPZBaseMatrix(const TPZBaseMatrix &cp) = default;
  /** @brief Move constructor*/
  TPZBaseMatrix(TPZBaseMatrix &&cp) = default;
  /** @brief Simple destructor */
  virtual ~TPZBaseMatrix() = default;

  /** @brief Copy assignment operator*/
  TPZBaseMatrix &operator=(const TPZBaseMatrix &) = default;
  /** @brief Move assignment operator*/
  TPZBaseMatrix &operator=(TPZBaseMatrix &&) = default;

  /**
   * @brief Returns the approximate size of the memory footprint (amount
   * of memory required to store this object).
   */
  virtual int64_t MemoryFootprint() const = 0;

  /** @brief Fill matrix storage with randomic values */
  virtual void AutoFill(int64_t nrow, int64_t ncol, SymProp sym) = 0;
    
  /** @brief Returns number of rows */
  inline int64_t Rows() const { return fRow; }
  /** @brief Returns number of cols */
  inline int64_t Cols() const { return fCol; }

  /** @brief Returns the dimension of the matrix if the matrix is square.*/
  /** If the matrix is not square, returns an error */
  inline virtual int64_t Dim() const {
    if (IsSquare())
      return Rows();
    PZError << "Cannot call TPZBaseMatrix::Dim() for a non-square ";
    PZError << "matrix! Aborting...";
    DebugStop();
    return 0; // compiler wont complain
  }

  /**
   * @brief Redimensions a matriz keeping the previous values
   * @param newRows Specifies the new number of rows in matrix
   * @param newCols Specifies the new number of Columns in matrix
   */
  virtual int Resize(const int64_t newRows, const int64_t newCols) = 0;
  /**
   * @brief Redimensions the matrix reinitializing it with zero
   * @param newRows Specifies the new number of rows in matrix.
   * @param newCols Specifies the new number of Columns in matrix.
   */
  virtual int Redim(const int64_t newRows, const int64_t newCols) = 0;

  /** @brief Zeroes the matrix */
  virtual int Zero() = 0;

  virtual void GetRowIndices(const int64_t i, TPZVec<int64_t> &indices) const;

  /** @brief Gets symmetry property of current matrix
      @note This flag is set by the user. Use VerifySymmetry for actually checking values.*/
  SymProp GetSymmetry() const{ return fSymProp; }

  /** @brief Performs a check to ensure if current matrix value is symmetric */
  virtual SymProp VerifySymmetry(REAL tol) const = 0;

  /** @brief Sets symmetry property of current matrix.
      If matrix is not square and symmetry property is not SymProp::NonSym, this function throws an error.
      Matrix classes that store only lower/upper triangular portion of the matrix will overload this method
      as non-symmetric formats are not allowed.
      @note Should only be called if the user is sure about the symmetry property. Can lead to inconsistent states.*/
  virtual void SetSymmetry(SymProp sp);
  /** @brief Checks if current matrix is square */
  inline int IsSquare() const {
    return fRow == fCol; }

  /** @brief Simetrizes copies upper plan to the lower plan, making its data
   * simetric */
//  virtual void Simetrize() = 0;

  /** @brief Sets current matrix as definite positive */
  void SetDefPositive(bool v) {
    fDefPositive = v;
  }
  
  /** @brief Checks if current matrix is definite positive */
  virtual bool IsDefPositive() const {
    return fDefPositive;
  }
  /** @brief Checks if current matrix is already decomposed */
  DecomposeType IsDecomposed() const {
    return fDecomposed; }



  /** @brief Sets current matrix to decomposed state */
  virtual void SetIsDecomposed(DecomposeType val) {
    fDecomposed = val; }


  /** @brief decompose the system of equations acording to the decomposition
   * scheme */
    virtual int Decompose(const DecomposeType dt) {
      std::cout << "TPZBaseMatrix::Decompose is not implemented\n";
      DebugStop();
      return 0;
    }

  /** @brief It prints the matrix data in a MatrixFormat Rows X Cols */
  virtual void Print(const char *name, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const = 0;
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

protected:
  /** @brief Number of rows in matrix */
  int64_t fRow;
  /** @brief Number of cols in matrix */
  int64_t fCol;
  /** @brief Decomposition type used to decompose the current matrix */
  DecomposeType fDecomposed;
  /** @brief Symmetry property of the matrix*/
  SymProp fSymProp{SymProp::NonSym};
  /** @brief Definite Posistiveness of current matrix */
  bool fDefPositive;
};

//! Convert enum to string
constexpr auto SymPropName(SymProp sp){
  switch(sp){
  case SymProp::NonSym: return "NonSym";
  case SymProp::Sym: return "Sym";
  case SymProp::Herm: return "Herm";
  }
  //avoids compiler warning
  unreachable();
}

//! Convert enum to string
constexpr auto DecomposeTypeName(DecomposeType dec){
  switch(dec){
  case DecomposeType::ELU: return "ELU";
  case DecomposeType::ELDLt: return "ELDLt";
  case DecomposeType::ELUPivot: return "ELUPivot";
  case DecomposeType::ECholesky: return "ECholesky";
  case DecomposeType::ENoDecompose: return "ENoDecompose";
  }
  //avoids compiler warning
  unreachable();
}

/** @} */

#endif

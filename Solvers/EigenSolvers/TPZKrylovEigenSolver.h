#ifndef TPZKRYLOVEIGENSOLVER_H
#define TPZKRYLOVEIGENSOLVER_H
#include "TPZEigenSolver.h"
#include "TPZSpectralTransform.h"

//! Sorting method for calculated eigenvalues
enum class TPZEigenSort{
  EAbsAscending,/*!< Ascending magnitude*/
  EAbsDescending,/*!< Descending magnitude*/
  ERealAscending,/*!< Ascending real part*/
  ERealDescending,/*!< Descending real part*/
  EImagAscending,/*!< Ascending imaginary part*/
  EImagDescending/*!< Descending imaginary part*/
};

/** @brief Solvers for eigenvalue problems using Krylov methods.
    The eigenvalue problem is solved in the projected Krylov subspace 
    obtained by an Arnoldi iteration. 
    See TPZSpectralTransform for possible spectral transformations*/
template<class TVar>
class TPZKrylovEigenSolver : public TPZEigenSolver<TVar>
{
public:
  using TPZEigenSolver<TVar>::TPZEigenSolver;

  TPZKrylovEigenSolver<TVar> * Clone() const override
  {return new TPZKrylovEigenSolver<TVar>(*this);}

  /** @name BasicUsage */
  /** @{*/  
  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns 1 if executed correctly
   */
  int SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @return Returns 1 if executed correctly
   */
  int SolveEigenProblem(TPZVec<CTVar> &w) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @param[out] eigenVectors Stores the correspondent eigenvectors
   * @return Returns 1 if executed correctly
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                   TPZFMatrix<CTVar> &eigenVectors) override;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and does not calculates the eigenvectors
   * @param[out] w Eigenvalues in ascending magnitude order
   * @return Returns 1 if executed correctly
   */
  int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w) override;

  //! Sets spectral transformation to be used (it creates an internal copy)
  void SetSpectralTransform(TPZSpectralTransform<TVar> &s);
  //! Gets spectral transformation to be used
  inline TPZAutoPointer<TPZSpectralTransform<TVar>> SpectralTransform();
  /** 
      @brief Sets number of Eigenpairs to calculate. 
      @note If the dimension of the Krylov subspace is insufficient, 
      it will be adjusted to `n`.
  */
  inline void SetNEigenpairs(int n);
  //! Gets number of Eigenpairs to calculate
  inline int NEigenpairs() const;
  /** @brief Sets the dimension of the Krylov subspace. 
      @note If not set, defaults to `10*nev`, where `nev` is the number 
      of sought eigenvalues.*/
  inline void SetKrylovDim(int d);
  //! Gets the dimension of the Krylov subspace
  inline int KrylovDim() const;
  //! Sets tolerance for norm of the Krylov vectors
  inline void SetTolerance(RTVar s);
  //! Gets tolerances
  inline RTVar Tolerance() const;
  /** @brief Sets the first vector to be used in the Krylov subspace
      @note The input vector will be normalized. If not set, a random 
      vector will be generated.
  */
  inline void SetKrylovInitialVector(TPZFMatrix<TVar> q);
  //! Gets the first vector used in the Krylov subspace
  inline TPZFMatrix<TVar> KrylovInitialVector() const;
  /** @brief Decides criterium for sorting the obtained eigenvalues. 
      @note By default it is set to TPZEigenSort::EAbsAscending .*/
  inline void SetEigenSorting(TPZEigenSort ord);
  //! Returns criterium for sorting the obtained eigenvalues
  inline TPZEigenSort EigenSorting() const;
  

  /** @brief Performs the Arnoldi iteration on a given matrix. 
      This iteration creates an orthonormal basis of the Krylov subspace of A. 
      The first vector of the Krylov basis can be set through 
      SetKrylovInitialVector method. The dimension of the Krylov subspace 
      can be set throught SetKrylovDim
      @param[in] A Matrix in which the algorithm is performed
      @param[out] Q Arnoldi vectors of A
      @param[out] H Representation of A in basis Q
      @return returns true if succeeded, false otherwise
  */
  [[nodiscard]] bool ArnoldiIteration(const TPZMatrix<TVar> &A,
                                      TPZVec<TPZAutoPointer<TPZFMatrix<TVar>>> &Q,
                                      TPZFMatrix<TVar> &H);
  /** @}*/
protected:
  //! Number of Eigenpairs to calculate
  int fNEigenpairs{1};
  //! Dimension of the Krylov subspace to calculate
  int fKrylovDim{-1};
  //! Initial vector to be used to create Krylov subspace
  TPZFMatrix<TVar> fKrylovVector;
  //! Tolerance
  RTVar fTolerance{std::numeric_limits<RTVar>::epsilon()};
  //! Spectral Transformation
  TPZAutoPointer<TPZSpectralTransform<TVar>> fST;
  //! Sorting order of the eigenvalues
  TPZEigenSort fEigenSort{TPZEigenSort::EAbsAscending};
  //! Implementation of Solve methods
  int SolveImpl(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors, bool computeVectors);
};

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetSpectralTransform(TPZSpectralTransform<TVar> &s)
{
  fST = s.Clone();
}

template<class TVar>
TPZAutoPointer<TPZSpectralTransform<TVar> >
TPZKrylovEigenSolver<TVar>::SpectralTransform()
{
  return fST;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetNEigenpairs(int n)
{
  if(n < 1) n = 1;
  fNEigenpairs = n;
  if(n > fKrylovDim){
    fKrylovDim = 10 * n;
    std::cout<< "Adjusted Krylov dim to "<< fKrylovDim<<std::endl;
  }
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::NEigenpairs() const
{
  return fNEigenpairs;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetKrylovDim(int k)
{
  if(k<2){
    k = 2;
    std::cout<< "Adjusted Krylov dim to "<< k<<std::endl;
  }
  fKrylovDim = k;
}

template<class TVar>
int TPZKrylovEigenSolver<TVar>::KrylovDim() const
{
  return fKrylovDim;
}


template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetTolerance(RTVar tol)
{
  fTolerance = tol;
}
  
template<class TVar>
RTVar TPZKrylovEigenSolver<TVar>::Tolerance() const
{
  return fTolerance;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetKrylovInitialVector(TPZFMatrix<TVar> q)
{
  fKrylovVector = q;
}

template<class TVar>
TPZEigenSort TPZKrylovEigenSolver<TVar>::EigenSorting() const
{
  return fEigenSort;
}

template<class TVar>
void TPZKrylovEigenSolver<TVar>::SetEigenSorting(TPZEigenSort ord)
{
  fEigenSort = ord;
}
template<class TVar>
TPZFMatrix<TVar> TPZKrylovEigenSolver<TVar>::KrylovInitialVector() const
{
  return fKrylovVector;
}

extern template class TPZKrylovEigenSolver<float>;
extern template class TPZKrylovEigenSolver<double>;
extern template class TPZKrylovEigenSolver<std::complex<float>>;
extern template class TPZKrylovEigenSolver<std::complex<double>>;
#endif
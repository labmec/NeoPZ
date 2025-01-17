/**
 * @file
 * @brief Contains the implementation of the TPZEigenSparseMatrix methods.
 */

#include "TPZEigenSparseMatrix.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#include <memory.h>
#include <string>
#include <map>
#include <thread>
#include <vector>
#include "tpzverysparsematrix.h"
#include "pzstack.h"

using namespace std;

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************
// ****************************************************************************
//
// Constructor
//
// ****************************************************************************


template<class TVar>
TPZEigenSparseMatrix<TVar>::TPZEigenSparseMatrix() : TPZRegisterClassId(&TPZEigenSparseMatrix::ClassId),
TPZSYsmpMatrix<TVar>()
{
}





template<class TVar>
TPZEigenSparseMatrix<TVar>::TPZEigenSparseMatrix(const int64_t rows,const int64_t cols ) :
TPZRegisterClassId(&TPZEigenSparseMatrix::ClassId),TPZSYsmpMatrix<TVar>(rows,cols) {
#ifdef CONSTRUCTOR
	cerr << "TPZEigenSparseMatrix(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZEigenSparseMatrix<TVar>::~TPZEigenSparseMatrix() {
	// Deletes everything associated with a TPZEigenSparseMatrix
#ifdef DESTRUCTOR
	cerr << "~TPZEigenSparseMatrix()\n";
#endif
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************


template<class TVar>
void TPZEigenSparseMatrix<TVar>::CheckTypeCompatibility(const TPZMatrix<TVar>*A, const TPZMatrix<TVar>*B)const
{
  auto incompatSparse = [](){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: incompatible matrices\n.Aborting...\n";
    DebugStop();
  };
	auto aPtr = dynamic_cast<const TPZEigenSparseMatrix<TVar>*>(A);
  auto bPtr = dynamic_cast<const TPZEigenSparseMatrix<TVar>*>(B);
  if(!aPtr || !bPtr){
    incompatSparse();
  }
	bool check{false};
	const auto nIA = aPtr->fIA.size();
	for(auto i = 0; i < nIA; i++){
		check = check || aPtr->fIA[i] != bPtr->fIA[i];
	}

	const auto nJA = aPtr->fJA.size();
	for(auto i = 0; i < nJA; i++){
		check = check || aPtr->fJA[i] != bPtr->fJA[i];
	}
	if(check) incompatSparse();
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
TPZEigenSparseMatrix<TVar> TPZEigenSparseMatrix<TVar>::operator+(const TPZEigenSparseMatrix<TVar>&mat) const
{
	CheckTypeCompatibility(this,&mat);
	auto res(*this);
    res += mat;
	return res;
}
template<class TVar>
TPZEigenSparseMatrix<TVar> TPZEigenSparseMatrix<TVar>::operator-(const TPZEigenSparseMatrix<TVar>&mat) const
{
	CheckTypeCompatibility(this,&mat);
	auto res(*this);
    res -= mat;
	return res;
}

template<class TVar>
TPZEigenSparseMatrix<TVar> TPZEigenSparseMatrix<TVar>::operator*(const TVar alpha) const
{
	auto res(*this);
	for(auto &el : res.fA) el *= alpha;
	return res;
}

template<class TVar>
TPZEigenSparseMatrix<TVar> &TPZEigenSparseMatrix<TVar>::operator+=(const TPZEigenSparseMatrix<TVar> &A )
{
	TPZEigenSparseMatrix<TVar> res((*this)+A);
	*this = res;
	return *this;
}
template<class TVar>
TPZEigenSparseMatrix<TVar> &TPZEigenSparseMatrix<TVar>::operator-=(const TPZEigenSparseMatrix<TVar> &A )
{
	TPZEigenSparseMatrix<TVar> res((*this)-A);
	*this = res;
	return *this;
}
template<class TVar>
TPZMatrix<TVar> &TPZEigenSparseMatrix<TVar>::operator*=(const TVar val)
{
	TPZEigenSparseMatrix<TVar> res((*this)*val);
	*this = res;
	return *this;
}


// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZEigenSparseMatrix<TVar>::Print(const char *title, ostream &out ,const MatrixOutputFormat form) const {
	// Print the matrix along with a identification title
    TPZSYsmpMatrix<TVar>::Print(title,out,form);
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************




static int64_t  FindCol(int64_t *colf, int64_t *coll, int64_t col)
{
	if(col == *colf) return 0;
	int64_t *begin = colf;
	int64_t *end = coll;
	while (begin != end)
	{
		int64_t dist = (end-begin)/2;
		int64_t *mid = begin+dist;
		if(*mid == col) return (mid-colf);
		else if(*mid > col) end=mid;
		else begin = mid;
	}
	return -1;
}

/// this is a class that doesn't implement direct decompostion
/** @brief decompose the system of equations acording to the decomposition
 * scheme */
template<class TVar>
int TPZEigenSparseMatrix<TVar>::Decompose(const DecomposeType dt)
{
//    Eigen::setNbThreads(16);
    std::cout << "\n=====> Starting factorization with Eigen..." << std::endl;
    const int n = Eigen::nbThreads( );
    std::cout << "Using n = " << n << " threads" << std::endl;
    
    if(!fEigenMatrix){
        fEigenMatrix = new EigenSparse(this->fRow,this->fCol,this->fJA.size(),&this->fIA[0],&this->fJA[0],&this->fA[0]);
    }
    if(dt == ECholesky) {
#ifndef USEACCELERATESUPPORT
        fCholesky = new EigenCholesky(*fEigenMatrix);        
        fCholesky->analyzePattern(*fEigenMatrix);
        fCholesky->factorize(*fEigenMatrix);
#endif
    }
    else if(dt == ELDLt) {
#if defined(MACOSX) && defined(USEACCELERATESUPPORT)
        fLDLT = new EigenLDLT(*fEigenMatrix);
#else
        fLDLT = new EigenLDLT(*fEigenMatrix);
#endif
        fLDLT->analyzePattern(*fEigenMatrix);
        fLDLT->factorize(*fEigenMatrix);
    }
    else if(dt == ELU) {
#ifndef USEACCELERATESUPPORT
        fLU = new EigenLU(*fEigenMatrix);
        fLU->analyzePattern(*fEigenMatrix);
        fLU->factorize(*fEigenMatrix);
#endif
    }
    else {
        DebugStop();
    }
    this->SetIsDecomposed(dt);
    return 0;
}
/**
 * @brief Solves the linear system using Direct methods
 * @param F The right hand side of the system and where the solution is stored.
 * @param dt Indicates type of decomposition
 */
template<class TVar>
int TPZEigenSparseMatrix<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt)
{
    Eigen::Map<Eigen::Matrix<TVar,Eigen::Dynamic,1>  > b(&F(0,0),F.Rows());
    if(dt == ECholesky) {
        if(fLDLT || fLU) DebugStop();
        if(!fCholesky) Decompose(dt);
        b = fCholesky->solve(b);
    }
    else if(dt == ELDLt) {
        if(fCholesky || fLU) DebugStop();
        if(!fLDLT) Decompose(dt);
        using namespace Eigen;
        //SimplicialLDLT<SparseMatrix<double,0,int64_t> > *solver = fLDLT;
        b = fLDLT->solve(b);
    }
    else if(dt == ELU) {
        if(fLDLT || fCholesky) DebugStop();
        if(!fLU) Decompose(dt);
        b = fLU->solve(b);
    }
    else{
        DebugStop();
    }
    return 0;
}

template<class TVar>
int TPZEigenSparseMatrix<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const
{
    // const Eigen::Map<Eigen::Matrix<TVar,Eigen::Dynamic,1>  > b(&F(0,0),F.Rows());
    // Eigen::Map<Eigen::Matrix<TVar,Eigen::Dynamic,1>  > bres(&F(0,0),F.Rows());

    // if(fCholesky) {
    //     if(dt != ECholesky) DebugStop();
    //     b = fCholesky->solve(b);
    // }
    // else if(fLDLT) {
    //     if(dt != ELDLt) DebugStop();
    //     b = fLDLT->solve(b);
    // }
    // else if(fLU) {
    //     if(dt != ELU) DebugStop();
    //     b = fLU->solve(b);
    // }
    // else{
    //     DebugStop();
    // }
    DebugStop();
    return 0;
}

template<class TVar>
void TPZEigenSparseMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                         const TVar alpha,const TVar beta, const int opt) const {
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        std::cout << "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" ;
        return;
    }
    if(!IsZero(beta) && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        std::cout << "TPZFMatrix::MultAdd matrix y with incompatible dimensions>";
        return;
    }
#endif
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0) {
        z.Zero();
        return;
    }
    
    DebugStop(); // Please implement me
    
    // z = beta * y + alpha * opt(this)*x
    
//    Eigen::Map<Eigen::Matrix<TVar,Eigen::Dynamic,1>  > b(&F(0,0),F.Rows());

}


template<class TVar>
int TPZEigenSparseMatrix<TVar>::ClassId() const{
    return Hash("TPZEigenSparseMatrix") ^ TPZSYsmpMatrix<TVar>::ClassId() << 1;
}

template class TPZEigenSparseMatrix<double>;
//template class TPZEigenSparseMatrix<float>;

#ifndef USEACCELERATESUPPORT
template class TPZEigenSparseMatrix<long double>;
template class TPZEigenSparseMatrix<std::complex<long double>>;
template class TPZEigenSparseMatrix<std::complex<double>>;
template class TPZEigenSparseMatrix<std::complex<float>>;
#endif

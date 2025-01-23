/**
 * @file
 * @brief Contains the implementation of the TPZFMatrix<>methods.
 */


#include "pzfmatrix.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <complex>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include "pzerror.h"
#include "pzaxestools.h"
#include "pzextractval.h"
#include "pzlog.h"
#include "pzmatrix.h"
#include "TPZSavable.h"
#include "pzvec.h"
#include "tpzverysparsematrix.h"
#include "TPZParallelUtils.h"
#include "tfad.h"
#include "fad.h"

class TPZStream;

#ifdef PZDEBUG
#define DEBUG2
#endif

#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzfmatrix");
static TPZLogger loggerCheck("pz.checkconsistency");
#endif

#ifdef USING_LAPACK
/** Maybe some calls just need BLAS. We need to check it. */
#include "TPZLapack.h"
#include "TPZLapackEigenSolver.h"
#define BLAS_MULT
#endif


//#define IsZero( a )  ( fabs(a) < 1.e-20)


using namespace std;

/********************/
/*** Constructors ***/

template <class TVar>
TPZFMatrix<TVar>::TPZFMatrix(const TPZMatrix<TVar> &mat) : 
TPZRegisterClassId(&TPZFMatrix::ClassId),
TPZMatrix<TVar>(mat), fElem(0),fGiven(0),fSize(0) {
    if(this->fRow*this->fCol) {
        
        fElem = new TVar[this->fRow*this->fCol];
        
        int64_t i,j;
        for(j=0; j<this->fCol; j++) {
            for(i=0; i<this->fRow; i++) {
                this->PutVal(i, j, mat.GetVal(i,j));
            }
        }
    }
}


/********************************/
/*** Constructor( TPZFMatrix<TVar> & ) ***/

template<class TVar>
TPZFMatrix<TVar>::TPZFMatrix(const TPZFMatrix<TVar> &A)
: TPZRegisterClassId(&TPZFMatrix::ClassId),
TPZMatrix<TVar>( A.fRow, A.fCol ), fElem(0), fGiven(0), fSize(0) {
    int64_t size = this->fRow * this->fCol;
    if(!size) return;
    fElem = new TVar[ size ] ;
#ifdef PZDEBUG2
    if ( size && fElem == NULL ) Error( "Constructor <memory allocation error>." );
#endif
    // Copia a matriz
    TVar * src = A.fElem;
    TVar * p = fElem;
    for(int64_t i = 0; i< size; i++) p[i]=src[i];
//    memcpy((void *)(p),(void *)(src),(size_t)size*sizeof(TVar));
    this->fSymProp = A.GetSymmetry();
}


template<class TVar>
TPZFMatrix<TVar>::TPZFMatrix(TPZFMatrix<TVar> &&A)
    : TPZMatrix<TVar>(A), fElem(0),
      fGiven(0),fSize(0),
      fPivot(0), fWork(0)
{
    // A has preallocated memory, we cant copy the pointer
    if(A.fElem == A.fGiven)
    {
        int64_t size = A.fRow * A.fCol;
        fElem = new TVar[size] ;
        for(int64_t i = 0; i< size; i++) fElem[i]=A.fElem[i];
//        memcpy((void *)(fElem),(void *)(A.fElem),(size_t)size*sizeof(TVar));
    }
    else {
        fElem = A.fElem;
        A.fElem=nullptr;
    }
    fPivot = A.fPivot;
    fWork = A.fWork;
}

/********************************/
/*** Constructor( TPZVerySparseMatrix<TVar> & ) ***/

template<class TVar>
TPZFMatrix<TVar>::TPZFMatrix(TPZVerySparseMatrix <TVar> const & A)
: TPZRegisterClassId(&TPZFMatrix::ClassId),
TPZMatrix<TVar>( A.Rows(), A.Cols() ), fElem(0), fGiven(0), fSize(0) {
    
    int64_t size = this->fRow * this->fCol;
    if(!size) return;
    fElem = new TVar[ size ] ;
    
#ifdef PZDEBUG2
    if ( size && fElem == NULL ) Error( "Constructor <memory allocation error>." );
#endif
    
    typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator it = A.MapBegin();
    typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator end = A.MapEnd();
    
    TVar * p = fElem;
    memset((void *)p, 0, (size_t)size*sizeof(TVar));
    
    for (; it != end; it++) {
        const std::pair<int64_t, int64_t>& key = it->first;
        PutVal(key.first, key.second, it->second);
    }
    
}

template<class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator=(const TPZFMatrix<TVar> &A ) {
    if(this == &A) return *this;
    // amount of storage needed
    int64_t size = A.fRow * A.fCol;
    
    TVar * newElem = fElem;
    if(fSize < size && size != this->fRow*this->fCol) {
        newElem = new TVar[size];
    } else if (fSize >= size) {
        newElem = fGiven;
    }
    
    if ( newElem == NULL && size > 0) Error( "Operator= <memory allocation error>." );
    if (fElem && fElem != newElem && fElem != fGiven) delete[]( fElem );
    this->fRow  = A.fRow;
    this->fCol  = A.fCol;
    fElem = newElem;
    
    // Copia a matriz
    for(int64_t i = 0; i<size; i++) fElem[i] = A.fElem[i];
//    memcpy((void *)(fElem),(void *)(A.fElem),(size_t)size*sizeof(TVar));
    
    TPZMatrix<TVar>::operator=(A);
    
    
    return *this;
}

template<class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator=(TPZFMatrix<TVar> &&A ) {
    if(this == &A) return *this;
    // if A points to preallocated memory, then this code is wrong
    if(A.fGiven && A.fElem == A.fGiven) {
        // amount of storage needed
        int64_t size = A.fRow * A.fCol;
        
        TVar * newElem = 0;
        if(size == this->fRow*this->fCol) {
            newElem = fElem;
        } else if(fSize < size) {
            // not enough preallocated space, allocate dynamically
            newElem = new TVar[size];
        } else {
            newElem = fGiven;
        }
        if ( newElem == NULL && size > 0) Error( "Operator= <memory allocation error>." );
        if (fElem && fElem != newElem && fElem != fGiven) delete[]( fElem );
        fElem = newElem;
        // Copia a matriz
        for(int64_t i = 0; i<size; i++) fElem[i] = A.fElem[i];
    } else {
        // A does not have or does not use preallocated memory
        if (fElem && fElem != A.fElem && fElem != fGiven) delete[]( fElem );
        fElem=A.fElem;
        fPivot = A.fPivot;
        fWork = A.fWork;
        
        A.fElem = nullptr;
    }
    TPZMatrix<TVar>::operator=(A);
    return *this;
}

/******** Operacoes com matrizes FULL  ********/

/******************/
/*** Operator = ***/

template< class TVar >
TPZFMatrix<TVar>& TPZFMatrix<TVar>::operator= (const std::initializer_list<TVar>& list) {
	Resize(list.size(), 1);

	auto it = list.begin();
	auto it_end = list.end();
	TVar* aux = fElem;
	for (; it != it_end; it++, aux++)
		*aux = *it;

	return *this;
}

template<class TVar>
void TPZFMatrix<TVar>::CopyFrom(const TPZMatrix<TVar> *mat)
{
    const auto r = mat->Rows();
    const auto c = mat->Cols();
    this->Resize(r,c);
    for(auto i = 0 ;i < r;i++){
        for(auto j = 0; j < c; j++){
            this->PutVal(i,j,mat->GetVal(i,j));
        }
    }
}

template< class TVar >
void TPZFMatrix<TVar>::InitializeEqualFromList (const std::initializer_list< std::initializer_list<TVar> >& list){
    size_t n_row = list.size();
    size_t n_col = 1;

    auto row_it = list.begin();
    auto row_it_end = list.end();

#ifdef PZDEBUG
    bool col_n_found = false;
    for (auto it = row_it; it != row_it_end; it++) {
        if (!col_n_found) {
            n_col = it->size();
            col_n_found = true;
        }
        else {
            if (n_col != it->size())
                this->Error("TPZFMatrix constructor: inconsistent number of columns in initializer list");
        }
    }
#else
    n_col = row_it->size();
#endif
    this->Resize(n_row, n_col);

    for (uint32_t row_n = 0; row_it != row_it_end; row_it++, row_n++) {
        auto col_it = row_it->begin();
        auto col_it_end = row_it->end();
        for (uint32_t col_n = 0; col_it != col_it_end; col_it++, col_n++) {
            this->fElem[col_n * this->fRow + row_n] = *col_it;
        }
    }
}

template< class TVar >
TPZFMatrix<TVar>& TPZFMatrix<TVar>::operator= (const std::initializer_list< std::initializer_list<TVar> >& list) {
	this->InitializeEqualFromList(list);

	return *this;
}


template <class TVar>
void TPZFMatrix<TVar>::AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = rhs.Rows();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
            operator()(destination[i],j) += rhs(i,j);
        }
    }
}

template<class TVar>
void TPZFMatrix<TVar>::AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols() && source.NElements()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = source.NElements();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
            operator()(destination[i],j) += rhs(source[i],j);
        }
    }
}

template<>
void TPZFMatrix<double>::AddFel(TPZFMatrix<double> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols() && source.NElements()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = source.NElements();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
          pzutils::AtomicAdd(operator()(destination[i],j),rhs(source[i],j));
        }
    }
}

template<>
void TPZFMatrix<float>::AddFel(TPZFMatrix<float> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols() && source.NElements()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = source.NElements();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
          pzutils::AtomicAdd(operator()(destination[i],j),rhs(source[i],j));
        }
    }
}

template<class TVar>
void TPZFMatrix<TVar>::AddFelNonAtomic(TPZFMatrix<TVar> &rhs,TPZVec<int64_t> &source, TPZVec<int64_t> &destination) {
    if(rhs.Cols() != this->Cols() && source.NElements()) {
        PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
        DebugStop();
        return;
    }
    int64_t ncol = this->Cols();
    int64_t nrow = source.NElements();
    int64_t i,j;
    for(j=0; j<ncol; j++) {
        for(i=0; i<nrow; i++) {
            operator()(destination[i],j) += rhs(source[i],j);
        }
    }
}


/*******************************/
/*** Operator+( TPZFMatrix>& ) ***/

template <class TVar>
TPZFMatrix<TVar>
TPZFMatrix<TVar>::operator+(const TPZFMatrix<TVar> &A ) const {
    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
        Error( "Operator+ <matrixs with different dimensions>" );
    
    auto res(*this);
    TVar * pa = A.fElem;
    const auto size = this->Rows() * this->Cols();
    TVar * pr = res.fElem, *plast = pr+size;
    
    while(pr < plast) *pr++ +=  (*pa++);
    
    return res;
}

/*******************************/
/*** Operator-( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar>
TPZFMatrix<TVar>::operator-(const TPZFMatrix<TVar> &A ) const {

    if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
        Error( "Operator- <matrixs with different dimensions>" );
    
    auto res (*this);
    TVar * pa = A.fElem;
    const auto size = this->Rows()*this->Cols();
    TVar * pr = res.fElem, *prlast =pr+size;
    
    while(pr < prlast) *pr++ -= (*pa++);
    return res;
}

template <class TVar>
void TPZFMatrix<TVar>::GramSchmidt(TPZFMatrix<TVar> &Orthog, TPZFMatrix<TVar> &TransfToOrthog)
{
    if constexpr (std::is_same<TVar, TFad<6, REAL>>::value ||
                  std::is_same<TVar, Fad<float>>::value ||
                  std::is_same<TVar, Fad<double>>::value ||
                  std::is_same<TVar, Fad<long double>>::value||
                  std::is_same<TVar,int>::value||
                  std::is_same<TVar,TPZFlopCounter>::value) {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" not implemented for this type\n";
        PZError<<"Aborting...";
        DebugStop();
    }
#ifdef PZ_LOG2
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        Print("GrSchmidt Entrada",sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    double scale = 1.;
    for(int64_t j = 0; j < this->Cols(); j++)
    {
        double norm = 0.;
        for(int64_t i = 0; i < this->Rows(); i++)
        {
            norm += fabs(TPZExtractVal::val(this->GetVal(i,j)*this->GetVal(i,j)));
        }
        norm = sqrt(norm);
        if(norm > 1.e-10)
        {
            if((1.)/norm > scale) scale = (1.)/norm;
        }
    }
    
    this->operator *=( scale );
    
    int64_t QTDcomp = this->Rows();
    int64_t QTDvec = this->Cols();
    Orthog.Resize(QTDcomp,QTDvec);
    Orthog.Zero();
    /// Making a copy of *this (Ortog = *this)
    for(int64_t r = 0; r < QTDcomp; r++)
    {
        for(int64_t c = 0; c < QTDvec; c++)
        {
            Orthog(r,c) = GetVal(r,c);
        }
    }
    
#ifdef PZDEBUG
    int check = 0;
    for(int64_t c = 0; c < QTDvec; c++)
    {
        TVar summ = 0.;
        for(int64_t r = 0; r < QTDcomp; r++)
        {
            summ += fabs(GetVal(r,c));
        }
        if(fabs(summ) < 0.00001)
        {
            std::stringstream sout;
            sout << "Null Vector on Gram-Schmidt Method! Col = " << c << "\n";
            LOGPZ_ERROR(logger,sout.str())
            check = 1;
        }
    }
#endif
    
    TVar dotUp, dotDown;
    for(int64_t c = 1; c < QTDvec; c++)
    {
        for(int64_t stop = 0; stop < c; stop++)
        {
            dotUp = 0.;
            dotDown = 0.;
            for(int64_t r = 0; r < QTDcomp; r++)
            {
                dotUp += GetVal(r,c)*Orthog(r,stop);
                dotDown += Orthog(r,stop)*Orthog(r,stop);
            }
            if(fabs(dotDown) < 1.E-8)
            {
#ifdef PZDEBUG
                if(check == 0)
                {
                    std::stringstream sout;
                    sout << "Parallel Vectors on Gram-Schmidt Method! Col = " << stop << "\n";
                    LOGPZ_ERROR(logger,sout.str())
                }
#endif
                
                for(int64_t r = 0; r < QTDcomp; r++)
                {
                    Orthog(r,stop) = 0.;
                }
            }
            else
            {
#ifdef PZ_LOG2
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "dotdown = " << dotDown << " dotup = " << dotUp;
                    LOGPZ_DEBUG(logger,sout.str())
                }
#endif
                
                for(int64_t r = 0; r < QTDcomp; r++)
                {
                    Orthog(r,c) -= dotUp*Orthog(r,stop)/dotDown;
                }
            }
        }
    }
    for(int64_t c = 0; c < QTDvec; c++)
    {
        dotUp = 0.;
        for(int64_t r = 0; r < QTDcomp; r++)
        {
            dotUp += Orthog(r,c)*Orthog(r,c);
        }
        if(fabs(dotUp) > 1.e-8)
        {
            for(int64_t r = 0; r < QTDcomp; r++)
            {
                Orthog(r,c) = Orthog(r,c)/sqrt(dotUp);
            }
        }
        else {
#ifdef PZ_LOG
            std::stringstream sout;
            sout << "Linearly dependent columns dotUp = " << dotUp;
            LOGPZ_ERROR(logger,sout.str())
#endif
            
            for(int64_t r = 0; r < QTDcomp; r++)
            {
                Orthog(r,c) = 0.;
            }
        }
        
    }
    Orthog.Multiply(*this,TransfToOrthog,1);
    
    this->operator*= ( 1./scale );
    TransfToOrthog.operator*= ( 1./scale );
    
#ifdef PZ_LOG2
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << endl;
        sout << "Output GS" << endl;
        Orthog.Print("Orthog matrix",sout);
        TransfToOrthog.Print("TransfToOrthog matrix",sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
#ifdef PZDEBUG
    TPZFNMatrix<9, TVar> OrthogT;
    Orthog.Transpose(&OrthogT);
    TPZAxesTools<TVar>::VerifyAxes(OrthogT);
#endif
}


template <class TVar>
void TPZFMatrix<TVar>::DeterminantInverse(TVar &determinant, TPZFMatrix<TVar> &inverse)
{
    TPZFNMatrix<100, TVar> copy(*this);
    inverse.Redim(this->Rows(),this->Rows());
    int64_t r;
    for(r=0; r<this->Rows(); r++) inverse(r,r) = 1.;
    copy.Solve_LU(&inverse);
    if(this->Rows() == 1) {
        determinant = GetVal(0, 0);
    } else if(this->Rows() == 2) {
        determinant = g(0,0)*g(1,1)-g(0,1)*g(1,0);
    } else if(this->Rows() == 3) {
        determinant = 0.;
        for(int i=0; i<3; i++) {
            int j = (i+1)%3;
            int k = (i+2)%3;
            determinant += g(0,i)*g(1,j)*g(2,k);
            determinant -= g(2,i)*g(1,j)*g(0,k);
        }
    } else {
        // not sure if this is correct
        determinant = 1.;
        for(int64_t i=0; i<this->Rows(); i++) determinant *= copy.g(i,i);
        if(copy.fPivot.size()) {
            for(int64_t i=0; i<this->Rows(); i++) if(copy.fPivot[i] != i+1) determinant *= -1.;
        }

    }
}


template <class TVar>
/** @brief Initialize pivot with i = i  */
void TPZFMatrix<TVar>::InitializePivot()
{
    const int nRows = this->Rows();
    fPivot.Resize(nRows);
    for(int64_t i = 0; i < nRows; i++){
        fPivot[i] = i+1; // Fortran based indexing
    }
}

template <class TVar>
void TPZFMatrix<TVar>::MultAdd(const TVar *ptr, int64_t rows, int64_t cols, const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                               const TVar alpha,const TVar beta ,const int opt)
{
    
    
    if ((!opt && cols != x.Rows()) || (opt && rows != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (TVar)0. && ((!opt && rows != y.Rows()) || (opt && cols != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
        return;
    }
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != rows) {
            z.Redim(rows,x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != cols) {
            z.Redim(cols,x.Cols());
        }
    }
    unsigned numeq = opt ? cols : rows;
    int64_t xcols = x.Cols();
    int64_t ic, c;
    if(!(rows*cols)) return;
    for (ic = 0; ic < xcols; ic++) {
        TVar *zp = &z(0,ic), *zlast = zp+numeq;
        if(beta != (TVar)0.) {
            const TVar *yp = &y.g(0,ic);
            if(&z != &y) {
                for (int64_t i = 0; i<numeq; i++) {
                    zp[i] = yp[i];
                }
//                memcpy((void *)zp,(void *)yp,numeq*sizeof(TVar));
            }
            for(int64_t i=0; i< numeq; i++) for(int64_t c=0; c<xcols; c++) z(i,c) *= beta;
        } else {
            while(zp != zlast) {
                *zp = 0.;
                zp ++;
            }
        }
    }
    
    
    for (ic = 0; ic < xcols; ic++) {
        if(!opt) {
            for ( c = 0; c<cols; c++) {
                TVar * zp = &z(0,ic), *zlast = zp+rows;
                const TVar * fp = ptr +rows*c;
                const TVar * xp = &x.g(c,ic);
                while(zp < zlast) {
                    *zp += alpha* *fp++ * *xp;
                    zp ++;
                }
            }
        } else {
            const TVar * fp = ptr;
            TVar *zp = &z(0,ic);
            for (c = 0; c<cols; c++) {
                TVar val = 0.;
                // bug correction philippe 5/2/97
                //					 REAL * xp = &x(0,ic), xlast = xp + numeq;
                const TVar *xp = &x.g(0,ic);
                const TVar *xlast = xp + rows;
                while(xp < xlast) {
                    val += *fp++ * *xp;
                    xp ++;
                }
                *zp += alpha *val;
                zp ++;
            }
        }
    }
    
    
    
    
    
}

#ifdef USING_LAPACK
template<>
void TPZFMatrix<double>::MultAdd(const TPZFMatrix<double> &x,const TPZFMatrix<double> &y, TPZFMatrix<double> &z,
                                 const double alpha,const double beta,const int opt) const {
    
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (double)0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
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
        if (beta != 0) {
            z = y;
            z *= beta;
        }
        return;
    }
    if (beta != (double)0.) {
        z = y;
    }
    if (Rows() == 0 || Cols() == 0 || x.Rows() == 0 || x.Cols() == 0) {
        return;
    }
    if (!opt) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->Rows(), x.Cols(), this->Cols(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    } else {
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->Cols(), x.Cols(), this->Rows(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    }
    
}
template<>
void TPZFMatrix<float>::MultAdd(const TPZFMatrix<float> &x,const TPZFMatrix<float> &y, TPZFMatrix<float> &z,
                                const float alpha,const float beta,const int opt) const {
    
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (float)0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
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
    }
    if (beta != (float)0.) {
        z = y;
    }
    if (Rows() == 0 || Cols() == 0 || x.Rows() == 0 || x.Cols() == 0) {
        return;
    }
    if (!opt) {
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->Rows(), x.Cols(), this->Cols(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    } else {
        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->Cols(), x.Cols(), this->Rows(),
                    alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), beta, z.fElem, z.Rows());
    }
    
}

template<>
void TPZFMatrix<std::complex<double> >::MultAdd(const TPZFMatrix<std::complex<double> > &x,const TPZFMatrix<std::complex<double> > &y, TPZFMatrix<std::complex<double> > &z,
                                                const std::complex<double> alpha,const std::complex<double> beta,const int opt) const {
    
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(abs(beta) != 0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
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
    }
    if (abs(beta) != 0.) {
        z = y;
    }

    if (Rows() == 0 || Cols() == 0 || x.Rows() == 0 || x.Cols() == 0) {
        return;
    }
    
    if (!opt) {
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, this->Rows(), x.Cols(), this->Cols(),
                    &alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), &beta, z.fElem, z.Rows());
    } else {
        cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans, this->Cols(), x.Cols(), this->Rows(),
                    &alpha, this->fElem, this->Rows(), x.fElem, x.Rows(), &beta, z.fElem, z.Rows());
    }
}
#endif // USING_LAPACK

/**
 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
 * @param x Is x on the above operation
 * @param y Is y on the above operation
 * @param z Is z on the above operation
 * @param alpha Is alpha on the above operation
 * @param beta Is beta on the above operation
 * @param opt Indicates if is Transpose or not
 */
template <class TVar>
void TPZFMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                               const TVar alpha,const TVar beta,const int opt) const {
    
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
        return;
    }
    if(beta != (TVar)0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
        return;
    }
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0)
    {
        z.Zero();
    }
    unsigned numeq = opt ? this->Cols() : this->Rows();
    int64_t rows = this->Rows();
    int64_t cols = this->Cols();
    int64_t xcols = x.Cols();
    int64_t ic, c;
    if (numeq)
    {
        for (ic = 0; ic < xcols; ic++) {
            TVar *zp = &z(0,ic), *zlast = zp+numeq;
            if(beta != (TVar)0.) {
                const TVar *yp = &y.g(0,ic);
                if(&z != &y) {
                    for(int64_t i = 0; i<numeq; i++) zp[i]=yp[i];
//                    memcpy((void *)zp,(void *)yp,numeq*sizeof(TVar));
                }
                for(int64_t i=0; i< numeq; i++) z(i,ic) *= beta;
                
            } else {
                while(zp != zlast) {
                    *zp = 0.;
                    zp ++;
                }
            }
        }
    }
    
    if(!(rows*cols)) return;
    
    for (ic = 0; ic < xcols; ic++) {
        if(!opt) {
            for ( c = 0; c<cols; c++) {
                TVar * zp = &z(0,ic), *zlast = zp+rows;
                TVar * fp = fElem +rows*c;
                const TVar * xp = &x.g(c,ic);
                while(zp < zlast) {
                    *zp += alpha* *fp++ * *xp;
                    zp ++;
                }
            }
        } else {
            TVar * fp = fElem,  *zp = &z(0,ic);
            for (c = 0; c<cols; c++) {
                TVar val = 0.;
                // bug correction philippe 5/2/97
                //					 REAL * xp = &x(0,ic), xlast = xp + numeq;
                const TVar *xp = &x.g(0,ic);
                const TVar *xlast = xp + rows;
                while(xp < xlast) {
                    val += *fp++ * *xp;
                    xp ++;
                }
                *zp += alpha *val;
                zp ++;
            }
        }
    }
}

template <class TVar>
void TPZFMatrix<TVar>::AddContribution(int64_t i, int64_t j, const TPZFMatrix<TVar> & A, bool transpA, const TPZFMatrix<TVar>& B, 
						 		       bool transpB, const TVar alpha)
{
    if constexpr (!std::is_floating_point_v<RTVar>)
        Error( "AddContribution must be used with a floating point type variable TVar!\n");

#ifdef USING_LAPACK
    lapack_int nrows = A.Rows();
    lapack_int ncols = B.Cols();
    lapack_int seconddimA = A.Cols();
    lapack_int seconddimB = B.Rows();
    CBLAS_TRANSPOSE transposeA = CblasNoTrans;
    CBLAS_TRANSPOSE transposeB = CblasNoTrans;

    if (transpA)
    {
        nrows = A.Cols();
        seconddimA = A.Rows();
        transposeA = CblasTrans;
    }
    if (transpB)
    {
        ncols = B.Rows();
        seconddimB = B.Cols();
        transposeB = CblasTrans;
    }

#else
    int64_t nrows = A.Rows();
    int64_t ncols = B.Cols();
    int64_t seconddimA = A.Cols();
    int64_t seconddimB = B.Rows();

    if (transpA)
    {
        nrows = A.Cols();
        seconddimA = A.Rows();
    }
    if (transpB)
    {
        ncols = B.Rows();
        seconddimB = B.Cols();
    }
#endif

#ifdef PZDEBUG
    if (seconddimA != seconddimB)
    {
        std::stringstream out;
        out << "AddContribution Matrices A and B have incompatible dimensions. ";
        out << "Second dimension of A: " << seconddimA << ". Second dimension of B: " << seconddimB;
        Error( out.str().c_str() );
        return;
    }
    else if ((i+nrows > this->Rows()) || (j+ncols > this->Cols()))
    {
        std::stringstream out;
        out << "AddContribution trying to add a submatrix out of the full matrix bounds. ";
        out << "Submatrix bigger row: " << i+nrows << ". Full matrix size: " << this->Rows();
        out << "Submatrix bigger column: " << j+ncols << ". Full matrix size: " << this->Cols();
        Error( "AddContribution trying to add a submatrix out of the original matrix bounds" );
        return;
    }
#endif

#ifdef USING_LAPACK
    TVar* ptr = &(*this)(i,j);
    const TVar beta = 1.0;
    if constexpr(std::is_same_v<TVar,double>){
        cblas_dgemm(CblasColMajor, transposeA, transposeB, nrows, ncols, seconddimA, alpha, A.fElem, A.Rows(), B.fElem, B.Rows(), beta, ptr, this->Rows());
        return;
    }
    else if constexpr(std::is_same_v<TVar,float>){
        cblas_sgemm(CblasColMajor, transposeA, transposeB, nrows, ncols, seconddimA, alpha, A.fElem, A.Rows(), B.fElem, B.Rows(), beta, ptr, this->Rows());
        return;
    }
    else if constexpr(std::is_same_v<TVar,std::complex<double>>){
        cblas_zgemm(CblasColMajor, transposeA, transposeB, nrows, ncols, seconddimA, &alpha, A.fElem, A.Rows(), B.fElem, B.Rows(), &beta, ptr, this->Rows());
        return;
    }
    else if constexpr(std::is_same_v<TVar,std::complex<float>>){
        cblas_cgemm(CblasColMajor, transposeA, transposeB, nrows, ncols, seconddimA, &alpha, A.fElem, A.Rows(), B.fElem, B.Rows(), &beta, ptr, this->Rows());
        return;
    }

#endif //USING LAPACK

TPZFMatrix<TVar> tempA = A;
tempA *= alpha;
TPZFMatrix<TVar> tempB = B;
if (transpA)
    tempA.Transpose();
if (transpB)
    tempB.Transpose();

for (int64_t row = 0; row < nrows; row++)
{
	for (int64_t col = 0; col < ncols; col++)
	{
		for (int64_t k = 0; k < seconddimA; k++)
		{
			(*this)(i+row,j+col) += tempA(row,k) * tempB(k,col);
		}
	}
}

}

/********************************/
/*** Operator+=( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> & TPZFMatrix<TVar>::operator+=(const TPZMatrix<TVar> &aBase ) {
    auto A = dynamic_cast<const TPZFMatrix<TVar>*>(&aBase);
    if(!A){
        Error(__PRETTY_FUNCTION__," incompatible argument");
    }
    if ( (A->Rows() != this->Rows())  ||  (A->Cols() != this->Cols()) )
        Error( "Operator+= <matrices with different dimensions>" );
    
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * pm = fElem, *pmlast=pm+size;
    TVar * pa = A->fElem;
    while(pm < pmlast) (*pm++) += (*pa++);
    return( *this );
}

/*******************************/
/*** Operator-=( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator-=(const TPZMatrix<TVar> &aBase ) {
    auto A = dynamic_cast<const TPZFMatrix<TVar>*>(&aBase);
    if(!A){
        Error(__PRETTY_FUNCTION__," incompatible argument");
    }
    if ( (A->Rows() != this->Rows())  ||  (A->Cols() != this->Cols()) )
        Error( "Operator-= <matrixs with different dimensions>" );
    
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * pm = fElem;
    TVar * pa = A->fElem;
    
    for ( int64_t i = 0; i < size; i++ ) *pm++ -= *pa++;
    
    return( *this );
}

#ifdef USING_LAPACK
template <>
void TPZFMatrix<double>::ZAXPY(const double alpha,const TPZFMatrix<double> &p) {
    // 	//Como definir o tamanho dos vetores
    int size  = (fRow*fCol) ;
    cblas_daxpy(size, alpha, &p.fElem[0], 1, &fElem[0], 1);
}
#endif

template <class TVar>
void TPZFMatrix<TVar>::ZAXPY(const TVar alpha,const TPZFMatrix<TVar> &p) {
    
    TVar * pt = fElem;
    TVar * pp = p.fElem;
    TVar * ptlast = fElem + this->fRow*this->fCol;
    while(pt < ptlast) *pt++ += alpha * *pp++;
}

#ifdef USING_LAPACK
template<>
void TPZFMatrix<double>::TimesBetaPlusZ(const double beta,const TPZFMatrix<double> &z)
{
    int size = fRow*fCol;
    cblas_dscal(size,beta,fElem,1);
    cblas_daxpy(size,1.,z.fElem,1,fElem,1);
}
#endif

template<class TVar>
void TPZFMatrix<TVar>::TimesBetaPlusZ(const TVar beta,const TPZFMatrix<TVar> &z) {
    // #ifdef USING_ATLAS
    // 	int size = fRow*fCol;
    // 	cblas_dscal(size,beta,fElem,1);
    // 	cblas_daxpy(size,1.,z.fElem,1,fElem,1);
    // #else
    
    TVar * pt = fElem,  *ptlast = fElem + this->fRow*this->fCol;
    TVar * pz = z.fElem;
    //	while(pt < ptlast) *pt++ *= (beta) + *pz++;
    while(pt < ptlast) {
        *pt *= (beta);
        *pt++ += *pz++;
    }
    // #endif
}

/******** Operacoes com MATRIZES GENERICAS ********/

/******************/
/*** Operator = ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator=(const TPZMatrix<TVar> &A ) {
    int64_t arows  = A.Rows();
    int64_t acols  = A.Cols();
    int64_t size = arows * acols;
    if(fElem != fGiven) {
        delete []fElem;
        fElem = 0;
    }
    this->fRow  =  arows;
    this->fCol  = acols;
    if(fSize < size) {
        fElem = new TVar[ arows * acols ] ;
    } else {
        fElem = fGiven;
    }
    TVar * dst = fElem;
    for ( int64_t c = 0; c < this->fCol; c++ )
        for ( int64_t r = 0; r < this->fRow; r++ )
            *dst++ = A.Get( r, c );
    return( *this );
}


/******** Operacoes com valores NUMERICOS ********/

/******************/
/*** Operator = ***/
template <class TVar>
TPZFMatrix<TVar>& TPZFMatrix<TVar>::operator=(const TVar value ) {
    int64_t size = ((int64_t)this->fRow) * this->fCol;
    TVar * dst   = fElem;
    for ( int64_t i = 0; i < size; i++ )
        *dst++ = value;
    this->fDecomposed = ENoDecompose;
    return *this;
}



/***************************/
/*** Operator+=( value ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator+=(const TVar value ) {
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    
    TVar * dst = fElem, *dstlast = dst+size;
    while ( dst < dstlast ) *dst++ += value;
    return( *this );
}



/**************************/
/*** Operator+( value ) ***/
template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator+(const TVar value ) const {
    TPZFMatrix<TVar> res( *this );
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    
    TVar * dst = res.fElem,  *dstlast = dst+size;
    while ( dst < dstlast )
        *dst++ += value;
    
    return( res );
}

template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator-  (const TVar val ) const {
    return operator+( -val );
}


/**************************/
/*** Operator*( value ) ***/

template <class TVar>
TPZFMatrix<TVar>
TPZFMatrix<TVar>::operator*(const TVar value ) const
{
    auto res(*this);
    res *= value;
    return res;
}

/***************************/
/*** Operator*=( value ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator*=( const TVar value ) {
    int64_t size = ((int64_t)this->Rows()) * this->Cols();
    TVar * dst = fElem, *dstlast = dst+size;
    while ( dst < dstlast ) *dst++ *= value;
    return( *this );
}

/**************/
/*** Resize ***/
template <class TVar>
int TPZFMatrix<TVar>::Resize(const int64_t newRows,const int64_t newCols) {
    if ( newRows == this->Rows() && newCols == this->Cols() ) return( 1 );
    int64_t newsize = ((int64_t)newRows)*newCols;
    TVar * newElem{nullptr};
    if(fGiven && fElem != fGiven && newsize <= fSize)
    {
        newElem = fGiven;
    } else if(fGiven && fElem == fGiven && newsize <= fSize && newCols == 1)
    {
        this->fRow = newRows;
        this->fCol = newCols;
        return( 1 );
    } else if(fGiven && fElem == fGiven && newsize <= fSize && newRows == this->fRow)
    {
        this->fRow = newRows;
        this->fCol = newCols;
        return( 1 );
    } else
    {
        newElem = new TVar[ newRows * newCols ] ;
    }
    if ( newElem == NULL )
        Error( "Resize <memory allocation error>." );
    
    int64_t minRow  = ( this->fRow < newRows ? this->fRow : newRows );
    int64_t minCol  = ( this->fCol < newCols ? this->fCol : newCols );
    TVar * src;
    TVar * dst;
    int64_t r, c;
    
    for ( c = 0; c < minCol; c++ ) {
        // Copia as linhas da matriz antiga para a nova.
        // Copia os elementos de uma linha.
        dst = newElem + c*newRows;
        src = fElem + c*this->fRow;
        for ( r = 0; r < minRow; r++ ) *dst++ = *src++;
        
        // Se a nova linha for maior (mais colunas), preenche o
        //  resto da linha com ZEROS.
        for ( ; r < newRows; r++ ) *dst++ = 0.0;
    }
    
    // Preenche as linha que sobrarem (se sobrarem) com ZEROS.
    for ( ;c < newCols; c++ )
    {
        dst = newElem + c*newRows;
        for (r = 0 ; r < newRows; r++ ) *dst++ = 0.0;
    }
    
    if (fElem && fElem != fGiven )delete[]( fElem );
    fElem = newElem;
    this->fRow  = newRows;
    this->fCol  = newCols;
    return( 1 );
}

template <class TVar>
int TPZFMatrix<TVar>::Remodel(const int64_t newRows,const int64_t newCols) {
    if(newRows*newCols != this->fRow*this->fCol) return -1;
    this->fRow = newRows;
    this->fCol = newCols;
    return 1;
}

/********************/
/*** Transpose () ***/
template <class TVar>
void TPZFMatrix<TVar>::Transpose(TPZMatrix<TVar> *const T) const{
    T->Resize( this->Cols(), this->Rows() );
    //Transposta por filas
    TVar * p = fElem;
    for ( int64_t c = 0; c < this->Cols(); c++ ) {
        for ( int64_t r = 0; r < this->Rows(); r++ ) {
            T->PutVal( c, r, *p++ );
            //            cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
        }
    }
}

template <class TVar>
void TPZFMatrix<TVar>::Transpose() {
    TPZFMatrix<TVar> temp;
    Transpose(&temp);
    *this = temp;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_LU(TPZVec<int> &index) {
    
    
    if (this->fDecomposed != ENoDecompose && this->fDecomposed != ELUPivot) DebugStop();
    
    if (this->fDecomposed != ENoDecompose) {
        return ELUPivot;
    }
    
    if ( this->Rows() != this->Cols() ) {
        cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
        return 0;
    }
    
    
    lapack_int nRows = this->Rows();
    lapack_int zero = 0;
    float b;lapack_int info;
    
    InitializePivot();
    TPZManVector<lapack_int,2000> pivot(nRows);
    for(int i = 0; i < nRows; i++){pivot[i]=fPivot[i];}
    //    int sgesv_(__CLPK_integer *__n, __CLPK_integer *__nrhs, __CLPK_real *__a,
    //               __CLPK_integer *__lda, __CLPK_integer *__ipiv, __CLPK_real *__b,
    //               __CLPK_integer *__ldb,
    //               __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                __IPHONE_4_0);
    
    
    sgesv_(&nRows,&zero,fElem,&nRows,pivot.begin(),&b,&nRows,&info);
    if(info){
        std:cout << "\ninfo " << info << std::endl;
        DebugStop();
    }
    index.Resize(nRows);
    for(int i = 0; i < nRows; i++){
        index[i] = fPivot[i] = pivot[i];
    }
    this->fDecomposed = ELUPivot;
    return 1;
}
template <>
int TPZFMatrix<double>::Decompose_LU(TPZVec<int> &index) {
    
    
    if (this->fDecomposed != ENoDecompose && this->fDecomposed != ELUPivot) DebugStop();
    
    if (this->fDecomposed != ENoDecompose) {
        return ELUPivot;
    }
    
    if ( this->Rows() != this->Cols() ) {
        cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
        return 0;
    }
    
    
    lapack_int nRows = this->Rows();
    if (nRows == 0) return 0;

    lapack_int zero = 0;
    double b;lapack_int info;
    
    // If the matrix is 1x1, the lapack function dgesv_ does not modify
    // fPivot. And, if fPivot is not initialized it can lead to problems
    // when the function Substitution() is called. That is why we are
    // now initializing fPivot before calling dgesv_
//    fPivot.Resize(nRows);
    InitializePivot();
    TPZManVector<lapack_int,2000> pivot(nRows);
    for(int i = 0; i < nRows; i++){pivot[i]=fPivot[i];}
    //    int sgesv_(__CLPK_integer *__n, __CLPK_integer *__nrhs, __CLPK_real *__a,
    //               __CLPK_integer *__lda, __CLPK_integer *__ipiv, __CLPK_real *__b,
    //               __CLPK_integer *__ldb,
    //               __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                __IPHONE_4_0);
    
    
    dgetrf_(&nRows,&nRows,fElem,&nRows,&pivot[0],&info);
    if(info){
        std:cout << "\ninfo " << info << std::endl;
        DebugStop();
    }
//    dgesv_(&nRows,&zero,fElem,&nRows,&fPivot[0],&b,&nRows,&info);
    index.Resize(nRows);
    for(int i = 0; i < nRows; i++){
        index[i] = fPivot[i] = pivot[i];
    }
    this->fDecomposed = ELUPivot;
    return 1;
}

template <>
int TPZFMatrix<std::complex<double>>::Decompose_LU(TPZVec<int> &index) {
    
    
    if (this->fDecomposed != ENoDecompose && this->fDecomposed != ELUPivot) DebugStop();
    
    if (this->fDecomposed != ENoDecompose) {
        return ELUPivot;
    }
    
    if ( this->Rows() != this->Cols() ) {
        cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
        return 0;
    }
    
    
    lapack_int nRows = this->Rows();
    if (nRows == 0) return 0;

    lapack_int zero = 0;
    double b;lapack_int info;
    
    // If the matrix is 1x1, the lapack function dgesv_ does not modify
    // fPivot. And, if fPivot is not initialized it can lead to problems
    // when the function Substitution() is called. That is why we are
    // now initializing fPivot before calling dgesv_
//    fPivot.Resize(nRows);
    InitializePivot();
    TPZManVector<lapack_int,2000> pivot(nRows);
    for(int i = 0; i < nRows; i++){pivot[i]=fPivot[i];}
    //    int sgesv_(__CLPK_integer *__n, __CLPK_integer *__nrhs, __CLPK_real *__a,
    //               __CLPK_integer *__lda, __CLPK_integer *__ipiv, __CLPK_real *__b,
    //               __CLPK_integer *__ldb,
    //               __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                __IPHONE_4_0);
    

    static_assert(sizeof(vardoublecomplex) == sizeof(std::complex<double>),
                  "Incompatible types");
    zgetrf_(&nRows,&nRows,(vardoublecomplex*)fElem,&nRows,&pivot[0],&info);
    if(info){
        std:cout << "\ninfo " << info << std::endl;
        DebugStop();
    }
//    dgesv_(&nRows,&zero,fElem,&nRows,&fPivot[0],&b,&nRows,&info);
    index.Resize(nRows);
    for(int i = 0; i < nRows; i++){
        index[i] = fPivot[i] = pivot[i];
    }
    this->fDecomposed = ELUPivot;
    return 1;
}
#endif //USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU(TPZVec<int> &index) {
    
    if (this->fDecomposed == ELUPivot) return 1;
    if (this->fDecomposed) return 0;
    
    if ( this->Rows() != this->Cols() ) {
        cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
        return 0;
    }
    
    int64_t i,j,k;
    TVar sum = 0.;
    int64_t nRows = this->Rows();
    int64_t nCols = this->Cols();
    
    index.Resize(nRows);
    //inicializo o vetor de índices para o caso de pivotamento
    for (i=0;i<nRows;i++) index[i] = i;
    
    //LU
    for (j=0;j<nCols;j++){
        // cout << "line..." << j << endl;
        for (i=0;i<=j;i++){
            sum = 0.;
            for (k=0;k<i;k++){
                sum += this->Get(i,k) * this->Get(k,j);
            }
            TVar aux = this->Get(i,j);
            PutVal(i,j,(aux - sum));
            //cout << "0_A[" << i << "," << j << "]= " << Get(i,j) << endl;
        }
        //Print(cout);
        TVar piv = this->Get(j,j);
        //  cout << "Pivo 1 =" << piv << endl;
        int64_t row = j;
        for (i=j+1;i<nRows;i++){
            sum = 0.;
            for (k=0;k<(j);k++){
                sum += this->Get(i,k) * this->Get(k,j);
            }
            TVar aux = this->Get(i,j);
            PutVal(i,j,(aux - sum));
            //cout << "1_A[" << i << "," << j << "]= " << Get(i,j) << endl;
            
            if (fabs(this->Get(i,j)) > fabs(piv)){
                piv = this->Get(i,j);
                //  cout << "Pivo 2 =" << piv << endl;
                row = i;
            }
        }
        //    Print(cout);
        if (row > j){
            for (k=0;k<nCols;k++){
                TVar aux = this->Get(j,k);
                PutVal(j,k,this->Get(row,k));
                //cout << "2_A[" << j << "," << k << "]= " << Get(j,k) << endl;
                PutVal(row,k,aux);
                //cout << "3_A[" << row << "," << k << "]= " << Get(row,k) << endl;
            }
            k = index[j];
            index[j] = index[row];
            index[row] = k;
        }
        //    cout << "Pivo = " << piv << endl;
        for (i=j+1;i<nRows;i++){
            if (fabs(piv) < fabs((TVar)1e-12)) {
                cout << "Pivot < 1e-12. Probably matrix is singular. " << piv << endl;
                DebugStop();
            }
            TVar aux = this->Get(i,j) / piv;
            PutVal(i,j,aux);
            //cout << "4_A[" << i << "," << j << "]= " << Get(i,j) << endl;
        }
        //Print(cout);
    }
    this->fDecomposed = ELUPivot;
    return 1;
}


/*****************/
/*** DecomposeLU ***/
template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU(std::list<int64_t> &singular) {
    //return Decompose_LU();
#ifndef USING_LAPACK
    if (  this->fDecomposed && this->fDecomposed != ELU)  Error( "Decompose_LU <Matrix already Decomposed with other scheme>" );
#else
    if (  this->fDecomposed && this->fDecomposed != ELUPivot)  Error( "Decompose_LU <Matrix already Decomposed with other scheme>" );
#endif
    if (this->fDecomposed) return 1;
    
    const int  min = ( this->Cols() < (this->Rows()) ) ? this->Cols() : this->Rows();
    const int nrows = this->Rows();
    const int ncols = this->Cols();
    
    for ( int k = 0; k < min ; k++ ) {
        TVar pivot = GetVal(k, k);
        if (IsZero( pivot )){
            singular.push_back(k);
            PutVal(k,k,1.);
            pivot = 1.;
        }
        for ( int i = k+1; i < nrows; i++ ) {
            this->operator()(i,k) /= pivot;
        }
        for ( int j = k+1; j < ncols; j++ ){
            int i = k+1;
            TVar * elemPtr = &(this->operator()(i,j));
            TVar * ikPtr = &(this->operator()(i,k));
            const TVar kjVal = GetVal(k,j);
            for ( ; i < nrows; i++, elemPtr++, ikPtr++ ) {
                (*elemPtr) -= (*ikPtr)*kjVal;
                //          this->operator()(i,j) -= GetVal( i, k )*GetVal(k,j);
            }
        }
    }
    this->fDecomposed=ELU;
#ifdef USING_LAPACK
    fPivot.resize(nrows);
    for (int i=0; i<nrows; i++) {
        fPivot[i] = i+1;
    }
    this->fDecomposed = ELUPivot;
#endif //USING_LAPACK
    return 1;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_LU() {
    
    
    return this->Decompose_LU(fPivot);
}
template <>
int TPZFMatrix<double>::Decompose_LU() {
    
    
    return this->Decompose_LU(fPivot);
}
template <>
int TPZFMatrix<std::complex<double>>::Decompose_LU() {
    
    
    return this->Decompose_LU(fPivot);
}
#endif //USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU() {
    
    return this->Decompose_LU(fPivot);
}


template <class TVar>
int TPZFMatrix<TVar>::Substitution(const TVar *ptr, int64_t rows, TPZFMatrix<TVar> *B)
{
    int64_t rowb = B->Rows();
    int64_t colb = B->Cols();
    if ( rowb != rows ) Error( "static::SubstitutionLU <incompatible dimensions>" );
    int64_t i,j;
    for ( i = 0; i < rowb; i++ ) {
        for ( int64_t col = 0; col < colb; col++ )
            for (j = 0; j < i; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, j) * GETVAL(B, rowb, j, col));
    }
    
    for (int64_t col=0; col<colb; col++){
        for ( i = rowb-1; i >= 0; i-- ) {
            for (j = i+1; j < rowb ; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, j) * GETVAL(B, rowb, j, col));
            if ( IsZero( SELECTEL(ptr, rows, i, i)/*GetVal(i, i)*/ ) ) {
                if (fabs(SELECTEL(ptr, rows, i, i)/*GetVal(i, i)*/) > fabs((TVar)0.)) {
                    TVar tmp = GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, i);/*B->GetVal(i, col) - GetVal(i, i)*/
                    if (fabs(tmp) > fabs(((TVar)1e-12))) {
                        Error( "static::BackSub(SubstitutionLU) <Matrix is singular even after Power Plus..." );
                    }
                }else  Error( "static::BackSub(SubstitutionLU) <Matrix is singular" );
            }
            PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col)/SELECTEL(ptr, rows, i, i));
            //B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
        }
    }
    return( 1 );
}
//#endif

#include "fadType.h"

/****************/
/*** Substitution ***/
template <class TVar>
int TPZFMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const {
#ifdef USING_LAPACK    
    if (this->fDecomposed != ELUPivot) {
        Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
    }
    return Substitution(B,fPivot);
#else
    if(this->fDecomposed == ELUPivot){
        return Substitution(B,fPivot);
    }
    else if(this->fDecomposed != ELU) {
        Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
    }
#endif
    int64_t rowb = B->Rows();
    int64_t colb = B->Cols();
    int64_t row = this->Rows();
    if ( rowb != this->Rows() ) Error( "SubstitutionLU <incompatible dimensions>" );
    
    
    int64_t i,j;
    for ( i = 0; i < rowb; i++ ) {
        for ( int64_t col = 0; col < colb; col++ )
            for (j = 0; j < i; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - GETVAL(this, row, i, j) * GETVAL(B, rowb, j, col));
    }
    
    for (int64_t col=0; col<colb; col++){
        for ( i = rowb-1; i >= 0; i-- ) {
            for (j = i+1; j < rowb ; j++ )
                //B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
                PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - GETVAL(this, row, i, j) * GETVAL(B, rowb, j, col));
            if ( IsZero( GETVAL(this, row, i, i)/*GetVal(i, i)*/ ) ) {
                if (fabs(GETVAL(this, row, i, i)/*GetVal(i, i)*/) > 0.){
                    TVar diff = GETVAL(B, rowb, i, col) - GETVAL(this, row, i, i)/*B->GetVal(i, col) - GetVal(i, i)*/;
                    if (fabs(diff) > 1e-12){
                        Error( "BackSub(SubstitutionLU) <Matrix is singular even after Power Plus..." );
                    }
                }else  Error( "BackSub(SubstitutionLU) <Matrix is singular" );
            }
            PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col)/GETVAL(this, row, i, i));
            //B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
        }
    }
    return( 1 );
    
}


#ifdef USING_LAPACK
template<>
int TPZFMatrix<float>::Substitution( TPZFMatrix<float> *B, const TPZVec<int> &index ) const{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    TPZFMatrix<float> &b = *B;
    
    if (!this->fDecomposed){
        PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
        return 0;
    }
    
    if (this->fDecomposed != ELUPivot){
        PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
    }
    
    //    int sgetrs_(char *__trans, __CLPK_integer *__n, __CLPK_integer *__nrhs,
    //                __CLPK_real *__a, __CLPK_integer *__lda, __CLPK_integer *__ipiv,
    //                __CLPK_real *__b, __CLPK_integer *__ldb,
    //                __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                 __IPHONE_4_0);
    lapack_int nRows = this->Rows();
    char notrans = 'N';
    lapack_int BCols = B->Cols();
    lapack_int info = 0;

    TPZManVector<lapack_int,2000> pivot(nRows);
    for(int i = 0; i < nRows; i++){pivot[i]=index[i];}
    sgetrs_(&notrans,&nRows,&BCols,fElem,&nRows,&pivot[0],B->fElem,&nRows,&info);
#ifdef PZDEBUG
    if(info != 0)
    {
        DebugStop();
    }
#endif
    
    return 1;
}

template<>
int TPZFMatrix<double>::Substitution( TPZFMatrix<double> *B, const TPZVec<int> &index ) const{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    
    if (!this->fDecomposed){
        PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
        return 0;
    }
    
    if (this->fDecomposed != ELUPivot){
        PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
    }
    
    //    int sgetrs_(char *__trans, __CLPK_integer *__n, __CLPK_integer *__nrhs,
    //                __CLPK_real *__a, __CLPK_integer *__lda, __CLPK_integer *__ipiv,
    //                __CLPK_real *__b, __CLPK_integer *__ldb,
    //                __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                 __IPHONE_4_0);
    lapack_int nRows = this->Rows();
    char notrans = 'N';
    lapack_int BCols = B->Cols();
    lapack_int info = 0;

    TPZManVector<lapack_int,2000> pivot(nRows);
    for(int i = 0; i < nRows; i++){pivot[i]=index[i];}
    dgetrs_(&notrans,&nRows,&BCols,fElem,&nRows,&pivot[0],B->fElem,&nRows,&info);
    
#ifdef PZDEBUG
    if(info != 0)
    {
        DebugStop();
    }
#endif
    
    return 1;
}

template<>
int TPZFMatrix<std::complex<double>>::Substitution( TPZFMatrix<std::complex<double>> *B, const TPZVec<int> &index ) const{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    
    if (!this->fDecomposed){
        PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
        return 0;
    }
    
    if (this->fDecomposed != ELUPivot){
        PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
    }
    
    //    int sgetrs_(char *__trans, __CLPK_integer *__n, __CLPK_integer *__nrhs,
    //                __CLPK_real *__a, __CLPK_integer *__lda, __CLPK_integer *__ipiv,
    //                __CLPK_real *__b, __CLPK_integer *__ldb,
    //                __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
    //                                                                 __IPHONE_4_0);
    lapack_int nRows = this->Rows();
    char notrans = 'N';
    lapack_int BCols = B->Cols();
    lapack_int info = 0;

    TPZManVector<lapack_int,2000> pivot(nRows);
    for(int i = 0; i < nRows; i++){pivot[i]=index[i];}
    static_assert(sizeof(std::complex<double>)==sizeof(vardoublecomplex), "incompatible sizes");
    zgetrs_(&notrans,&nRows,&BCols,(vardoublecomplex*)fElem,&nRows,
            &pivot[0],(vardoublecomplex*)B->fElem,&nRows,&info);
    
#ifdef PZDEBUG
    if(info != 0)
    {
        DebugStop();
    }
#endif
    
    return 1;
}

#endif //USING_LAPACK


template<class TVar>
int TPZFMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B, const TPZVec<int> &index ) const{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    TPZFMatrix<TVar> &b = *B;
    
    if (!this->fDecomposed){
        PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
        return 0;
    }
    
    if (this->fDecomposed != ELUPivot){
        PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
    }
    
    int nRows = this->Rows();
    
    if (index.NElements() != nRows || b.Rows() != nRows)
    {
        cout << "TMatrix::Substituicao ERRO : vetores com dimensões incompativeis:\n"
        << "this->fIndex = " << index.NElements() << "  b = " << b.Rows() << endl;
        return 0;
    }
    
    int64_t ncols = B->Cols();
    for(int64_t ic = 0; ic<ncols; ic++)
    {
        int64_t i,j;
        TVar sum = 0;
        
        TPZVec<TVar> v(nRows);
        
        
        for (i=0;i<nRows;i++)
        {
            v[i] = b(index[i],ic);
        }
        
        //Ly=b
        for (i=0;i<nRows;i++)
        {
            sum = 0.;
            for (j=0;j<(i);j++) sum +=this->Get(i,j) * v[j];
            v[i] -= sum;
        }
        
        //Ux=y
        for (i=(nRows-1);i>-1;i--)
        {
            sum = 0.;
            for (j=(i+1);j<nRows;j++) sum += this->Get(i,j) * v[j];
            v[i] = (v[i] - sum) / this->Get(i,i);
        }
        
        for (i=0;i<nRows;i++) b(i,ic) = v[i];
    }
    return 1;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Substitution( TPZFMatrix<float> *B ) const {
    
    return this->Substitution(B,fPivot);
}
template <>
int TPZFMatrix<double>::Substitution( TPZFMatrix<double> *B ) const {
    
    return this->Substitution(B,fPivot);
}
template <>
int TPZFMatrix<std::complex<double>>::Substitution( TPZFMatrix<std::complex<double>> *B ) const {
    
    return this->Substitution(B,fPivot);
}
#endif //USING_LAPACK


#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_Cholesky() {
    if (  this->fDecomposed && this->fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows() != this->Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
    lapack_int dim=this->Dim();
    
    TPZFMatrix<float> B(*this);
    lapack_int nrhs = 0;
    float *A = fElem;
    char uplo = 'U';
    lapack_int info;
    //    sposv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    spotrf_(&uplo, &dim, A, &dim, &info);
    this->fDecomposed = ECholesky;
    
    if(info){
        std:cout << "\ninfo " << info << std::endl;
        DebugStop();
    }
    return 1;
}
template <>
int TPZFMatrix<double>::Decompose_Cholesky() {
    if (  this->fDecomposed && this->fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows() != this->Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
    lapack_int dim=this->Dim();
    
    double B;
    lapack_int nrhs = 0;
    double *A = fElem;
    char uplo = 'U';
    lapack_int info;
    dpotrf_(&uplo, &dim, A, &dim, &info);
    this->fDecomposed = ECholesky;
    
    if(info){
        std:cout << "\ninfo " << info << std::endl;
        DebugStop();
    }
    return 1;
}
#endif //USING_LAPACK


template <class TVar>
int TPZFMatrix<TVar>::Decompose_Cholesky(){
    if (  this->fDecomposed && this->fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows() != this->Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
    //return 0;
    
    int dim=this->Dim();
    
    for (int i=0 ; i<dim; i++) {
        

        TVar &diagII = this->g(i,i);
        for(int k=0; k<i; k++) { //diagonal elements
            TVar sum = 0;
            if constexpr (is_complex<TVar>::value){
                sum += this->operator()(i,k)*std::conj(this->operator()(i,k));
            }else{
                sum+= this->operator()(i,k)*this->operator()(i,k);
            }
            diagII -= sum;
        }
        
        
        if( IsZero(diagII) ){
            DebugStop();
        }
        if constexpr(!is_complex<TVar>::value) if (diagII < 0.) DebugStop();
        
        diagII = sqrt(diagII);

        for (int j=i+1;j<dim; j++) {//off-diagonal elements
            TVar sum = 0.;
            int k = 0;
            TVar * ikPtr = &(this->g(k,i));///(k,i) = (i,k) given that the matrix is symmetric, but the alignment will speed up execution
            TVar * kjPtr = &(this->g(k,j));
            for(; k<i; k++, kjPtr++, ikPtr++) {
                if constexpr(is_complex<TVar>::value){
                    sum += std::conj(*ikPtr)*(*kjPtr);
                }else{
                    sum += (*ikPtr)*(*kjPtr);
                }
            }
            
            const TVar val = (this->GetVal(i,j) - sum)/(diagII);
            this->PutVal(i,j,val);
            if constexpr (is_complex<TVar>::value) this->PutVal(j,i,std::conj(val));
            else this->PutVal(j,i,val);
        }
    }
    
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    this->fDecomposed = ECholesky;
    return ECholesky;
}


template <class TVar>
int TPZFMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular) {
    
    if (  this->fDecomposed && this->fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows() != this->Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
    //return 0;
    
    int dim=this->Dim();
    
    for (int i=0 ; i<dim; i++) {
        

        TVar &diagII = this->g(i,i);
        for(int k=0; k<i; k++) { //diagonal elements
            TVar sum = 0;
            if constexpr (is_complex<TVar>::value){
                sum += this->operator()(i,k)*std::conj(this->operator()(i,k));
            }else{
                sum+= this->operator()(i,k)*this->operator()(i,k);
            }
            diagII -= sum;
        }
        
        
        if( IsZero(diagII) ){
            singular.push_back(i);
            (diagII) = 1.;
        }
        diagII = sqrt(diagII);

        for (int j=i+1;j<dim; j++) {//off-diagonal elements
            TVar sum = 0.;
            int k = 0;
            TVar * ikPtr = &(this->g(k,i));///(k,i) = (i,k) given that the matrix is symmetric, but the alignment will speed up execution
            TVar * kjPtr = &(this->g(k,j));
            for(; k<i; k++, kjPtr++, ikPtr++) {
                if constexpr(is_complex<TVar>::value){
                    sum += std::conj(*ikPtr)*(*kjPtr);
                }else{
                    sum += (*ikPtr)*(*kjPtr);
                }
            }
            
            const TVar val = (this->GetVal(i,j) - sum)/(diagII);
            this->PutVal(i,j,val);
            if constexpr (is_complex<TVar>::value) this->PutVal(j,i,std::conj(val));
            else this->PutVal(j,i,val);
        }
    }
    
    //    std::cout << __PRETTY_FUNCTION__ << std::endl;
    this->fDecomposed = ECholesky;
    return ECholesky;
}

template <class TVar>
int TPZFMatrix<TVar>::Substitution(const TVar *ptr, int64_t rows, TPZFMatrix<TVar> *B, const TPZVec<int> &index )
{
    
    if(!B){
        PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
        return 0;
    }
    
    TPZFMatrix<TVar> &b = *B;
    
    
    
    if (index.NElements() != rows || b.Rows() != rows)
    {
        cout << "TMatrix::Substituicao ERRO : vetores com dimensões incompativeis:\n"
        << "this->fIndex = " << index.NElements() << "  b = " << b.Rows() << endl;
        return 0;
    }
    
    int64_t i,j;
    TVar sum = 0;
    
    TPZVec<TVar> v(rows);
    
    
    for (i=0;i<rows;i++)
    {
        v[i] = b(index[i]);
    }
    
    //Ly=b
    for (i=0;i<rows;i++)
    {
        sum = 0.;
        for (j=0;j<(i);j++) sum += SELECTEL(ptr,rows,i,j) * v[j];
        v[i] -= sum;
    }
    
    //Ux=y
    for (i=(rows-1);i>-1;i--)
    {
        sum = 0.;
        for (j=(i+1);j<rows;j++) sum += SELECTEL(ptr,rows,i,j) * v[j];
        v[i] = (v[i] - sum) / SELECTEL(ptr,rows,i,i);
    }
    
    for (i=0;i<rows;i++) b(i) = v[i];
    return 1;
}

#ifdef USING_LAPACK
template <>
int TPZFMatrix<float>::Decompose_LDLt() {
    
    if (  this->fDecomposed && this->fDecomposed != ELDLt) {
        Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
    } else if(this->fDecomposed ) {
        return ELDLt;
    }
    if ( this->Rows()!=this->Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
    char uplo = 'U';
    lapack_int dim = Rows();
    lapack_int nrhs = 0;
    fPivot.Resize(dim,0);
    float B  = 0.;
    lapack_int worksize = 3*dim;
    fWork.Resize(worksize);
    lapack_int info;
    
    if (dim == 0) {
        this->fDecomposed  = ELDLt;
        this->fDefPositive = 0;
        return( 1 );
    }
    
    //    ssysv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_real *__work#>, <#__CLPK_integer *__lwork#>, <#__CLPK_integer *__info#>)

    InitializePivot();
    TPZManVector<lapack_int,2000> pivot(dim);
    for(int i = 0 ; i < dim ; i++){pivot[i] = fPivot[i];}
    ssysv_(&uplo, &dim, &nrhs, fElem, &dim, &pivot[0], &B, &dim, &fWork[0], &worksize, &info);
    
    if(info){
        std:cout << "\ninfo " << info << std::endl;
        DebugStop();
    }
    
    for(int i = 0; i < dim; i++){fPivot[i] = pivot[i];}
    fDecomposed = ELDLt;
    
    return 1;
}

template <>
int TPZFMatrix<double>::Decompose_LDLt() {
    
    if (  this->fDecomposed && this->fDecomposed != ELDLt) {
        Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
    } else if(this->fDecomposed ) {
        return ELDLt;
    }
    if ( this->Rows()!=this->Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
    char uplo = 'L';
    lapack_int dim = Rows();
    lapack_int nrhs = 0;
    fPivot.Resize(dim,0);
    double B  = 0.;
    lapack_int worksize = 3*dim;
    fWork.Resize(worksize);
    lapack_int info;
    
    if (dim == 0) {
        this->fDecomposed  = ELDLt;
        this->fDefPositive = 0;
        return( 1 );
    }
    
    //    ssysv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_real *__work#>, <#__CLPK_integer *__lwork#>, <#__CLPK_integer *__info#>)
    TPZManVector<lapack_int,2000> pivot(dim);
    for(int i = 0 ; i < dim ; i++){pivot[i] = fPivot[i];}
    dsysv_(&uplo, &dim, &nrhs, fElem, &dim, &pivot[0], &B, &dim, &fWork[0], &worksize, &info);
    
    if(info){
        std:cout << "\ninfo " << info << std::endl;
        DebugStop();
    }
    
    for(int i = 0; i < dim; i++){fPivot[i] = pivot[i];}
    fDecomposed = ELDLt;
    return 1;
}

#endif //USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LDLt() {
    if (  this->fDecomposed && this->fDecomposed != ELDLt) {
        Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
    } else if(this->fDecomposed ) {
        return ELDLt;
    }
    if ( this->Rows()!=this->Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
    
    const int dim=this->Rows();
    
    for (auto j = 0; j < dim; j++ ) {
        TVar sum = 0;
        for (auto k=0; k<j; k++) {
            if constexpr(is_complex<TVar>::value){
                sum +=GetVal(k,k)*std::conj(GetVal(k,j))*GetVal(k,j);
            }else{
                sum +=GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
            }
        }
        PutVal( j,j,GetVal(j,j) - sum );
        TVar tmp = GetVal(j,j);
        if ( IsZero(tmp) ) Error( "Decompose_LDLt <Zero on diagonal>" );
        for(auto l=j+1; l<dim;l++) {
            TVar sum = 0;
            for (auto k=0; k<j; k++) {
                if constexpr(is_complex<TVar>::value){
                    sum += GetVal(k,k)*std::conj(GetVal(j,k))*GetVal(l,k);
                }else{
                    sum += GetVal(k,k)*GetVal(j,k)*GetVal(l,k);
                }
            }
            TVar val = (GetVal(l,j) - sum)/tmp;
            PutVal(l,j,val);
            if constexpr(is_complex<TVar>::value){val=std::conj(val);}
            PutVal(j,l,val);
        }
    }
    this->fDecomposed  = ELDLt;
    this->fDefPositive = 0;
    return( 1 );
}


#ifdef USING_LAPACK
/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_Forward( TPZFMatrix<float>* b ) const
{
    if (fDecomposed == ECholesky) {
        //        CALL strsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
        //                   179      $               one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        //        b->Print("before =",std::cout,EMathematicaInput);
        cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        TPZMatrix<float>::Subst_Forward(b);
        //        b->Print("after =",std::cout,EMathematicaInput);
        return 1;
    }
    else
    {
        DebugStop();
        return 0;
    }
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_Forward( TPZFMatrix<double>* b ) const
{
    if (fDecomposed == ECholesky) {
        //        CALL strsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
        //                   179      $               one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        //        b->Print("before =",std::cout,EMathematicaInput);
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        TPZMatrix<float>::Subst_Forward(b);
        //        b->Print("after =",std::cout,EMathematicaInput);
        return 1;
    }
    else
    {
        DebugStop();
        return 0;
    }
}
#endif 
/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>* B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed != ECholesky)
        return( 0 );
    for ( int64_t r = 0; r < this->Dim(); r++ ) {
        TVar pivot = GetVal( r, r );
        for ( int64_t c = 0; c < B->Cols();  c++ ) {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = 0, ..., r-1.
            //
            TVar sum = 0.0;
            for ( int64_t i = 0; i < r; i++ ) sum += GetVal(r, i) * B->GetVal(i, c);
            
            // Faz B[r,c] = (B[r,c] - sum) / A[r,r].
            //
            B->PutVal( r, c, (B->GetVal(r, c) - sum) / pivot );
        }
    }
    return( 1 );
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar>* B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed != ECholesky) return( 0 );
    for ( int64_t r = this->Dim()-1;  r >= 0;  r-- ) {
        TVar pivot = GetVal( r, r );
        for ( int64_t c = 0; c < B->Cols(); c++ ) {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = N, ..., r+1.
            //
            TVar sum = 0.0;
            for ( int64_t i = this->Dim()-1; i > r; i-- ) sum += GetVal(r, i) * B->GetVal(i, c);
            // Faz B[r,c] = (B[r,c] - sum) / A[r,r].
            //
            B->PutVal( r, c, (B->GetVal(r, c) - sum) / pivot );
        }
    }
    return( 1 );
}

#ifdef USING_LAPACK
template<>
int TPZFMatrix<float>::Subst_Backward( TPZFMatrix<float>* b ) const
{
    if (fDecomposed == ECholesky) {
        
        //        spotrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
        //        CALL strsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
        //                   184      $               nrhs, one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        int info = 0;
        //        b->Print("before =",std::cout,EMathematicaInput);
        //        TPZMatrix<float>::Subst_Backward(b);
        cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        spotrs_(upper, &dim, &ncol, fElem, &dim, b->fElem, &ncol, &info);
        //        b->Print("after =",std::cout,EMathematicaInput);
        
        return 1;
        
        if (info != 0) {
            DebugStop();
        }
    }
    else
    {
        DebugStop();
        return 0;
    }
}

template<>
int TPZFMatrix<double>::Subst_Backward( TPZFMatrix<double>* b ) const
{
    if (fDecomposed == ECholesky) {
        
        //        spotrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
        //        CALL strsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
        //                   184      $               nrhs, one, a, lda, b, ldb )
        char left[]= "Left", upper[] = "Upper", transpose[] = "Transpose", non_unit[] = "Non-unit";
        int ncol = b->Cols();
        int dim = Rows();
        int info = 0;
        //        b->Print("before =",std::cout,EMathematicaInput);
        //        TPZMatrix<float>::Subst_Backward(b);
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, dim, ncol, 1., fElem, dim, b->fElem, dim);
        //        spotrs_(upper, &dim, &ncol, fElem, &dim, b->fElem, &ncol, &info);
        //        b->Print("after =",std::cout,EMathematicaInput);
        
        return 1;
        
        if (info != 0) {
            DebugStop();
        }
    }
    else
    {
        DebugStop();
        return 0;
    }
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_LForward( TPZFMatrix<float>* b ) const
{
    //    ssytrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    
    if (fDecomposed != ELDLt) {
        DebugStop();
    }
    
    char uplo = 'U';
    lapack_int dim = Rows();
    lapack_int nrhs = b->Cols();
    float B  = 0.;
    lapack_int info;
    
    //    ssytrs_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__a#>, <#__CLPK_integer *__lda#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    TPZManVector<lapack_int,2000> pivot(dim);
    for(int i = 0 ; i < dim ; i++){pivot[i] = fPivot[i];}
    ssytrs_(&uplo, &dim, &nrhs, fElem, &dim, &pivot[0], b->fElem, &dim, &info);
    for(int i = 0; i < dim; i++){fPivot[i] = pivot[i];}
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_LForward( TPZFMatrix<double>* b ) const
{
   
    if (fDecomposed != ELDLt) {
        DebugStop();
    }
    
    char uplo = 'L';
    lapack_int dim = Rows();
    lapack_int nrhs = b->Cols();
    double B  = 0.;
    lapack_int info;
    if (dim == 0 || nrhs == 0) {
        return 0;
    }
    TPZManVector<lapack_int,2000> pivot(dim);
    for(int i = 0 ; i < dim ; i++){pivot[i] = fPivot[i];}
    dsytrs_(&uplo, &dim, &nrhs, fElem, &dim, &pivot[0], b->fElem, &dim, &info);
    for(int i = 0; i < dim; i++){fPivot[i] = pivot[i];}
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_LBackward( TPZFMatrix<float>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_LBackward( TPZFMatrix<double>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<float>::Subst_Diag( TPZFMatrix<float>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<>
int TPZFMatrix<double>::Subst_Diag( TPZFMatrix<double>* b ) const
{
    return 1;
}

#endif
/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar>* B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed != ELDLt) {
        Error("TPZMatrix::Subst_LForward incompatible dimensions\n");
    }
    for ( int64_t r = 0; r < this->Dim(); r++ ) {
        for ( int64_t c = 0; c < B->Cols();  c++ )    {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = 0, ..., r-1.
            //
            TVar sum = 0.0;
            for ( int64_t i = 0; i < r; i++ ) sum += GetVal(r, i) * B->GetVal(i, c);
            
            // Faz B[r,c] = (B[r,c] - sum) / A[r,r].
            //
            B->PutVal( r, c, B->GetVal(r, c) - sum );
        }
    }
    return( 1 );
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar>* B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed != ELDLt){
        Error("TPZMatrix::Subst_LBackward incompatible dimensions \n");
    }
    
    for ( int64_t r = this->Dim()-1;  r >= 0;  r-- ) {
        for ( int64_t c = 0; c < B->Cols(); c++ ) {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = N, ..., r+1.
            //
            TVar sum = 0.0;
            for ( int64_t i = this->Dim()-1; i > r; i-- ) sum += GetVal(r, i) * B->GetVal(i, c);
            
            // Faz B[r,c] = B[r,c] - sum.
            //
            B->PutVal( r, c, B->GetVal(r, c) - sum );
        }
    }
    return( 1 );
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZFMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar>* B ) const
{
    if ( (B->Rows() != this->Dim())) {
        Error("TPZMatrix::Subst_Diag incompatible dimensions\n");
    }
    for ( int64_t r = 0; r < this->Dim(); r++ ) {
        TVar pivot = GetVal( r, r );
        for ( int64_t c = 0; c < B->Cols(); c++ ) {
            B->PutVal( r, c, B->GetVal( r, c ) / pivot );
        }
    }
    return( 1 );
}

/** @brief Implement dot product for matrices */
template<class TVar>
TVar Dot(const TPZFMatrix<TVar> &A, const TPZFMatrix<TVar> &B) {
    int64_t size = (A.Rows())*A.Cols();
    TVar result = 0.;
    if(!size) return result;
    const TVar *fpA = &A.g(0,0), *fpB = &B.g(0,0);
    const TVar *fpLast = fpA+size;
    while(fpA < fpLast)
    {
        if constexpr (is_complex<TVar>::value){
            result += *fpA++ * std::conj(*fpB++);
        }
        else{
            result += (*fpA++ * *fpB++);
        }
    }
    return result;
    // #endif
}

template
std::complex<float> Dot(const TPZFMatrix< std::complex<float> > &A, const TPZFMatrix< std::complex<float> > &B);

template
std::complex<double> Dot(const TPZFMatrix< std::complex<double> > &A, const TPZFMatrix< std::complex<double> > &B);

template
std::complex<long double> Dot(const TPZFMatrix< std::complex<long double> > &A, const TPZFMatrix< std::complex<long double> > &B);

template
long double Dot(const TPZFMatrix<long double> &A, const TPZFMatrix<long double> &B);

template
double Dot(const TPZFMatrix<double> &A, const TPZFMatrix<double> &B);

template
float Dot(const TPZFMatrix<float> &A, const TPZFMatrix<float> &B);

template
int64_t Dot(const TPZFMatrix<int64_t> &A, const TPZFMatrix<int64_t> &B);

template
int Dot(const TPZFMatrix<int> &A, const TPZFMatrix<int> &B);

template
Fad<float> Dot(const TPZFMatrix<Fad<float> > &A, const TPZFMatrix<Fad<float> > &B);

template
Fad<double> Dot(const TPZFMatrix<Fad<double> > &A, const TPZFMatrix<Fad<double> > &B);

template
Fad<long double> Dot(const TPZFMatrix<Fad<long double> > &A, const TPZFMatrix<Fad<long double> > &B);

template
TPZFlopCounter Dot(const TPZFMatrix<TPZFlopCounter> &A, const TPZFMatrix<TPZFlopCounter> &B);

/** @brief Increments value over all entries of the matrix A. */
template <class TVar>
TPZFMatrix<TVar> operator+(const TVar value, const TPZFMatrix<TVar> &A ) {
    return( A + value );
}

/** @brief Decrements value over all entries of the matrix A. */
template <class TVar>
TPZFMatrix<TVar> operator-(const TVar value, const TPZFMatrix<TVar> &A ) {
    return( A - value );
}

/************************** Private **************************/

/*************/
/*** Error ***/
template<class TVar>
int TPZFMatrix<TVar>::Error(const char *msg1,const char *msg2 ) {
    ostringstream out;
    out << "TPZFMatrix::" << msg1;
    if(msg2) out << msg2;
    out << ".\n";
    LOGPZ_ERROR (logger, out.str().c_str());
    DebugStop();
    return 0;
}


/*************/
/*** Clear ***/
template <class TVar>
int TPZFMatrix<TVar>::Clear() {
    if(fElem && fElem != fGiven) delete[]( fElem );
    fElem = NULL;
    this->fRow  = this->fCol = 0;
    return( 1 );
}


template <class TVar>
void TPZFMatrix<TVar>::Read( TPZStream &buf, void *context ){ //ok
    if constexpr (std::is_same<TVar, TFad<6, REAL>>::value ||
                  std::is_same<TVar, Fad<float>>::value ||
                  std::is_same<TVar, Fad<double>>::value ||
                  std::is_same<TVar, Fad<long double>>::value||
                  std::is_same<TVar, TPZFlopCounter>::value) {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" not implemented for this type\n";
        PZError<<"Aborting...";
        DebugStop();
    }
    TPZMatrix<TVar>::Read(buf,context);
    int64_t row = this->fRow;
    int64_t col = this->fCol;
    //this is odd, but necessary.
    this->fRow = this->fCol = 0;
    Resize(row,col);
    buf.Read(fElem,this->fRow*this->fCol);
}

//template <class TVar>
//void TPZFMatrix<TVar>::Write( TPZStream &buf, int withclassid ) {
//    const TPZFMatrix<TVar> *cp = this;
//    cp->Write(buf,withclassid);
//    //    const Write(buf, withclassid);
//    //	TPZMatrix<TVar>::Write(buf,withclassid);
//    //	buf.Write(fElem,this->fRow*this->fCol);
//}

template <class TVar>
void TPZFMatrix<TVar>::Write( TPZStream &buf, int withclassid ) const { //ok
    if constexpr (std::is_same<TVar, TFad<6, REAL>>::value ||
                  std::is_same<TVar, Fad<float>>::value ||
                  std::is_same<TVar, Fad<double>>::value ||
                  std::is_same<TVar, Fad<long double>>::value||
                  std::is_same<TVar, TPZFlopCounter>::value) {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" not implemented for this type\n";
        PZError<<"Aborting...";
        DebugStop();
    }
    TPZMatrix<TVar>::Write(buf,withclassid);
    buf.Write(fElem,this->fRow*this->fCol);
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template<class TVar>
bool TPZFMatrix<TVar>::Compare(TPZSavable *copy, bool override)
{
    TPZFMatrix<TVar> *fmat = dynamic_cast<TPZFMatrix<TVar> *> (copy);
    if(!fmat) return false;
    
    bool matresult = TPZMatrix<TVar>::Compare(copy,false);
    int64_t nel = this->fRow*this->fCol;
    TVar diff=0.;
    int64_t numdif = 0;
    int64_t iel;
    for(iel=0; iel<nel; iel++)
    {
        if(fElem[iel] != fmat->fElem[iel])
        {
            matresult = false;
            numdif++;
        }
        TVar exp = fElem[iel]-fmat->fElem[iel];
        diff += fabs(exp);
    }
    if(!matresult)
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " did not compare ";
        sout << " number different terms " << numdif << " number terms " << this->fRow*this->fCol;
        sout << " difference in norm L1 " << diff;
        LOGPZ_ERROR(loggerCheck,sout.str())
    }
    if(!matresult && override)
    {
        this->operator=(*fmat);
    }
    return matresult;
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template<class TVar>
bool TPZFMatrix<TVar>::Compare(TPZSavable *copy, bool override) const
{
    TPZFMatrix<TVar> *fmat = dynamic_cast<TPZFMatrix<TVar> *> (copy);
    if(!fmat) return false;
    
    bool matresult = TPZMatrix<TVar>::Compare(copy,false);
    int64_t nel = this->fRow*this->fCol;
    TVar diff=0.;
    int64_t numdif = 0;
    int64_t iel;
    for(iel=0; iel<nel; iel++)
    {
        if(fElem[iel] != fmat->fElem[iel])
        {
            matresult = false;
            numdif++;
        }
        TVar exp = fElem[iel]-fmat->fElem[iel];
        diff += fabs(exp);
    }
    if(!matresult)
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " did not compare ";
        sout << " number different terms " << numdif << " number terms " << this->fRow*this->fCol;
        sout << " difference in norm L1 " << diff;
        LOGPZ_ERROR(loggerCheck,sout.str())
    }
    if(!matresult && override)
    {
        DebugStop();
    }
    return matresult;
}

template <>
void TPZFMatrix<TFad<6,REAL> >::PrintStatic(const TFad<6,REAL> *ptr, int64_t rows, int64_t cols, const char *name, std::ostream& out,const MatrixOutputFormat form)
{
    DebugStop();
}

template <class TVar>
void TPZFMatrix<TVar>::PrintStatic(const TVar *ptr, int64_t rows, int64_t cols, const char *name, std::ostream& out,const MatrixOutputFormat form){
    
    if(form == EFormatted) {
        out << "Writing matrix '";
        if(name) out << name;
        out << "' (" << rows << " x " << cols << "):\n";
        
        for ( int64_t row = 0; row < rows; row++) {
            out << "\t";
            for ( int64_t col = 0; col < cols; col++ ) {
                out << SELECTEL(ptr,rows, row, col) << "  ";
            }
            out << "\n";
        }
        out << "\n";
    } else if (form == EInputFormat) {
        out << rows << " " << cols << endl;
        for ( int64_t row = 0; row < rows; row++) {
            for ( int64_t col = 0; col < cols; col++ ) {
                TVar val = SELECTEL(ptr,rows,row, col);
                if(val != (TVar)0.) out << row << ' ' << col << ' ' << val << std::endl;
            }
        }
        out << "-1 -1 0.\n";
    } else if( form == EMathematicaInput)
    {
        char number[32];
        out << name << "\n{ ";
        for ( int64_t row = 0; row < rows; row++) {
            out << "\n{ ";
            for ( int64_t col = 0; col < cols; col++ ) {
                TVar val = SELECTEL(ptr,rows,row, col);
                snprintf(number, 32, "%16.16lf", (double)fabs(val));
                out << number;
                if(col < cols-1)
                    out << ", ";
                if((col+1) % 6 == 0)out << std::endl;
            }
            out << " }";
            if(row < rows-1)
                out << ",";
        }
        
        out << " }\n";
        
    }
    
}

template<class TVar>
int TPZFMatrix<TVar>::ClassId() const{
    return Hash("TPZFMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}

template <class TVar>
int TPZFMatrix<TVar>::SetSize(const int64_t newRows,const int64_t newCols) {
    int64_t newsize = ((int64_t)newRows)*newCols;
    int64_t oldsize = this->fRow*this->fCol;
    if(newsize == oldsize) return 1;
    if(fElem && fElem != fGiven)
    {
        delete []fElem;
        fElem = 0;
    }
    if(fGiven && newsize <= fSize)
    {
        fElem = fGiven;
    } else
    {
        fElem = new TVar[ newRows * newCols ] ;
    }
    if (newsize && fElem == NULL )
        Error( "Resize <memory allocation error>." );
    
    return( 1 );
}

#ifdef USING_LAPACK

template <class TVar>
int TPZFMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &eigenvalues)
{
    if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                  && is_arithmetic_pz<TVar>::value){
        TPZLapackEigenSolver<TVar> solver;
        return solver.SolveEigenProblem(*this,eigenvalues);
    }else{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Incompatible types.\nAborting...\n";
        DebugStop();
        return -1;
    }
}

template <class TVar>
int TPZFMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &eigenvalues,
                                        TPZFMatrix < CTVar > &eigenvectors)
{
    if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                  && is_arithmetic_pz<TVar>::value){
        TPZLapackEigenSolver<TVar> solver;
        return solver.SolveEigenProblem(*this,eigenvalues,eigenvectors);
    }else{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Incompatible types.\nAborting...\n";
        DebugStop();
        return -1;
    }
}


template< class TVar>
int
TPZFMatrix<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &B ,
                                               TPZVec < CTVar > &w,
                                               TPZFMatrix < CTVar > &eigenvectors)
{
    if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                  && is_arithmetic_pz<TVar>::value){
        TPZLapackEigenSolver<TVar> solver;
        return solver.SolveGeneralisedEigenProblem(*this,B,w,eigenvectors);
    }else{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Incompatible types.\nAborting...\n";
        DebugStop();
        return -1;
    }
}
template< class TVar>
int
TPZFMatrix<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &B ,
                                               TPZVec <CTVar> &w)
{
    if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                  && is_arithmetic_pz<TVar>::value){
        TPZLapackEigenSolver<TVar> solver;
        return solver.SolveGeneralisedEigenProblem(*this,B,w);
    }else{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Incompatible types.\nAborting...\n";
        DebugStop();
        return -1;
    }
}


template<typename TVar>
int TPZFMatrix<TVar>::SingularValueDecomposition(TPZFMatrix<TVar>& U, TPZFMatrix<TVar>& S, TPZFMatrix<TVar>& VT,char jobU, char jobVT){
    PZError<<__PRETTY_FUNCTION__
           <<"\n\tNot implemented for your type yet."
           <<"\n\tYou're welcome to pull-request your implementation to github.com/labmec/neopz \n";
    DebugStop();
    return -1;
}

template<>
int TPZFMatrix<double>::SingularValueDecomposition(TPZFMatrix<double>& U, TPZFMatrix<double>& S, TPZFMatrix<double>& VT,char jobU, char jobVT){
    // Setup matrix sizes. This section of the code should happen for all specializations
    lapack_int m = this->Rows();
    lapack_int n = this->Cols();
    lapack_int min = std::min(m,n);
    switch(jobU){
        case 'A': U.Resize(m,m);    break;
        case 'S': U.Resize(m,min);  break;
        case 'N': U.Resize(1,1);    break;
        default: PZError << "\nInvalid input jobU = \'" << jobU << "\' unrecognized;\n"; DebugStop();
    }
    switch(jobVT){
        case 'A': VT.Resize(n,n);   break;
        case 'S': VT.Resize(min,n); break;
        case 'N': VT.Resize(1,1);   break;
        default: PZError << "\nInvalid input jobVT = \'" << jobVT << "\' unrecognized;\n"; DebugStop();
    }
    S.Resize(min,1);

    // Get matrices pointers
    double* A_ptr = &(this->operator()(0,0));
    double* U_ptr = &U(0,0);
    double* S_ptr = &S(0,0);
    double* VT_ptr = &VT(0,0);
    // Setup auxiliar variables
    lapack_int nrows = fRow, ncols = fCol;
    lapack_int info = 0;
    lapack_int lda = nrows;
    lapack_int ldu = 1, ldvt = 1;
    switch(jobU){
        case 'A':
        case 'S': ldu = nrows; break;
        default: ldu = 1; break;
    }
    switch(jobVT){
        case 'A': ldvt = ncols; break;
        case 'S': ldvt = min; break;
        default:  ldvt = 1; break;
    }
    double work_opt;
    lapack_int lwork = -1; //<-- Pass -1 to tell Lapack to compute it for you
    // first do a pseudo-run to compute optimal work size
    dgesvd_(&jobU,&jobVT,&nrows,&ncols,A_ptr,&lda,S_ptr,U_ptr,&ldu,VT_ptr,&ldvt,&work_opt,&lwork,&info);
    lwork = (lapack_int)work_opt;
    TPZVec<double> work(lwork,0.);
    // then do actual computation of SVD
    dgesvd_(&jobU,&jobVT,&nrows,&ncols,A_ptr,&lda,S_ptr,U_ptr,&ldu,VT_ptr,&ldvt,&work[0],&lwork,&info);
    if(info  ==  0){/*all ok*/}
    else if(info>0){
        PZError<<"\nSingularValueDecomposition failed to converge\n";
        PZError<<"\tCheck matrix before running lapack::dgesvd since it overwrites the matrix\n";
        DebugStop();
        // return info; // in case you don't want a DebugStop()
    }
    else if(info<0){
        PZError<<"\nThe "<<-info<<"th parameter passed to lapack::dgesvd is a bad parameter\n";
        DebugStop();
    }

    const int success = 1;
    return success;
}


#endif // USING_LAPACK

/******************/
/*** Destructor ***/
template<class TVar>
TPZFMatrix<TVar>::~TPZFMatrix() {
    if(fElem && fElem != fGiven) delete[]( fElem );
    fElem = 0;
    fSize = 0;
}

/** @brief Returns the norm of the matrix A */
template<>
TFad<6,REAL> Norm(const TPZFMatrix<TFad<6,REAL> > &A)
{
    DebugStop();
    TFad<6,REAL> res;
    return res;
}
template<>
Fad<REAL> Norm(const TPZFMatrix<Fad<REAL> > &A)
{
    DebugStop();
    Fad<REAL> res;
    return res;
}


template<class TVar>
void TPZFMatrix<TVar>::Symetrize() {
  
  if ( this->Rows() != this->Cols() ) {
    Error( "Simetrize only work for square matrices" );
  }
  
  int64_t row,col;
  int64_t fDim1 = this->Rows();
  const bool must_conj =
    is_complex<TVar>::value && this->GetSymmetry() == SymProp::Herm;
  for(row=0; row<fDim1; row++) {
    if(must_conj){
      if constexpr (is_complex<TVar>::value){
        for(col=row+1; col<fDim1; col++) {
          this->s(col,row) = std::conj(this->s(row,col));
        }
      }else{
        DebugStop();//unreachable
      }
    }else{
      for(col=row+1; col<fDim1; col++) {
        this->s(col,row) = this->s(row,col);
      }
    }
    
  }
  
}



#ifndef USING_LAPACK    
#define NON_LAPACK \
    PZError<<__PRETTY_FUNCTION__<<" requires Lapack\n";             \
    PZError<<" Set either USING_LAPACK=ON or USING_MKL=ON on CMake ";   \
    PZError<<" when configuring NeoPZ library"<<std::endl;              \
    DebugStop();\
    return -1;

template<class TVar>
int TPZFMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &w, TPZFMatrix < CTVar > &eigenVectors){NON_LAPACK}

template<class TVar>
int TPZFMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &w){NON_LAPACK}

template<class TVar>
int TPZFMatrix<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix< TVar > &B , TPZVec < CTVar > &w, TPZFMatrix < CTVar > &eigenVectors){NON_LAPACK}

template<class TVar>
int TPZFMatrix<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix< TVar > &B , TPZVec < CTVar > &w){NON_LAPACK}        

template<typename TVar>
int TPZFMatrix<TVar>::SingularValueDecomposition(TPZFMatrix<TVar>& U, TPZFMatrix<TVar>& S, TPZFMatrix<TVar>& VT,char jobU, char jobVT){NON_LAPACK}
#undef NON_LAPACK
#endif

#include <complex>

template class TPZFMatrix<int >;
template class TPZFMatrix<int64_t >;
template class TPZFMatrix<float >;
template class TPZFMatrix<double >;
template class TPZFMatrix<long double>;

template class TPZFMatrix< std::complex<float> >;
template class TPZFMatrix< std::complex<double> >;
template class TPZFMatrix< std::complex<long double> >;

template class TPZFMatrix<TPZFlopCounter>;

template class TPZRestoreClass< TPZFMatrix<int> >;
template class TPZRestoreClass< TPZFMatrix<int64_t> >;
template class TPZRestoreClass< TPZFMatrix<double> >;
template class TPZRestoreClass< TPZFMatrix<float> >;
template class TPZRestoreClass< TPZFMatrix<long double> >;

template class TPZRestoreClass< TPZFMatrix<std::complex<float> > >;
template class TPZRestoreClass< TPZFMatrix<std::complex<double> > >;
template class TPZRestoreClass< TPZFMatrix<std::complex<long double> > >;
template class TPZRestoreClass< TPZFMatrix<TPZFlopCounter > >;

template class TPZFMatrix<TFad<6,REAL> >;
template class TPZFMatrix<Fad<double> >;
template class TPZFMatrix<Fad<float> >;
template class TPZFMatrix<Fad<long double> >;

template class TPZRestoreClass<TPZFMatrix<TFad<6,REAL> >>;
template class TPZRestoreClass<TPZFMatrix<Fad<double> >>;
template class TPZRestoreClass<TPZFMatrix<Fad<float> >>;
template class TPZRestoreClass<TPZFMatrix<Fad<long double> >>;

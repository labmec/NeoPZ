/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrix methods.
 */

#include <memory.h>

#include "pzsysmp.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "TPZParallelUtils.h"

// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix() : TPZRegisterClassId(&TPZSYsmpMatrix::ClassId),
TPZMatrix<TVar>() {
#ifdef USING_MKL
    fPardisoControl.SetMatrix(this);
#endif
    
#ifdef CONSTRUCTOR
    cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}
template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix(const TPZSYsmpMatrix<TVar> &cp) : 
    TPZRegisterClassId(&TPZSYsmpMatrix::ClassId),
    TPZMatrix<TVar>(cp), fIA(cp.fIA), fJA(cp.fJA), fA(cp.fA), fDiag(cp.fDiag)
#ifdef USING_MKL
    , fPardisoControl(cp.fPardisoControl)
#endif
    {
#ifdef USING_MKL
        fPardisoControl.SetMatrix(this);
#endif
    }
template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix(const int64_t rows,const int64_t cols ) : TPZRegisterClassId(&TPZSYsmpMatrix::ClassId),
TPZMatrix<TVar>(rows,cols) {

#ifdef USING_MKL
    fPardisoControl.SetMatrix(this);
#endif

#ifdef CONSTRUCTOR
	cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZSYsmpMatrix<TVar>::~TPZSYsmpMatrix() {
	// Deletes everything associated with a TPZSYsmpMatrix
#ifdef DESTRUCTOR
	cerr << "~TPZSYsmpMatrix()\n";
#endif
}

template<class TVar>
TPZSYsmpMatrix<TVar> &TPZSYsmpMatrix<TVar>::operator=(const TPZSYsmpMatrix<TVar> &copy) 
{
    TPZMatrix<TVar>::operator=(copy);
    fIA =copy.fIA;
    fJA = copy.fJA;
    fA = copy.fA;
    fDiag = copy.fDiag;
#ifdef USING_MKL
    fPardisoControl = copy.fPardisoControl;
    fPardisoControl.SetMatrix(this);
#endif
    return *this;
}

template <class TVar> int TPZSYsmpMatrix<TVar>::Zero() {
  fA.Fill(0.);
  fDiag.Fill(0.);
#ifndef USING_MKL
  TPZMatrix<TVar>::fDecomposed = ENoDecompose;
#endif
  return 0;
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

template<class TVar>
const TVar &TPZSYsmpMatrix<TVar>::GetVal(const int64_t r,const int64_t c ) const {
	// Get the matrix entry at (row,col) without bound checking
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }
	for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
		if ( fJA[ic] == col ) return fA[ic];
	}
	return this->gZero;
}

/** @brief Put values without bounds checking \n
 *  This method is faster than "Put" if DEBUG is defined.
 */
template<class TVar>
int TPZSYsmpMatrix<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & val )
{
    // Get the matrix entry at (row,col) without bound checking
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }
    for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col )
        {
            fA[ic] = val;
            return 0;
        }
    }
    if (val != (TVar(0.))) {
        DebugStop();
    }
    return 0;
    
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::AddKelAtomic(TPZFMatrix<TVar>&elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex){
    
    int64_t nelem = sourceindex.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    double prevval;
    int64_t row, col;
    
    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destinationindex[icoef];
        ieqs = sourceindex[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destinationindex[jcoef];
            jeqs = sourceindex[jcoef];
            {
                row = ieq; col = jeq;
                if (row > col) {
                    int64_t temp = row;
                    row = col;
                    col = temp;
                }
                int ic;
                for(ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
                    if ( fJA[ic] == col )
                    {
                        pzutils::AtomicAdd(fA[ic],elmat(ieqs,jeqs));
                        break;
                    }
                }
                if (ic == fIA[row+1]) {
                    if (elmat(ieqs,jeqs) != (TVar(0.))){
                        DebugStop();
                        
                    }
                    
                }
            }

        }
    
    }
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
							 TPZFMatrix<TVar> &z,
							 const TVar alpha,const TVar beta,const int opt) const {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	int64_t  ir, ic;
	int64_t  r = (opt) ? this->Cols() : this->Rows();
	
	// Determine how to initialize z
	if(IsZero(beta)) {
        z = y*beta;
	} else {
        z.Zero();
	}
	
	// Compute alpha * A * x
    int64_t ncols = x.Cols();
    for (int64_t col=0; col<ncols; col++)
    {
        for(int64_t ir=0; ir<this->Rows(); ir++) {
            for(int64_t ic=fIA[ir]; ic<fIA[ir+1]; ic++) {
                int64_t jc = fJA[ic];
                z(ir,col) += alpha * fA[ic] * x.g(jc,col);
                if(jc != ir)
                {
                    z(jc,col) += alpha * fA[ic] * x.g(ir,col);
                }
            }
        }
    }
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrix<TVar>::Print(const char *title, std::ostream &out ,const MatrixOutputFormat form) const {
	// Print the matrix along with a identification title
	if(form == EInputFormat) {
		out << "\nTSYsmpMatrix Print: " << title << '\n'
        << "\tNon zero elements    = " << fA.size()  << '\n'
		<< "\tRows    = " << this->Rows()  << '\n'
		<< "\tColumns = " << this->Cols() << '\n';
		int i;
		out << "\tIA\tJA\tA\n"
		<< "\t--\t--\t-\n";
		for(i=0; i<=this->Rows(); i++) {
			out << i      << '\t'
			<< fIA[i] << '\t'
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
		for(i=this->Rows()+1; i<fIA[this->Rows()]-1; i++) {
			out << i      << "\t\t"
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
	} else {
		TPZMatrix<TVar>::Print(title,out,form);
	}
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrix<TVar>::ComputeDiagonal() {
	if(!fDiag.size()) fDiag.resize(this->Rows());
	for(int ir=0; ir<this->Rows(); ir++) {
		fDiag[ir] = GetVal(ir,ir);
	}
}

/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZSYsmpMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric)
{
    if (!symmetric || nrow != ncol) {
        DebugStop();
    }
    TPZFMatrix<TVar> orig;
    orig.AutoFill(nrow,ncol,symmetric);
    
    TPZVec<int64_t> IA(nrow+1);
    TPZStack<int64_t> JA;
    TPZStack<TVar> A;
    IA[0] = 0;
    TPZVec<std::set<int64_t> > eqs(nrow);
    for (int64_t row=0; row<nrow; row++) {
        eqs[row].insert(row);
        for (int64_t col = 0; col<ncol; col++) {
            REAL test = rand()*1./RAND_MAX;
            if (test > 0.5) {
                eqs[row].insert(col);
                if (symmetric) {
                    eqs[col].insert(row);
                }
            }
        }
    }
    int64_t pos=0;
    for (int64_t row=0; row< nrow; row++) {
        for (std::set<int64_t>::iterator col = eqs[row].begin(); col != eqs[row].end(); col++) {
            if(*col >= row)
            {
                JA.Push(*col);
                A.Push(orig(row,*col));
            }
        }
        IA[row+1] = JA.size();
    }
    TPZMatrix<TVar>::Resize(nrow,ncol);
    SetData(IA, JA, A);
}

#ifdef USING_MKL

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
    Decompose_LDLt();
    return 1;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt()
{
    if(this->IsDecomposed() == ELDLt) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    fPardisoControl.SetMatrixType(TPZPardisoControl<TVar>::ESymmetric,TPZPardisoControl<TVar>::EIndefinite);
    fPardisoControl.Decompose();
    this->SetIsDecomposed(ELDLt);
    return 1;
    
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky()
{
    if(this->IsDecomposed() == ECholesky) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }

    fPardisoControl.SetMatrixType(TPZPardisoControl<TVar>::ESymmetric,TPZPardisoControl<TVar>::EPositiveDefinite);
    fPardisoControl.Decompose();

    this->SetIsDecomposed(ECholesky);
    return 1;
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    return Decompose_Cholesky();
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar>* b ) const
{
    TPZFMatrix<TVar> x(*b);
    fPardisoControl.Solve(*b,x);
    *b = x;
    return 1;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar>* b ) const
{
    return 1;
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>* b ) const
{
    TPZFMatrix<TVar> x(*b);
    fPardisoControl.Solve(*b,x);
    *b = x;
    return 1;
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

#else
//perhaps we could default to a less eficient implementation for
//solving. Perhaps removing the DebugStop() on PutVal() would be enough?
#define NOMKL \
    PZError<<__PRETTY_FUNCTION__<<" is not available if NeoPZ ";\
    PZError<<"was not configured with MKL. Aborting..."<<std::endl;\
    DebugStop();\
    return -1;


template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{    
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt()
{
    NOMKL
    
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky()
{
    NOMKL
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}
#endif


template<class TVar>
int TPZSYsmpMatrix<TVar>::ClassId() const{
    return Hash("TPZSYsmpMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template class TPZSYsmpMatrix<double>;
template class TPZSYsmpMatrix<float>;
template class TPZSYsmpMatrix<long double>;
template class TPZSYsmpMatrix<std::complex<float>>;
template class TPZSYsmpMatrix<std::complex<double>>;
template class TPZSYsmpMatrix<std::complex<long double>>;
